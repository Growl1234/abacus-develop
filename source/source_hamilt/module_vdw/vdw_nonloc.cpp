#include "vdw_nonloc.h"

#include "source_base/constants.h"
#include "source_base/module_parallel/para_mpi_func.h"
#include "source_base/timer.h"
#include "source_base/tool_quit.h"
#include "source_basis/module_pw/pw_basis.h"
#include "source_hamilt/module_xc/xc_functional.h"

#include <vdwxc.h>
#if defined(__MPI) && defined(VDWXC_ENABLE_MPI)
#include <vdwxc_mpi.h>
#endif

#include <algorithm>
#include <cassert>
#include <limits>
#include <numeric>

namespace vdw
{
// With lattice axes (c,a,b), both codes own physical z slabs. A block is the
// intersection with one peer. Host storage is [field][xy][z], library and
// packed storage are [field][z][xy].
struct Nonloc::Slabs
{
    struct Block
    {
        int peer;
        int start;
        int peer_start;
        int nz;
    };
    struct Side
    {
        std::vector<Block> blocks;
        std::vector<int> count;
        std::vector<int> disp;
        std::vector<double> buffer;

        void plan(
            int start,
            int nz,
            const int *lengths,
            const Parallel::ParaWorld &pool,
            size_t plane_size)
        {
            size_t elements = 0;
            int first = 0;
            for (int p = 0; p < pool.size(); ++p)
            {
                const int end = first + lengths[p];
                const int lo = std::max(start, first);
                const int n = std::max(0, std::min(start + nz, end) - lo);
                if (n)
                {
                    blocks.push_back({p, lo - start, lo - first, n});
                }
                const size_t msg = p == pool.rank() ? 0 : n * plane_size;
                if (elements + msg > std::numeric_limits<int>::max())
                {
                    ModuleBase::WARNING_QUIT("vdw::Nonloc", "Slab exchange exceeds the MPI int count limit");
                }
                disp.push_back(static_cast<int>(elements));
                count.push_back(static_cast<int>(msg));
                elements += msg;
                first = end;
            }
            buffer.resize(std::max<size_t>(1, elements));
        }
    };

    const Parallel::ParaWorld pool;
    const ModulePW::PW_Basis &pw;
    const int fields;
    const size_t stride;
    std::vector<double> input;
    Side host;
    Side lib;

    Slabs(
        const Parallel::ParaWorld &pool,
        const ModulePW::PW_Basis &pw,
        int fields,
        int z0,
        int nz)
        : pool(pool), pw(pw), fields(fields), stride(std::max<size_t>(1, static_cast<size_t>(nz) * pw.nxy)),
          input(fields * stride)
    {
        std::vector<int> lengths(pool.size());
        Parallel::allgather_int(nz, lengths.data(), pool);
        if (std::accumulate(lengths.begin(), lengths.end(), 0) != pw.nz ||
            (nz && z0 != std::accumulate(lengths.begin(), lengths.begin() + pool.rank(), 0)))
        {
            ModuleBase::WARNING_QUIT("vdw::Nonloc", "Expected libvdwxc z slabs in rank order");
        }
        if (std::equal(lengths.begin(), lengths.end(), pw.numz))
        {
            if (nz)
            {
                host.blocks.push_back({pool.rank(), 0, 0, nz});
            }
            return;
        }
        const size_t plane_size = static_cast<size_t>(fields) * pw.nxy;
        host.plan(pw.startz_current, pw.nplane, lengths.data(), pool, plane_size);
        lib.plan(z0, nz, pw.numz, pool, plane_size);
    }

    // Every block writes a disjoint set of planes. A thread team can handle all
    // local and outgoing blocks without a barrier between them.
    template <bool forward>
    void transpose(double *host_fields, double *lib_fields)
    {
        const size_t xy = pw.nxy;
#pragma omp parallel
        {
            for (const Block &b : host.blocks)
            {
                const bool local = b.peer == pool.rank();
                double *slab = local ? lib_fields + b.peer_start * xy : host.buffer.data() + host.disp[b.peer];
                const size_t field_stride = local ? stride : b.nz * xy;
#pragma omp for collapse(2) schedule(static) nowait
                for (int f = 0; f < fields; ++f)
                {
                    for (int i = 0; i < pw.nxy; ++i)
                    {
                        double *row = host_fields + (f * xy + i) * pw.nplane + b.start;
                        double *col = slab + f * field_stride + i;
                        for (int z = 0; z < b.nz; ++z)
                        {
                            if (forward)
                            {
                                col[z * xy] = row[z];
                            }
                            else
                            {
                                row[z] = col[z * xy];
                            }
                        }
                    }
                }
            }
        }
    }

    template <bool forward>
    void copy_remote(double *lib_fields)
    {
        const size_t xy = pw.nxy;
        for (const Block &b : lib.blocks)
        {
            if (b.peer == pool.rank())
            {
                continue;
            }
            for (int f = 0; f < fields; ++f)
            {
                double *plane = lib_fields + f * stride + b.start * xy;
                double *packed = lib.buffer.data() + lib.disp[b.peer] + f * xy * b.nz;
                std::copy_n(forward ? packed : plane, b.nz * xy, forward ? plane : packed);
            }
        }
    }

    template <bool forward>
    void transfer(double *host_fields, double *lib_fields)
    {
        ModuleBase::timer::Scope timer("Nonloc", "slab_transfer");
        if (forward)
        {
            transpose<true>(host_fields, lib_fields);
        }
        if (!host.count.empty())
        {
            if (!forward)
            {
                copy_remote<false>(lib_fields);
            }
            const Side &send = forward ? host : lib;
            Side &recv = forward ? lib : host;
            Parallel::alltoallv_double(
                send.buffer.data(), send.count.data(), send.disp.data(), recv.buffer.data(), recv.count.data(),
                recv.disp.data(), pool);
            if (forward)
            {
                copy_remote<true>(lib_fields);
            }
        }
        if (!forward)
        {
            transpose<false>(host_fields, lib_fields);
        }
    }
};

Nonloc::Nonloc(
    const std::string &method,
    double b,
    double c,
    int nspin,
    const ModulePW::PW_Basis &pw)
    : basis_(pw), channels_(method == "rvv10" ? 1 : nspin), total_spin_(nspin == 2 && channels_ == 1)
{
    const int kernel = method == "rvv10"      ? VDWXC_RVV10
                       : method == "vdw_dfcx" ? VDWXC_DFCX
                       : method == "vdw_df2"  ? VDWXC_DF2
                                              : VDWXC_DF1;
    assert(is_nonloc(method) && (nspin == 1 || nspin == 2));
    handle_ = channels_ == 2 ? vdwxc_new_spin(kernel) : vdwxc_new(kernel);
    if (kernel == VDWXC_RVV10 && vdwxc_set_rvv10_parameters(handle_, b, c) != 0)
    {
        ModuleBase::WARNING_QUIT("vdw::Nonloc", "rVV10 requires finite b > 0 and C >= 0");
    }
    const bool distributed = pw.poolnproc > 1;
    const int n[] = {pw.nx, pw.ny, pw.nz};
    // libvdwxc uses a signed determinant for dV. A Cartesian inversion preserves all distances.
    const ModuleBase::Matrix3 a = (pw.latvec.Det() < 0 ? -pw.lat0 : pw.lat0) * pw.latvec;
    const double v[][3] = {{a.e11, a.e12, a.e13}, {a.e21, a.e22, a.e23}, {a.e31, a.e32, a.e33}};
    const int i = distributed ? 2 : 0;
    const int j = (i + 1) % 3;
    const int k = (j + 1) % 3;
    vdwxc_set_unit_cell(
        handle_, n[i], n[j], n[k], v[i][0], v[i][1], v[i][2], v[j][0], v[j][1], v[j][2], v[k][0], v[k][1], v[k][2]);
    const int nf = 2 * channels_;
    if (distributed)
    {
#if defined(__MPI) && defined(VDWXC_ENABLE_MPI)
        vdwxc_init_mpi(handle_, pw.pool_world);
        int shape[3];
        int offset[3];
        vdwxc_get_local_grid(handle_, shape, offset);
        if (shape[1] != pw.nx || shape[2] != pw.ny || offset[1] || offset[2])
        {
            ModuleBase::WARNING_QUIT("vdw::Nonloc", "Expected a libvdwxc slab spanning the xy plane");
        }
        slabs_.reset(new Slabs(Parallel::ParaWorld::make_mpi("pw", pw.pool_world), pw, nf, offset[0], shape[0]));
#else
        ModuleBase::WARNING_QUIT("vdw::Nonloc", "Multi-rank pools require an MPI-enabled libvdwxc");
#endif
    }
    else
    {
        vdwxc_init_serial(handle_);
    }
    const size_t stride = std::max(1, pw.nrxx);
    fields_.resize(nf * stride);
    deriv_.resize(slabs_ ? slabs_->input.size() : fields_.size());
    gradient_.resize(channels_ * stride);
    rhog_.resize(std::max(1, pw.npw));
}

Nonloc::~Nonloc() { vdwxc_finalize(&handle_); }

void Nonloc::request_stress(bool needed)
{
    const int n = needed ? 3 : 0;
    stress_.create(n, n, false);
}

const ModuleBase::matrix &Nonloc::stress() const
{
    if (!stress_.c)
    {
        ModuleBase::WARNING_QUIT("vdw::Nonloc", "Stress was not requested for the current XC potential");
    }
    return stress_;
}

void Nonloc::add_xc(
    const double *const *valence,
    const double *core,
    double &etxc,
    double &vtxc,
    ModuleBase::matrix &v)
{
    ModuleBase::timer::Scope timer("Nonloc", "add_xc");
    const int nr = basis_.nrxx;
    const size_t stride = std::max(1, nr);
    const double core_fraction = 1.0 / channels_;
    for (int s = 0; s < channels_; ++s)
    {
        double *rho = fields_.data() + s * stride;
        double *sigma = fields_.data() + (channels_ + s) * stride;
        auto *grad = gradient_.data() + s * stride;
#pragma omp parallel for
        for (int r = 0; r < nr; ++r)
        {
            rho[r] = valence[s][r] + core_fraction * core[r];
            if (total_spin_)
            {
                rho[r] += valence[1][r];
            }
        }
        basis_.real2recip(rho, rhog_.data());
        XC_Functional::grad_rho(rhog_.data(), grad, &basis_, basis_.tpiba);
#pragma omp parallel for
        for (int r = 0; r < nr; ++r)
        {
            sigma[r] = grad[r].norm2();
        }
    }

    double *in = fields_.data();
    double *out = deriv_.data();
    if (slabs_)
    {
        in = slabs_->input.data();
        slabs_->transfer<true>(fields_.data(), in);
    }
    const size_t m = slabs_ ? slabs_->stride : stride;
    double *tensor = stress_.c;
    double e;
    if (channels_ == 1)
    {
        e = tensor ? vdwxc_stress(handle_, tensor, in, in + m, out, out + m)
                   : vdwxc_calculate(handle_, in, in + m, out, out + m);
    }
    else
    {
        e = tensor ? vdwxc_stress_spin(
                         handle_, tensor, in, in + m, in + 2 * m, in + 3 * m, out, out + m, out + 2 * m, out + 3 * m)
                   : vdwxc_calculate_spin(
                         handle_, in, in + m, in + 2 * m, in + 3 * m, out, out + m, out + 2 * m, out + 3 * m);
    }
    if (slabs_)
    {
        slabs_->transfer<false>(fields_.data(), out);
        out = fields_.data();
    }

    double nv = 0;
    double xx = 0;
    double xy = 0;
    double xz = 0;
    double yy = 0;
    double yz = 0;
    double zz = 0;
    for (int s = 0; s < channels_; ++s)
    {
        auto *grad = gradient_.data() + s * stride;
        const double *dn = out + s * stride;
        double *dsigma = out + (channels_ + s) * stride;
        if (tensor)
        {
            // Strain changes the Cartesian density gradient. Contract it with the same flux used for Vxc.
#pragma omp parallel for reduction(+ : xx, xy, xz, yy, yz, zz)
            for (int r = 0; r < nr; ++r)
            {
                const auto g = grad[r];
                grad[r] *= 2.0 * dsigma[r];
                const auto &flux = grad[r];
                xx += flux.x * g.x;
                xy += flux.x * g.y;
                xz += flux.x * g.z;
                yy += flux.y * g.y;
                yz += flux.y * g.z;
                zz += flux.z * g.z;
            }
        }
        else
        {
#pragma omp parallel for
            for (int r = 0; r < nr; ++r)
            {
                grad[r] *= 2.0 * dsigma[r];
            }
        }
        // Once the flux is formed, reuse its sigma-derivative buffer for the divergence.
        XC_Functional::grad_dot(grad, dsigma, &basis_, basis_.tpiba);
#pragma omp parallel for reduction(+ : nv)
        for (int r = 0; r < nr; ++r)
        {
            const double potential = ModuleBase::e2 * (dn[r] - dsigma[r]);
            v(s, r) += potential;
            nv += valence[s][r] * potential;
            if (total_spin_)
            {
                v(1, r) += potential;
                nv += valence[1][r] * potential;
            }
        }
    }
    // libvdwxc integrates its rank-local energy in Ha; band-energy subtraction uses only valence charge.
    double terms[] = {ModuleBase::e2 * e, nv * basis_.omega / basis_.nxyz, xx, xy, xz, yy, yz, zz};
    if (tensor)
    {
        // The host XC stress already contains -E/volume. Remove that volume derivative from the library tensor.
        int k = 2;
        for (int i = 0; i < 3; ++i)
        {
            for (int j = i; j < 3; ++j, ++k)
            {
                const double kernel = tensor[3 * i + j] - (i == j ? e : 0.0);
                terms[k] = ModuleBase::e2 * (terms[k] / basis_.nxyz - kernel / basis_.omega);
            }
        }
    }
    if (slabs_)
    {
        Parallel::reduce_all(terms, tensor ? 8 : 2, slabs_->pool);
    }
    if (tensor)
    {
        int k = 2;
        for (int i = 0; i < 3; ++i)
        {
            for (int j = i; j < 3; ++j)
            {
                stress_(i, j) = stress_(j, i) = terms[k++];
            }
        }
    }
    etxc += terms[0];
    vtxc += terms[1];
}

void Nonloc::add_core_potential(double *vcore) const
{
    // add_xc leaves dedn and the flux divergence in the host layout.
    const double *out = slabs_ ? fields_.data() : deriv_.data();
    const size_t stride = std::max(1, basis_.nrxx);
    const double weight = ModuleBase::e2 / channels_;
    for (int s = 0; s < channels_; ++s)
    {
        const double *dn = out + s * stride;
        const double *div = out + (channels_ + s) * stride;
#pragma omp parallel for
        for (int r = 0; r < basis_.nrxx; ++r)
        {
            vcore[r] += weight * (dn[r] - div[r]);
        }
    }
}
} // namespace vdw
