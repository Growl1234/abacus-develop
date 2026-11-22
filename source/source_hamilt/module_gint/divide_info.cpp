#include "divide_info.h"

#include "source_base/constants.h"
#include "source_base/module_parallel/para_mpi_func.h"
#include "source_base/module_parallel/para_tag.h"
#include "source_base/parallel_common.h"
#include "source_base/timer.h"
#include "source_basis/module_pw/pw_basis_big.h"
#include "source_cell/unitcell.h"
#include <algorithm>
#include <cmath>
#include <complex>
#include <functional>
#include <limits>

namespace ModuleGint
{

namespace
{
int slab_width(int planes, int ranks, int rank) { return planes / ranks + (rank < planes % ranks); }

bool can_transfer(const ModulePW::PW_Basis_Big &grid, SlabAxis axis)
{
    const int ranks = grid.poolnproc;
    const bool along_a = axis == SlabAxis::a;
    const int blocks = along_a ? grid.nbx : grid.nby;
    if (axis == SlabAxis::c || ranks == 1 || blocks < ranks || grid.nbz < ranks)
    {
        return false;
    }
    // Alltoallv counts and offsets must fit in an int on every rank.
    const int width = slab_width(blocks, ranks, 0) * (along_a ? grid.bx : grid.by);
    const long long lines = static_cast<long long>(width) * (along_a ? grid.ny : grid.nx);
    const int most_z = *std::max_element(grid.numz, grid.numz + ranks);
    const int limit = std::numeric_limits<int>::max();
    return lines <= limit / grid.nz && static_cast<long long>(grid.nx) * grid.ny <= limit / most_z;
}

// Fraction of the volume of a unit sphere below a plane at signed distance z.
double sphere_fraction(double z)
{
    z = std::max(-1.0, std::min(1.0, z));
    return 0.5 + 0.25 * z * (3.0 - z * z);
}

struct SlabWork
{
    double total = 0;
    double peak = 0;
};

SlabWork projected_work(
    const UnitCell &cell,
    SlabAxis axis,
    int planes,
    int ranks)
{
    const int d = static_cast<int>(axis);
    const double metric[] = {cell.GGT.e11, cell.GGT.e22, cell.GGT.e33};
    const double planes_per_bohr = planes * std::sqrt(metric[d]) / cell.lat0;
    const double plane_volume = cell.omega / planes;
    std::vector<double> orbitals(planes, 0.0);
    for (int it = 0; it < cell.ntype; ++it)
    {
        const Atom &atom = cell.atoms[it];
        const double radius = atom.Rcut * planes_per_bohr;
        if (radius <= 0 || atom.nw == 0)
        {
            continue;
        }
        // Bin equivalent atoms first. Convolve each occupied plane with this species' cutoff only once.
        std::vector<int> centers(planes, 0);
        for (const auto &pos : atom.taud)
        {
            const double f = pos[d] - std::floor(pos[d]);
            const int index = std::min(static_cast<int>(f * planes), planes - 1);
            ++centers[index];
        }
        const int reach = static_cast<int>(std::floor(radius + 0.5));
        const int length = std::min(planes, 2 * reach + 1);
        std::vector<double> kernel(length, 0.0);
        const double weight = (ModuleBase::FOUR_PI / 3.0) * atom.nw * atom.Rcut * atom.Rcut * atom.Rcut / plane_volume;
        double left = 0.0;
        for (int k = -reach; k <= reach; ++k)
        {
            const double right = sphere_fraction((k + 0.5) / radius);
            kernel[(k + reach) % planes] += weight * (right - left);
            left = right;
        }
        for (int i = 0; i < planes; ++i)
        {
            if (centers[i] == 0)
            {
                continue;
            }
            int plane = (i + planes - reach % planes) % planes;
            for (double contribution : kernel)
            {
                orbitals[plane] += centers[i] * contribution;
                if (++plane == planes)
                {
                    plane = 0;
                }
            }
        }
    }
    // Gint's pair contractions scale roughly as mesh points times overlapping orbitals squared.
    SlabWork work;
    int plane = 0;
    for (int rank = 0; rank < ranks; ++rank)
    {
        double local = 0.0;
        const int end = plane + slab_width(planes, ranks, rank);
        for (; plane < end; ++plane)
        {
            local += orbitals[plane] * orbitals[plane];
        }
        work.peak = std::max(work.peak, local);
        work.total += local;
    }
    return work;
}

SlabAxis lowest_cost_axis(const ModulePW::PW_Basis_Big &grid, const UnitCell &cell)
{
    const int ranks = grid.poolnproc;
    const bool allowed[] = {can_transfer(grid, SlabAxis::a), can_transfer(grid, SlabAxis::b)};
    if (!allowed[0] && !allowed[1])
    {
        return SlabAxis::c;
    }
    const SlabWork native = projected_work(cell, SlabAxis::c, grid.nbz, ranks);
    // Leave a 10% margin on the mean rank work for redistribution and uncertainty in the projection.
    const double margin = 0.10 / ranks;
    if (native.total == 0 || native.peak <= native.total * (1.0 / ranks + margin))
    {
        return SlabAxis::c;
    }
    SlabAxis best = SlabAxis::c;
    double best_cost = native.peak / native.total;
    const double mean_pairs = native.total / grid.nbz;
    const double fft_share = static_cast<double>(slab_width(grid.nbz, ranks, 0)) / grid.nbz;
    for (SlabAxis axis : {SlabAxis::a, SlabAxis::b})
    {
        if (!allowed[static_cast<int>(axis)])
        {
            continue;
        }
        const int planes = axis == SlabAxis::a ? grid.nbx : grid.nby;
        const double gint_share = static_cast<double>(slab_width(planes, ranks, 0)) / planes;
        // Charge for moving the full scalar field even when most of the grid is vacuum.
        // A b slab also packs strided rows on the FFT side.
        const double traffic = (gint_share + (axis == SlabAxis::a ? fft_share : 2 * fft_share)) / mean_pairs;
        const double penalty = std::max(margin, traffic);
        if (1.0 / ranks + penalty < best_cost)
        {
            const SlabWork candidate = projected_work(cell, axis, planes, ranks);
            const double cost = candidate.peak / candidate.total + penalty;
            if (cost < best_cost)
            {
                best = axis;
                best_cost = cost;
            }
        }
    }
    return best;
}
} // namespace

SlabAxis DivideInfo::select_axis(const ModulePW::PW_Basis_Big &grid, const UnitCell &cell)
{
    if (grid.poolnproc == 1)
    {
        return SlabAxis::c;
    }
    ModuleBase::timer::start("DivideInfo", "select_axis");
    int axis = static_cast<int>(SlabAxis::c);
    if (grid.poolrank == 0)
    {
        axis = static_cast<int>(lowest_cost_axis(grid, cell));
    }
#ifdef __MPI
    const auto pool = Parallel::ParaWorld::make_mpi(Parallel::ParaTag::pw, grid.pool_world);
    Parallel::bcast_int(axis, pool);
#endif
    ModuleBase::timer::end("DivideInfo", "select_axis");
    return static_cast<SlabAxis>(axis);
}

void DivideInfo::Layout::append(Rows runs)
{
    const int count = runs.number * runs.length;
    contiguous = contiguous && runs.first == size && (runs.number == 1 || runs.stride == runs.length);
    rows.push_back(runs);
    offsets.push_back(size);
    counts.push_back(count);
    size += count;
}

DivideInfo::DivideInfo(const ModulePW::PW_Basis_Big &grid, SlabAxis axis) : grid_(grid)
{
    if (!can_transfer(grid, axis))
    {
        block_start_.z = grid.nbzp_start;
        block_count_.set(grid.nbx, grid.nby, grid.nbzp);
        return;
    }
    const int ranks = grid.poolnproc;

    const bool along_a = axis == SlabAxis::a;
    const int blocks = along_a ? grid.nbx : grid.nby;
    const int d = static_cast<int>(axis);
    const int block_width = along_a ? grid.bx : grid.by;
    const int per_rank = blocks / ranks;
    const int extra = blocks % ranks;
    block_start_[d] = grid.poolrank * per_rank + std::min(grid.poolrank, extra);
    block_count_.set(grid.nbx, grid.nby, grid.nbz);
    block_count_[d] = slab_width(blocks, ranks, grid.poolrank);
    const int local_lines = block_count_.x * grid.bx * block_count_.y * grid.by;

    int first = 0;
    for (int p = 0; p < ranks; ++p)
    {
        const int width = slab_width(blocks, ranks, p) * block_width;
        // An a slab is contiguous on the FFT side. A b slab has one run per x.
        if (along_a)
        {
            const int plane = grid.ny * grid.nplane;
            fft_.append({first * plane, 1, 0, width * plane});
        }
        else
        {
            fft_.append({first * grid.nplane, grid.nx, grid.ny * grid.nplane, width * grid.nplane});
        }
        // On the Gint side each local (x,y) line has the full z extent.
        gint_.append({grid.startz[p], local_lines, grid.nz, grid.numz[p]});
        first += width;
    }
}

template <typename T>
void DivideInfo::transfer(
    const T *src,
    const Layout &send,
    T *dst,
    const Layout &recv,
    bool add) const
{
#ifdef __MPI
    ModuleBase::timer::start("DivideInfo", "transfer");
    std::vector<T> packed;
    if (!send.contiguous)
    {
        packed.resize(send.size);
        for (std::size_t p = 0; p < send.rows.size(); ++p)
        {
            const Rows &r = send.rows[p];
            for (int i = 0; i < r.number; ++i)
            {
                std::copy_n(src + r.first + i * r.stride, r.length, packed.data() + send.offsets[p] + i * r.length);
            }
        }
        src = packed.data();
    }

    std::vector<T> received(recv.size);
    Parallel_Common::alltoallv(
        src, send.counts.data(), send.offsets.data(), received.data(), recv.counts.data(), recv.offsets.data(),
        grid_.pool_world);
    for (std::size_t p = 0; p < recv.rows.size(); ++p)
    {
        const Rows &r = recv.rows[p];
        for (int i = 0; i < r.number; ++i)
        {
            const T *in = received.data() + recv.offsets[p] + i * r.length;
            T *out = dst + r.first + i * r.stride;
            if (add)
            {
                std::transform(in, in + r.length, out, out, std::plus<T>());
            }
            else
            {
                std::copy_n(in, r.length, out);
            }
        }
    }
    ModuleBase::timer::end("DivideInfo", "transfer");
#endif
}

const double *DivideInfo::from_fft(const double *field, std::vector<double> &scratch) const
{
    if (!needs_transfer())
    {
        return field;
    }
    scratch.resize(gint_.size);
    transfer(field, fft_, scratch.data(), gint_, false);
    return scratch.data();
}

int DivideInfo::fft_size() const { return grid_.nrxx; }

template <typename T>
void DivideInfo::copy_to_fft(const T *local, T *field) const
{
    if (needs_transfer())
    {
        transfer(local, gint_, field, fft_, false);
    }
    else if (grid_.nrxx > 0 && local != field)
    {
        std::copy_n(local, grid_.nrxx, field);
    }
}

template void DivideInfo::copy_to_fft(const double *, double *) const;
template void DivideInfo::copy_to_fft(const std::complex<double> *, std::complex<double> *) const;

void DivideInfo::add_to_fft(const double *contribution, double *field) const
{
    if (needs_transfer())
    {
        transfer(contribution, gint_, field, fft_, true);
    }
    else if (grid_.nrxx > 0)
    {
        std::transform(contribution, contribution + grid_.nrxx, field, field, std::plus<double>());
    }
}

} // namespace ModuleGint
