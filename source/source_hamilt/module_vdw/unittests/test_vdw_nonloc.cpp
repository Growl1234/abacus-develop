#include "../vdw_nonloc.h"

#include "source_base/constants.h"
#include "source_base/matrix.h"
#include "source_base/module_parallel/para_mpi_func.h"
#include "source_basis/module_pw/pw_basis.h"

#include <algorithm>
#include <cmath>
#include <gtest/gtest.h>

namespace
{
const ModuleBase::Matrix3 cell(
    1,
    .21,
    -.08,
    .13,
    1.08,
    .23,
    -.11,
    .14,
    .92);

struct Grid
{
    ModulePW::PW_Basis pw;
    Grid(MPI_Comm comm, const ModuleBase::Matrix3 &lattice)
    {
        int size;
        int rank;
        MPI_Comm_size(comm, &size);
        MPI_Comm_rank(comm, &rank);
        pw.initmpi(size, rank, comm);
        // Three ranks split z as (4,3,3) in ABACUS and (4,4,2) in libvdwxc.
        pw.initgrids(8.0, lattice, 12, 14, 10);
        pw.initparameters(false, 8.0);
        pw.setuptransform();
        pw.collect_local_pw();
    }
};

struct Result
{
    double energy = 0;
    double nv = 0;
    double variation = 0;
    double core_force = 0;
    ModuleBase::matrix v;
    Result(int ns, int nr) : v(ns, nr) {}
};

Result evaluate(
    vdw::Nonloc &nonloc,
    const ModulePW::PW_Basis &pw,
    int ns,
    double step,
    double core_scale,
    double core_shift)
{
    Result result(ns, pw.nrxx);
    ModuleBase::matrix rho(ns, pw.nrxx);
    ModuleBase::matrix direction(ns, pw.nrxx);
    std::vector<double> core(pw.nrxx);
    std::vector<double> dcore(pw.nrxx);
    for (int r = 0; r < pw.nrxx; ++r)
    {
        const double x = ModuleBase::TWO_PI * (r / pw.nplane / pw.ny) / pw.nx;
        const double y = ModuleBase::TWO_PI * (r / pw.nplane % pw.ny) / pw.ny;
        const double z = ModuleBase::TWO_PI * (r % pw.nplane + pw.startz_current) / pw.nz;
        const double c = .012 + .002 * std::cos(z);
        core[r] = core_scale * (.012 + .002 * std::cos(z - core_shift));
        dcore[r] = core_scale * .002 * std::sin(z - core_shift);
        for (int s = 0; s < ns; ++s)
        {
            // The mixed x-Nyquist modes enter the nonlocal convolution even though the
            // ABACUS gradient cutoff projects them out. Core transfer holds the XC density fixed.
            const double n =
                s == 0 ? .055 + .008 * std::sin(x + y) + .004 * std::cos(y + z) + .003 * std::cos(6 * x + y + z)
                       : .038 + .006 * std::cos(x - z) + .003 * std::sin(y + z) + .002 * std::cos(6 * x - y);
            direction(s, r) = s == 0 ? .011 * std::cos(y + z) + .004 * std::sin(x + y) + .001 * std::cos(6 * x + y + z)
                                     : .008 * std::sin(y + z) + .006 * std::cos(x - z);
            direction(s, r) += .006 * std::sin(z);
            rho(s, r) = n + (1.0 - core_scale) * c / ns + step * direction(s, r);
        }
    }
    const double *rows[2] = {rho.c, ns == 2 ? rho.c + pw.nrxx : nullptr};
    nonloc.add_xc(rows, core.data(), result.energy, result.nv, result.v);
    std::vector<double> vcore(pw.nrxx);
    nonloc.add_core_potential(vcore.data());
    double moments[3] = {};
    for (int r = 0; r < pw.nrxx; ++r)
    {
        moments[2] -= vcore[r] * dcore[r];
    }
    for (int i = 0; i < ns * pw.nrxx; ++i)
    {
        moments[0] += result.v.c[i] * rho.c[i];
        moments[1] += result.v.c[i] * direction.c[i];
    }
    const Parallel::ParaWorld pool = Parallel::ParaWorld::make_mpi("pw", pw.pool_world);
    Parallel::reduce_all(moments, 3, pool);
    EXPECT_NEAR(result.nv, moments[0] * pw.omega / pw.nxyz, 2e-13);
    result.variation = moments[1] * pw.omega / pw.nxyz;
    result.core_force = moments[2] * pw.omega / pw.nxyz;
    return result;
}

TEST(Nonloc, EnergyPotentialAndCoreCharge)
{
    Grid distributed(MPI_COMM_WORLD, cell);
    Grid serial(MPI_COMM_SELF, cell);
    const auto &pw = distributed.pw;
    const char *methods[] = {"vdw_df1", "vdw_df2", "vdw_dfcx", "rvv10"};
    // Independent serial C API energies in Ry from the analytic gradients of the retained modes.
    const double ref[][2] = {{-.00139600507024414, -.00424528438696795}, {-.00139379282377613, -.00424340065568855},
        {-.00139600507024414, -.00424528438696795}, {.000543203615197579, .000615238510137339}};
    for (int k = 0; k < 4; ++k)
    {
        for (int ns = 1; ns <= 2; ++ns)
        {
            SCOPED_TRACE(::testing::Message() << methods[k] << " nspin=" << ns);
            vdw::Nonloc adapter(methods[k], 8.5, .011, ns, pw);
            vdw::Nonloc serial_adapter(methods[k], 8.5, .011, ns, serial.pw);
            adapter.request_stress(true);
            serial_adapter.request_stress(true);
            const Result initial = evaluate(adapter, pw, ns, 0, 1, 0);
            const Result serial_ref = evaluate(serial_adapter, serial.pw, ns, 0, 1, 0);
            EXPECT_NEAR(initial.energy, ref[k][ns - 1], 2e-12);
            EXPECT_NEAR(initial.nv, serial_ref.nv, 2e-12);

            for (int i = 0; i < 9; ++i)
            {
                EXPECT_NEAR(adapter.stress().c[i], serial_adapter.stress().c[i], 2e-13);
            }
            adapter.request_stress(false);
            const double step = 1e-4;
            const double plus = evaluate(adapter, pw, ns, step, 1, 0).energy;
            const double minus = evaluate(adapter, pw, ns, -step, 1, 0).energy;
            EXPECT_NEAR((plus - minus) / (2 * step), initial.variation, 2e-9);

            // Break z symmetry in the valence density, then translate only the core.
            const double offset = 0.25;
            const double force = evaluate(adapter, pw, ns, offset, 1, 0).core_force;
            const double core_plus = evaluate(adapter, pw, ns, offset, 1, step).energy;
            const double core_minus = evaluate(adapter, pw, ns, offset, 1, -step).energy;
            EXPECT_GT(std::abs(force), 2e-8);
            EXPECT_NEAR(force, (core_minus - core_plus) / (2 * step), 2e-9);

            const Result moved = evaluate(adapter, pw, ns, 0, 0.5, 0);
            EXPECT_NEAR(moved.energy, initial.energy, 2e-13);
            EXPECT_GT(std::abs(moved.nv - initial.nv), 1e-8);
            double error = 0;
            for (int s = 0; s < ns; ++s)
            {
                for (int r = 0; r < pw.nrxx; ++r)
                {
                    const int global = (r / pw.nplane) * pw.nz + pw.startz_current + r % pw.nplane;
                    error = std::max(error, std::abs(initial.v(s, r) - serial_ref.v(s, global)));
                    error = std::max(error, std::abs(initial.v(s, r) - moved.v(s, r)));
                    if (k == 3 && ns == 2)
                    {
                        error = std::max(error, std::abs(initial.v(s, r) - initial.v(0, r)));
                    }
                }
            }
            EXPECT_LT(error, 5e-12);
        }
    }
}

ModuleBase::Matrix3 strained_cell(int i, int j, double step)
{
    auto a = cell.to_matrix();
    for (int row = 0; row < 3; ++row)
    {
        a(row, i) += step * a(row, j);
    }
    return ModuleBase::Matrix3(a(0, 0), a(0, 1), a(0, 2), a(1, 0), a(1, 1), a(1, 2), a(2, 0), a(2, 1), a(2, 2));
}

TEST(Nonloc, StressFromStrainedCells)
{
    Grid central(MPI_COMM_WORLD, cell);
    const double h = 1e-5;
    for (int i = 0; i < 3; ++i)
    {
        for (int j = i; j < 3; ++j)
        {
            Grid plus(MPI_COMM_WORLD, strained_cell(i, j, h));
            Grid minus(MPI_COMM_WORLD, strained_cell(i, j, -h));
            for (const char *method : {"vdw_df1", "vdw_df2", "vdw_dfcx", "rvv10"})
            {
                for (int ns : {1, 2})
                {
                    SCOPED_TRACE(::testing::Message() << method << " nspin=" << ns << " strain=" << i << j);
                    vdw::Nonloc adapter(method, 8.5, .011, ns, central.pw);
                    adapter.request_stress(true);
                    const double e = evaluate(adapter, central.pw, ns, 0, 0, 0).energy;
                    vdw::Nonloc pos(method, 8.5, .011, ns, plus.pw);
                    vdw::Nonloc neg(method, 8.5, .011, ns, minus.pw);
                    const double ep = evaluate(pos, plus.pw, ns, 0, 0, 0).energy;
                    const double em = evaluate(neg, minus.pw, ns, 0, 0, 0).energy;
                    // With density amplitudes fixed, the only host term is the XC energy on the diagonal.
                    const double sigma = adapter.stress()(i, j) - (i == j ? e / central.pw.omega : 0.0);
                    EXPECT_NEAR(sigma, (em - ep) / (2 * h * central.pw.omega), 2e-12);
                    EXPECT_DOUBLE_EQ(adapter.stress()(i, j), adapter.stress()(j, i));
                }
            }
        }
    }
}
} // namespace

int main(int argc, char **argv)
{
    MPI_Init(&argc, &argv);
    testing::InitGoogleTest(&argc, argv);
    const int result = RUN_ALL_TESTS();
    MPI_Finalize();
    return result;
}
