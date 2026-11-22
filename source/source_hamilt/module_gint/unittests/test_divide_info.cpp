#include "../divide_info.h"
#include "source_base/parallel_comm.h"
#include "source_basis/module_pw/pw_basis_big.h"
#include "source_cell/unitcell.h"

#include <algorithm>
#include <complex>
#include <gtest/gtest.h>
#include <limits>

namespace
{
using ModuleGint::DivideInfo;
using ModuleGint::SlabAxis;

void check_fields(ModulePW::PW_Basis_Big &grid, int pool_id)
{
    grid.setbxyz(2, 3, 2);
    grid.nbx = 5;
    grid.nby = 7;
    grid.nbz = 11;
    grid.nx = grid.nbx * grid.bx;
    grid.ny = grid.nby * grid.by;
    grid.nz = grid.nbz * grid.bz;
    grid.nxy = grid.nx * grid.ny;
    grid.distribute_r();

    // A distinct global point ID in each pool detects misplaced rows and cross-pool traffic.
    const auto id = [&](int x, int y, int z) { return ((pool_id * grid.nx + x) * grid.ny + y) * grid.nz + z; };
    std::vector<double> potential(grid.nrxx);
    for (int x = 0; x < grid.nx; ++x)
    {
        for (int y = 0; y < grid.ny; ++y)
        {
            for (int z = 0; z < grid.nplane; ++z)
            {
                potential[(x * grid.ny + y) * grid.nplane + z] = id(x, y, z + grid.startz_current);
            }
        }
    }

    for (SlabAxis axis : {SlabAxis::a, SlabAxis::b, SlabAxis::c})
    {
        SCOPED_TRACE(static_cast<int>(axis));
        const DivideInfo layout(grid, axis);
        const int d = static_cast<int>(axis);
        const int blocks[] = {grid.nbx, grid.nby, grid.nbz};
        const bool exchange =
            grid.poolnproc > 1 && axis != SlabAxis::c && blocks[d] >= grid.poolnproc && grid.nbz >= grid.poolnproc;
        const int split = exchange ? d : 2;
        const auto &count = layout.block_count();
        const auto &start = layout.block_start();
        for (int i = 0; i < 3; ++i)
        {
            const int base = blocks[i] / grid.poolnproc;
            const int rem = blocks[i] % grid.poolnproc;
            EXPECT_EQ(start[i], i == split ? grid.poolrank * base + std::min(grid.poolrank, rem) : 0);
            EXPECT_EQ(count[i], i == split ? base + (grid.poolrank < rem) : blocks[i]);
        }

        std::vector<double> scratch;
        const double *local = layout.from_fft(potential.data(), scratch);
        EXPECT_EQ(layout.needs_transfer(), exchange);
        EXPECT_EQ(local == potential.data(), !exchange);
        EXPECT_EQ(scratch.empty(), !exchange);

        std::vector<double> contribution;
        for (int x = 0; x < count.x * grid.bx; ++x)
        {
            for (int y = 0; y < count.y * grid.by; ++y)
            {
                for (int z = 0; z < count.z * grid.bz; ++z)
                {
                    const int point = id(start.x * grid.bx + x, start.y * grid.by + y, start.z * grid.bz + z);
                    EXPECT_EQ(local[contribution.size()], point);
                    contribution.push_back(2.0 * point + 1);
                }
            }
        }
        // GPU outputs overwrite the destination. Wavefunction cubes also
        // overwrite it, with independent real and imaginary parts at each point.
        std::vector<double> overwritten(grid.nrxx, -1);
        layout.copy_to_fft(contribution.data(), overwritten.data());
        std::vector<std::complex<double>> wave(contribution.size());
        for (size_t i = 0; i < wave.size(); ++i)
        {
            wave[i] = {contribution[i], -local[i]};
        }
        std::vector<std::complex<double>> fft_wave(grid.nrxx);
        layout.copy_to_fft(wave.data(), fft_wave.data());
        std::vector<double> rho(grid.nrxx, 0.25);
        layout.add_to_fft(contribution.data(), rho.data());
        for (int i = 0; i < grid.nrxx; ++i)
        {
            EXPECT_EQ(overwritten[i], 2.0 * potential[i] + 1);
            EXPECT_EQ(fft_wave[i], std::complex<double>(overwritten[i], -potential[i]));
            EXPECT_EQ(rho[i], 2.0 * potential[i] + 1.25);
        }
    }
}

TEST(DivideInfo, MapsEveryPointInBothDirections)
{
    ModulePW::PW_Basis_Big grid("cpu", "double");
#ifdef __MPI
    MPICommGroup world(MPI_COMM_WORLD);
    grid.initmpi(world.gsize, world.grank, world.parent_comm);
#endif
    check_fields(grid, 0);
}

// A sheet spans two lattice directions. Give different pools different sheets
// so that using the wrong communicator for the axis decision is observable.
void check_selection(ModulePW::PW_Basis_Big &grid, int pool_id)
{
    const int planes = 8 * grid.poolnproc;
    grid.nbx = grid.nx = planes;
    grid.nby = grid.ny = planes;
    grid.nbz = grid.nz = planes;
    grid.nxy = planes * planes;
    grid.distribute_r();

    Atom atom;
    UnitCell cell;
    cell.atoms = &atom;
    cell.ntype = 1;
    cell.lat0 = 1.0;
    cell.omega = 1000.0;
    cell.GGT = ModuleBase::Matrix3(0.01, 0, 0, 0, 0.01, 0, 0, 0, 0.01);
    atom.Rcut = 1.5;
    atom.nw = 10;
    for (int x = 0; x < planes; ++x)
    {
        // Positions outside the unit interval represent the same periodic sheet.
        for (int y = 0; y < planes; ++y)
        {
            atom.taud.emplace_back(1.0 + (x + 0.5) / planes, -1.0 + (y + 0.5) / planes, 0.125);
        }
    }
    atom.na = atom.taud.size();
    cell.nat = atom.na;
    if (pool_id != 0)
    {
        for (auto &pos : atom.taud)
        {
            std::swap(pos.x, pos.z);
        }
    }
    SlabAxis expected = pool_id == 0 && grid.poolnproc > 1 ? SlabAxis::a : SlabAxis::c;
    EXPECT_EQ(DivideInfo::select_axis(grid, cell), expected);

    // Only b has distributed work when both a and c are localized.
    for (auto &pos : atom.taud)
    {
        pos.x = pos.z = 0.125;
    }
    expected = grid.poolnproc == 1 ? SlabAxis::c : SlabAxis::b;
    EXPECT_EQ(DivideInfo::select_axis(grid, cell), expected);

    // Periodic images of a wide cutoff cover every slab.
    atom.Rcut = 30.0;
    EXPECT_EQ(DivideInfo::select_axis(grid, cell), SlabAxis::c);

    // A very dilute orbital distribution cannot pay for moving the FFT field.
    atom.Rcut = 0.05;
    EXPECT_EQ(DivideInfo::select_axis(grid, cell), SlabAxis::c);
}

TEST(DivideInfo, SelectsAxisFromProjectedOrbitalWork)
{
    ModulePW::PW_Basis_Big grid("cpu", "double");
#ifdef __MPI
    MPICommGroup world(MPI_COMM_WORLD);
    grid.initmpi(world.gsize, world.grank, world.parent_comm);
#endif
    check_selection(grid, 0);
}

#ifdef __MPI
TEST(DivideInfo, RejectsSlabsLargerThanTheMpiCountLimit)
{
    MPICommGroup world(MPI_COMM_WORLD);
    if (world.gsize == 1)
    {
        GTEST_SKIP() << "Requires multiple ranks";
    }
    ModulePW::PW_Basis_Big grid("cpu", "double");
    grid.initmpi(world.gsize, world.grank, world.parent_comm);
    grid.nbx = grid.nby = grid.poolnproc + 1;
    const int block_width = 4096 / grid.nbx;
    if (block_width == 0)
    {
        GTEST_SKIP() << "Fixture uses at most 4096 grid points per direction";
    }
    grid.setbxyz(block_width, block_width, 1);
    grid.nx = grid.nbx * grid.bx;
    grid.ny = grid.nby * grid.by;
    grid.nxy = grid.nx * grid.ny;
    // The native slab fits, but rank 0 would own two full planes along a or b.
    grid.nbz = grid.nz = grid.poolnproc * (std::numeric_limits<int>::max() / grid.nxy);
    grid.distribute_r();
    UnitCell cell;
    EXPECT_FALSE(DivideInfo(grid, SlabAxis::a).needs_transfer());
    EXPECT_FALSE(DivideInfo(grid, SlabAxis::b).needs_transfer());
    EXPECT_EQ(DivideInfo::select_axis(grid, cell), SlabAxis::c);
}

TEST(DivideInfo, KeepsNativeLayoutWhenThereAreTooFewPlanes)
{
    MPICommGroup world(MPI_COMM_WORLD);
    if (world.gsize == 1)
    {
        GTEST_SKIP() << "Requires multiple ranks";
    }
    ModulePW::PW_Basis_Big grid("cpu", "double");
    grid.initmpi(world.gsize, world.grank, world.parent_comm);
    grid.nbx = 1;
    grid.nby = 3;
    grid.nbz = world.gsize + 1;
    grid.nxy = grid.nbx * grid.nby;
    grid.distribute_r();
    const DivideInfo layout(grid, SlabAxis::a);
    EXPECT_FALSE(layout.needs_transfer());
    EXPECT_EQ(layout.block_start().z, grid.nbzp_start);
    EXPECT_EQ(layout.block_count().z, grid.nbzp);
}

TEST(DivideInfo, UsesTheGridCommunicator)
{
    MPICommGroup pool(MPI_COMM_WORLD);
    if (pool.gsize < 4)
    {
        GTEST_SKIP() << "Requires two pools with at least two ranks each";
    }
    pool.divide_group_comm(2, false);
    ModulePW::PW_Basis_Big grid("cpu", "double");
    grid.initmpi(pool.nprocs_in_group, pool.rank_in_group, pool.group_comm);
    check_fields(grid, pool.my_group);
    check_selection(grid, pool.my_group);
}
#endif
} // namespace
