#ifndef ABACUS_SOURCE_HAMILT_MODULE_GINT_DIVIDE_INFO_H
#define ABACUS_SOURCE_HAMILT_MODULE_GINT_DIVIDE_INFO_H

#include "source_base/vector3.h"
#include <vector>

class UnitCell;

namespace ModulePW
{
class PW_Basis_Big;
}

namespace ModuleGint
{

enum class SlabAxis
{
    a,
    b,
    c
};

// Ownership of whole integration blocks and transfers to/from the FFT's c slabs.
class DivideInfo
{
  public:
    DivideInfo(const ModulePW::PW_Basis_Big &grid, SlabAxis axis);
    // Collective on the grid's pool. Call before constructing a layout for grid integration.
    static SlabAxis select_axis(const ModulePW::PW_Basis_Big &grid, const UnitCell &cell);

    const ModuleBase::Vector3<int> &block_start() const { return block_start_; }
    const ModuleBase::Vector3<int> &block_count() const { return block_count_; }
    bool needs_transfer() const { return !gint_.counts.empty(); }

    // In the FFT layout return the original field, leaving scratch unused.
    const double *from_fft(const double *field, std::vector<double> &scratch) const;
    template <typename T>
    void copy_to_fft(const T *local, T *field) const;
    int fft_size() const;
    void add_to_fft(const double *contribution, double *field) const;

  private:
    // A peer owns equally spaced contiguous runs within the local field.
    struct Rows
    {
        int first;
        int number;
        int stride;
        int length;
    };
    struct Layout
    {
        void append(Rows runs);
        std::vector<Rows> rows;
        std::vector<int> counts;
        std::vector<int> offsets;
        int size = 0;
        bool contiguous = true;
    };

    template <typename T>
    void transfer(
        const T *src,
        const Layout &send,
        T *dst,
        const Layout &recv,
        bool add) const;
    const ModulePW::PW_Basis_Big &grid_;
    ModuleBase::Vector3<int> block_start_;
    ModuleBase::Vector3<int> block_count_;
    Layout fft_;
    Layout gint_;
};

} // namespace ModuleGint
#endif
