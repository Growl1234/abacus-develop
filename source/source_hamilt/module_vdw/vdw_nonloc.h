#ifndef VDW_NONLOC_H
#define VDW_NONLOC_H

#include <string>

namespace vdw
{
inline bool is_nonloc(const std::string &method)
{
    return method == "vdw_df1" || method == "vdw_df2" || method == "vdw_dfcx" || method == "rvv10";
}
} // namespace vdw

#ifdef __LIBVDWXC
#include "source_base/matrix.h"
#include "source_base/vector3.h"

#include <complex>
#include <memory>
#include <vector>

struct vdwxc_data_obj;
namespace ModulePW
{
class PW_Basis;
}

namespace vdw
{
// One fixed-cell SCF calculation. Both the FFT plan and its work arrays survive density updates.
class Nonloc
{
  public:
    Nonloc(
        const std::string &method,
        double b,
        double c,
        int nspin,
        const ModulePW::PW_Basis &basis);
    ~Nonloc();

    // Valence spin densities and unpolarized NLCC core density on the charge grid.
    // Add the cell energy, valence potential integral and local SCF potential in Ry.
    void add_xc(
        const double *const *rho,
        const double *core,
        double &etxc,
        double &vtxc,
        ModuleBase::matrix &v);

    // Request the Cartesian stress correction (Ry/Bohr^3) with the next potential update.
    void request_stress(bool needed);
    const ModuleBase::matrix &stress() const;

    // Add the spin-averaged nonlocal potential in Ry on the last charge grid.
    void add_core_potential(double *vcore) const;

  private:
    struct Slabs;
    const ModulePW::PW_Basis &basis_;
    const int channels_;
    const bool total_spin_;
    vdwxc_data_obj *handle_;
    std::unique_ptr<Slabs> slabs_;
    ModuleBase::matrix stress_;
    std::vector<double> fields_;
    std::vector<double> deriv_;
    std::vector<ModuleBase::Vector3<double>> gradient_;
    std::vector<std::complex<double>> rhog_;
};
} // namespace vdw
#endif
#endif
