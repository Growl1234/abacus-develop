#include "gint_interface.h"
#include "gint_drho.h"
#include "gint_fvl.h"
#include "gint_fvl_meta.h"
#include "gint_rho.h"
#include "gint_tau.h"
#include "gint_vl.h"
#include "gint_vl_metagga.h"
#include "gint_vl_mgga_n4.h"
#include "gint_vl_nspin4.h"

#ifdef __CUDA
#include "gint_fvl_gpu.h"
#include "gint_fvl_meta_gpu.h"
#include "gint_rho_gpu.h"
#include "gint_tau_gpu.h"
#include "gint_vl_gpu.h"
#include "gint_vl_metagga_gpu.h"
#include "gint_vl_mgga_n4_gpu.h"
#include "gint_vl_nspin4_gpu.h"
#endif

namespace ModuleGint
{
#ifndef __CUDA
// CPU-only builds name the GPU kernels in run_gint but never instantiate them.
class Gint_vl_gpu;
class Gint_vl_nspin4_gpu;
class Gint_vl_metagga_gpu;
class Gint_vl_metagga_nspin4_gpu;
class Gint_fvl_gpu;
class Gint_fvl_meta_gpu;
#endif

namespace
{
template <typename CpuKernel, typename GpuKernel, typename... Args>
void run_gint(const Args &...args)
{
#ifdef __CUDA
    if (Gint::get_gint_info().use_gpu())
    {
        GpuKernel(args...).cal_gint();
        return;
    }
#endif
    CpuKernel(args...).cal_gint();
}

// Both CPU and GPU kernels take host buffers in Gint order.
class InputFields
{
  public:
    InputFields(const std::vector<const double *> &fft, size_t count) : fft_(fft)
    {
        const auto &layout = Gint::get_gint_info().layout();
        if (layout.needs_transfer())
        {
            buffers_.resize(count);
            for (size_t i = 0; i < count; ++i)
            {
                local_.push_back(layout.from_fft(fft[i], buffers_[i]));
            }
        }
    }

    const std::vector<const double *> &data() const { return local_.empty() ? fft_ : local_; }

  private:
    const std::vector<const double *> &fft_;
    std::vector<std::vector<double>> buffers_;
    std::vector<const double *> local_;
};

class OutputFields
{
  public:
    OutputFields(double **fft, int count, const GintInfo &info) : fft_(fft), info_(info)
    {
        if (info_.layout().needs_transfer())
        {
            const int points = info_.get_local_mgrid_num();
            buffers_.resize(count);
            for (auto &buffer : buffers_)
            {
                buffer.resize(points);
                local_.push_back(buffer.data());
            }
        }
    }

    double **data() { return local_.empty() ? fft_ : local_.data(); }

    void add_to_fft() const
    {
        for (size_t i = 0; i < local_.size(); ++i)
        {
            info_.layout().add_to_fft(local_[i], fft_[i]);
        }
    }

    void copy_to_fft() const
    {
        for (size_t i = 0; i < local_.size(); ++i)
        {
            info_.layout().copy_to_fft(local_[i], fft_[i]);
        }
    }

  private:
    double **fft_;
    const GintInfo &info_;
    std::vector<std::vector<double>> buffers_;
    std::vector<double *> local_;
};
} // namespace

void cal_gint_vl(const double *vr_eff, HContainer<double> *hR)
{
    std::vector<double> vr_eff_storage;
    vr_eff = Gint::get_gint_info().layout().from_fft(vr_eff, vr_eff_storage);

    run_gint<Gint_vl, Gint_vl_gpu>(vr_eff, hR);
}

// nspin == 4 case
void cal_gint_vl(std::vector<const double *> vr_eff, HContainer<std::complex<double>> *hR)
{
    InputFields local_vr_eff(vr_eff, vr_eff.size());

    run_gint<Gint_vl_nspin4, Gint_vl_nspin4_gpu>(local_vr_eff.data(), hR);
}

void cal_gint_vl_metagga(const double *vr_eff, const double *vfork, HContainer<double> *hR)
{
    const auto &layout = Gint::get_gint_info().layout();
    std::vector<double> vr_eff_storage;
    std::vector<double> vfork_storage;
    vr_eff = layout.from_fft(vr_eff, vr_eff_storage);
    vfork = layout.from_fft(vfork, vfork_storage);

    run_gint<Gint_vl_metagga, Gint_vl_metagga_gpu>(vr_eff, vfork, hR);
}

// nspin == 4 case
void cal_gint_vl_metagga(
    std::vector<const double *> vr_eff, std::vector<const double *> vofk, HContainer<std::complex<double>> *hR)
{
    InputFields local_vr_eff(vr_eff, vr_eff.size());
    InputFields local_vofk(vofk, vofk.size());

    run_gint<Gint_vl_metagga_nspin4, Gint_vl_metagga_nspin4_gpu>(local_vr_eff.data(), local_vofk.data(), hR);
}

void cal_gint_rho(
    const std::vector<HContainer<double> *> &dm_vec,
    const int nspin,
    double **rho,
    bool is_dm_symm)
{
    const auto &info = Gint::get_gint_info();
    OutputFields local_rho(rho, nspin, info);

#ifdef __CUDA
    if (info.use_gpu())
    {
        Gint_rho_gpu gint_rho(dm_vec, nspin, local_rho.data(), is_dm_symm);
        gint_rho.cal_gint();
        local_rho.copy_to_fft(); // Device kernels overwrite the host field.
        return;
    }
#endif
    Gint_rho gint_rho(dm_vec, nspin, local_rho.data(), is_dm_symm);
    gint_rho.cal_gint();
    local_rho.add_to_fft();
}

void cal_gint_drho(
    const std::vector<HContainer<double> *> &dm_vec,
    const int nspin,
    double **drho_x,
    double **drho_y,
    double **drho_z)
{
    const auto &info = Gint::get_gint_info();
    OutputFields local_drho_x(drho_x, nspin, info);
    OutputFields local_drho_y(drho_y, nspin, info);
    OutputFields local_drho_z(drho_z, nspin, info);

    // CPU/fp64 only: set_phi_dphi provides double-precision gradients (no GPU/fp32 path).
    Gint_drho gint_drho(dm_vec, nspin, local_drho_x.data(), local_drho_y.data(), local_drho_z.data());
    gint_drho.cal_gint();
    local_drho_x.add_to_fft();
    local_drho_y.add_to_fft();
    local_drho_z.add_to_fft();
}

void cal_gint_tau(const std::vector<HContainer<double> *> &dm_vec, const int nspin, double **tau)
{
    const auto &info = Gint::get_gint_info();
    OutputFields local_tau(tau, nspin, info);

#ifdef __CUDA
    if (info.use_gpu())
    {
        Gint_tau_gpu gint_tau(dm_vec, nspin, local_tau.data());
        gint_tau.cal_gint();
        local_tau.copy_to_fft(); // Device kernels overwrite the host field.
        return;
    }
#endif
    Gint_tau gint_tau(dm_vec, nspin, local_tau.data());
    gint_tau.cal_gint();
    local_tau.add_to_fft();
}

void cal_gint_fvl(
    const int nspin,
    const std::vector<const double *> &vr_eff,
    const std::vector<HContainer<double> *> &dm_vec,
    const bool isforce,
    const bool isstress,
    ModuleBase::matrix *fvl,
    ModuleBase::matrix *svl)
{
    InputFields local_vr_eff(vr_eff, nspin);

    run_gint<Gint_fvl, Gint_fvl_gpu>(nspin, local_vr_eff.data(), dm_vec, isforce, isstress, fvl, svl);
}

void cal_gint_fvl_meta(
    const int nspin,
    const std::vector<const double *> &vr_eff,
    const std::vector<const double *> &vofk,
    const std::vector<HContainer<double> *> &dm_vec,
    const bool isforce,
    const bool isstress,
    ModuleBase::matrix *fvl,
    ModuleBase::matrix *svl)
{
    InputFields local_vr_eff(vr_eff, nspin);
    InputFields local_vofk(vofk, nspin);

    run_gint<Gint_fvl_meta, Gint_fvl_meta_gpu>(
        nspin, local_vr_eff.data(), local_vofk.data(), dm_vec, isforce, isstress, fvl, svl);
}

} // namespace ModuleGint
