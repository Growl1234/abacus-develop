#include "gint_env_k.h"

#include "gint_common.h"
#include "phi_operator.h"

#include <algorithm>
#include <cassert>

namespace ModuleGint
{

Gint_env_k::Gint_env_k(
    const std::complex<double> *psid,
    const Parallel_Orbitals *pv,
    const std::vector<Vec3d> &kvec_d,
    const int nbands,
    const int nlocal,
    const int ik,
    const int npol,
    std::complex<double> *wfc)
    : kvec_d_(kvec_d), ik_(ik), npol_(npol), wfc_(wfc)
{
    assert(npol_ == 1 || npol_ == 2);
    assert(wfc_ != nullptr);
    wfc_gint_.resize(nbands * gint_info_->get_lgd());
    wfc_2d_to_gint(psid, nbands, nlocal, *pv, wfc_gint_.data(), *gint_info_);
}

void Gint_env_k::cal_env_band(const int iband)
{
    ModuleBase::TITLE("Gint", "cal_gint_env");
    ModuleBase::timer::start("Gint", "cal_gint_env");
    const int local_mgrid_num = gint_info_->get_local_mgrid_num();
    const auto &layout = gint_info_->layout();
    const bool redistribute = layout.needs_transfer();
    std::vector<std::complex<double>> local(redistribute ? npol_ * local_mgrid_num : 0);
    auto *values = redistribute ? local.data() : wfc_;
    std::fill_n(values, npol_ * local_mgrid_num, std::complex<double>{});
    const std::complex<double> *wfc_gint_band =
        wfc_gint_.empty() ? nullptr : wfc_gint_.data() + iband * gint_info_->get_lgd();
#pragma omp parallel
    {
        PhiOperator phi_op;
        std::vector<double> phi;
#pragma omp for schedule(dynamic)
        for (int i = 0; i < gint_info_->get_bgrids_num(); i++)
        {
            const auto &biggrid = gint_info_->get_biggrids()[i];
            if (biggrid->get_atoms().size() == 0)
            {
                continue;
            }
            phi_op.set_bgrid(biggrid);
            const int phi_len = phi_op.get_rows() * phi_op.get_cols();
            phi.resize(phi_len);
            phi_op.set_phi(phi.data());
            phi_op.cal_env_k(
                phi.data(), wfc_gint_band, gint_info_->get_trace_lo(), ik_, npol_, kvec_d_, local_mgrid_num, values);
        }
    }
    if (redistribute)
    {
        for (int spinor = 0; spinor < npol_; ++spinor)
        {
            layout.copy_to_fft(values + spinor * local_mgrid_num, wfc_ + spinor * layout.fft_size());
        }
    }
    ModuleBase::timer::end("Gint", "cal_gint_env");
}

} // namespace ModuleGint
