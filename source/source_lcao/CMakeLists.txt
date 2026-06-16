add_subdirectory(module_rt)
add_subdirectory(module_deepks)
add_subdirectory(module_dftu)
add_subdirectory(module_deltaspin)

if(ENABLE_LCAO)
    add_subdirectory(module_operator_lcao)
    add_subdirectory(module_extrap)
    set(hamilt_lcao_sources
        hamilt_lcao.cpp
        module_extrap/wf_history_lcao.cpp
        module_extrap/wf_orthonormalize_lcao.cpp
        setup_dftu_lcao.cpp
        pulay_fs_center2.cpp
        force_stress_lcao.cpp
        force_lcao_gamma.cpp
        force_lcao_k.cpp
        stress_tools.cpp
        edm.cpp
        spar_dh.cpp
        spar_exx.cpp
        spar_hsr.cpp
        spar_st.cpp
        spar_u.cpp
<<<<<<< HEAD
        lcao_set.cpp
		lcao_set_fs.cpp
		lcao_set_st.cpp
		lcao_nl_mu.cpp	
        lcao_set_zero.cpp
        lcao_allocate.cpp
        lcao_set_mat2d.cpp
        lcao_init_basis.cpp
=======
        LCAO_set.cpp
        LCAO_set_fs.cpp
        LCAO_set_st.cpp
        LCAO_nl_mu.cpp
        LCAO_set_zero.cpp
        LCAO_allocate.cpp
        LCAO_set_mat2d.cpp
        LCAO_init_basis.cpp
>>>>>>> 2a092828f (CMake cleanup of source_lcao)
        setup_nonlocal.cpp
        setup_exx.cpp
        setup_deepks.cpp
        setup_dm.cpp
        rho_tau_lcao.cpp
<<<<<<< HEAD
		record_adj.cpp
        center2orb.cpp
        center2orb_orb11.cpp
        center2orb_orb21.cpp
        center2orb_orb22.cpp
=======
        record_adj.cpp
        center2_orb.cpp
        center2_orb-orb11.cpp
        center2_orb-orb21.cpp
        center2_orb-orb22.cpp
>>>>>>> 2a092828f (CMake cleanup of source_lcao)
        wavefunc_in_pw.cpp
    )

    add_library(
        hamilt_lcao
        OBJECT
        ${hamilt_lcao_sources}
    )

    if(ENABLE_COVERAGE)
    add_coverage(hamilt_lcao)
    endif()

    IF (BUILD_TESTING)
        add_subdirectory(test)
    endif()

endif()


