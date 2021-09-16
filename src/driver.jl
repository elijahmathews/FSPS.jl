const libfp = joinpath(@__DIR__, "../deps/fsps.so")

"""
    _is_setup()::Bool

Check the status of the `is_setup` variable.
"""
function _is_setup()
    unsafe_load(cglobal((:__driver_MOD_is_setup, libfp), Bool))
end


"""
    _setup(; compute_vega_mags = false, vactoair = false)

Setup FSPS.
"""
function _setup(; compute_vega_mags = false, vactoair = false)
    temp_compute_vega_mags = Cint[compute_vega_mags]
    temp_vactoair = Cint[vactoair]
    ccall(
        (:__driver_MOD_setup, libfp),
        Cvoid,
        (Ref{Cint}, Ref{Cint}),
        temp_compute_vega_mags,
        temp_vactoair,
    )
end


"""
    _set_ssp_params(...)
"""
function _set_ssp_params(; imf_type::Integer=2, imf_upper_limit::Real=120.0, imf_lower_limit::Real=0.08,
                         imf1::Real=1.3, imf2::Real=2.3, imf3::Real=2.3, vdmc::Real=0.08, mdave::Real=0.5,
                         dell::Real=0.0, delt::Real=0.0, sbss::Real=0.0, fbhb::Real=0.0, pagb::Real=1.0,
                         add_stellar_remnants::Integer=true, tpagb_norm_type::Integer=2,
                         add_agb_dust_model::Integer=true, agb_dust::Real=1.0, redgb::Real=1.0, agb::Real=1.0,
                         masscut::Real=150.0, fcstar::Real=1.0, evtype::Real=-1, use_wr_spectra::Integer=true,
                         logt_wmb_hot::Real=0.0, smooth_lsf::Integer=false)
    temp_imf_type = Cint[imf_type]
    temp_imf_upper_limit = Cdouble[imf_upper_limit]
    temp_imf_lower_limit = Cdouble[imf_lower_limit]
    temp_imf1 = Cdouble[imf1]
    temp_imf2 = Cdouble[imf2]
    temp_imf3 = Cdouble[imf3]
    temp_vdmc = Cdouble[vdmc]
    temp_mdave = Cdouble[mdave]
    temp_dell = Cdouble[dell]
    temp_delt = Cdouble[delt]
    temp_sbss = Cdouble[sbss]
    temp_fbhb = Cdouble[fbhb]
    temp_pagb = Cdouble[pagb]
    temp_add_stellar_remnants = Cint[add_stellar_remnants]
    temp_tpagb_norm_type = Cint[tpagb_norm_type]
    temp_add_agb_dust_model = Cint[add_agb_dust_model]
    temp_agb_dust = Cdouble[agb_dust]
    temp_redgb = Cdouble[redgb]
    temp_agb = Cdouble[agb]
    temp_masscut = Cdouble[masscut]
    temp_fcstar = Cdouble[fcstar]
    temp_evtype = Cdouble[evtype]
    temp_use_wr_spectra = Cint[use_wr_spectra]
    temp_logt_wmb_hot = Cdouble[logt_wmb_hot]
    temp_smooth_lsf = Cint[smooth_lsf]
    ccall(
        (:__driver_MOD_set_ssp_params, libfp),
        Cvoid,
        (
            Ref{Cint},    # imf_type
            Ref{Cdouble}, # imf_upper_limit
            Ref{Cdouble}, # imf_lower_limit
            Ref{Cdouble}, # imf1
            Ref{Cdouble}, # imf2
            Ref{Cdouble}, # imf3
            Ref{Cdouble}, # vdmc
            Ref{Cdouble}, # mdave
            Ref{Cdouble}, # dell
            Ref{Cdouble}, # delt
            Ref{Cdouble}, # sbss
            Ref{Cdouble}, # fbhb
            Ref{Cdouble}, # pagb
            Ref{Cint},    # add_stellar_remnants
            Ref{Cint},    # tpagb_norm_type
            Ref{Cint},    # add_agb_dust_model
            Ref{Cdouble}, # agb_dust
            Ref{Cdouble}, # redgb
            Ref{Cdouble}, # agb
            Ref{Cdouble}, # masscut
            Ref{Cdouble}, # fcstar
            Ref{Cdouble}, # evtype
            Ref{Cint},    # use_wr_spectra
            Ref{Cdouble}, # logt_wmb_hot
            Ref{Cint},    # smooth_lsf
        ),
        temp_imf_type,
        temp_imf_upper_limit,
        temp_imf_lower_limit,
        temp_imf1,
        temp_imf2,
        temp_imf3,
        temp_vdmc,
        temp_mdave,
        temp_dell,
        temp_delt,
        temp_sbss,
        temp_fbhb,
        temp_pagb,
        temp_add_stellar_remnants,
        temp_tpagb_norm_type,
        temp_add_agb_dust_model,
        temp_agb_dust,
        temp_redgb,
        temp_agb,
        temp_masscut,
        temp_fcstar,
        temp_evtype,
        temp_use_wr_spectra,
        temp_logt_wmb_hot,
        temp_smooth_lsf,
    )
end


#
# SKIP set_csp_params
#


#
# SKIP ssps
#


#
# SKIP ssp
#


#
# SKIP compute_zdep
#


#
# SKIP get_spec
#


#
# SKIP get_mags
#


#
# SKIP interp_ssp
#


#
# SKIP smooth_spectrum
#


#
# SKIP stellar_spectrum
#


#
# SKIP get_ssp_spec
#


#
# SKIP set_sfh_tab
#


#
# SKIP set_ssp_lsf
#


"""
    _get_setup_vars()
"""
function _get_setup_vars()
    temp_compute_vega_mags = zeros(Cint, 1) 
    temp_vactoair = zeros(Cint, 1)
    ccall(
        (:__driver_MOD_get_setup_vars, libfp),
        Cvoid,
        (Ref{Cint}, Ref{Cint}),
        temp_compute_vega_mags,
        temp_vactoair,
    )
    compute_vega_mags = Bool(temp_compute_vega_mags[1])
    vactoair = Bool(temp_vactoair[1])
    Dict(
        :compute_vega_mags => compute_vega_mags,
        :vactoair => vactoair,
    )
end


"""
    _get_nz()
"""
function _get_nz()
    temp_nz = zeros(Cint, 1)
    ccall(
        (:__driver_MOD_get_nz, libfp),
        Cvoid,
        (Ref{Cint},),
        temp_nz,
    )
    nz = temp_nz[1]
end


"""
    _get_zlegend()
"""
function _get_zlegend()
    nz = _get_nz()
    zlegend = zeros(Cdouble, nz)
    ccall(
        (:__driver_MOD_get_zlegend, libfp),
        Cvoid,
        (Ref{Cint}, Ref{Cdouble}),
        [nz],
        zlegend,
    )
    return zlegend
end


"""
    _get_timefull()
"""
function _get_timefull()
    ntfull = _get_ntfull()
    timefull = zeros(Cdouble, ntfull)
    ccall(
        (:__driver_MOD_get_timefull, libfp),
        Cvoid,
        (Ref{Cint}, Ref{Cdouble}),
        [ntfull],
        timefull,
    )
    return timefull
end


"""
    _get_ntfull()
"""
function _get_ntfull()
    temp_ntfull = zeros(Cint, 1)
    ccall(
        (:__driver_MOD_get_ntfull, libfp),
        Cvoid,
        (Ref{Cint},),
        temp_ntfull,
    )
    ntfull = temp_ntfull[1]
end


"""
    _get_nspec()
"""
function _get_nspec()
    temp_nspec = zeros(Cint, 1)
    ccall(
        (:__driver_MOD_get_nspec, libfp),
        Cvoid,
        (Ref{Cint},),
        temp_nspec,
    )
    nspec = temp_nspec[1]
end


"""
    _get_nbands()
"""
function _get_nbands()
    temp_nbands = zeros(Cint, 1)
    ccall(
        (:__driver_MOD_get_nbands, libfp),
        Cvoid,
        (Ref{Cint},),
        temp_nbands,
    )
    nbands = temp_nbands[1]
end


"""
    _get_nemline()
"""
function _get_nemline()
    temp_nemline = zeros(Cint, 1)
    ccall(
        (:__driver_MOD_get_nemline, libfp),
        Cvoid,
        (Ref{Cint},),
        temp_nemline,
    )
    nemline = temp_nemline[1]
end


"""
    _get_emlambda()
"""
function _get_emlambda()
    nemline = _get_nemline()
    emlambda = zeros(Cdouble, nemline)
    ccall(
        (:__driver_MOD_get_emlambda, libfp),
        Cvoid,
        (Ref{Cint}, Ref{Cdouble}),
        [nemline],
        emlambda,
    )
    return emlambda
end


"""
    _get_lambda()
"""
function _get_lambda()
    nspec = _get_nspec()
    lambda = zeros(Cdouble, nspec)
    ccall(
        (:__driver_MOD_get_lambda, libfp),
        Cvoid,
        (Ref{Cint}, Ref{Cdouble}),
        [nspec],
        lambda,
    )
    return lambda
end


"""
    _get_libraries()
"""
function _get_libraries()
    isoc_type = zeros(UInt8, 4)
    spec_type = zeros(UInt8, 5)
    str_dustem = zeros(UInt8, 6)
    ccall(
        (:__driver_MOD_get_libraries, libfp),
        Cvoid,
        (Ref{UInt8}, Ref{UInt8}, Ref{UInt8}),
        isoc_type,
        spec_type,
        str_dustem,
    )
    Dict(
        :isoc_type => String(isoc_type),
        :spec_type => String(spec_type),
        :str_dustem => String(str_dustem),
    )
end


"""
    _get_isochrone_dimensions()
"""
function _get_isochrone_dimensions()
    temp_nt = zeros(Cint, 1)
    temp_n_mass = zeros(Cint, 1)
    ccall(
        (:__driver_MOD_get_isochrone_dimensions, libfp),
        Cvoid,
        (Ref{Cint}, Ref{Cint}),
        temp_nt,
        temp_n_mass,
    )
    nt = temp_nt[1]
    n_mass = temp_n_mass[1]
    return (nt, n_mass)
end


"""
    _get_nmass_isochrone(zz, tt)
"""
function _get_nmass_isochrone(zz::Integer, tt::Integer)
    temp_zz = Cint[zz]
    temp_tt = Cint[tt]
    temp_nmass = zeros(Cint, 1)
    ccall(
        (:__driver_MOD_get_nmass_isochrone, libfp),
        Cvoid,
        (Ref{Cint}, Ref{Cint}, Ref{Cint}),
        temp_zz,
        temp_tt,
        temp_nmass,
    )
    nmass = temp_nmass[1]
end


"""
    _get_stats()
"""
function _get_stats()
    n_age = _get_ntfull()
    nline = _get_nemline()
    age = zeros(Cdouble, n_age)
    mass_csp = zeros(Cdouble, n_age)
    lbol_csp = zeros(Cdouble, n_age)
    sfr = zeros(Cdouble, n_age)
    mdust = zeros(Cdouble, n_age)
    mformed = zeros(Cdouble, n_age)
    emlines = zeros(Cdouble, n_age, nline)
    ccall(
        (:__driver_MOD_get_stats, libfp),
        Cvoid,
        (
            Ref{Cint},    # n_age
            Ref{Cint},    # nline
            Ref{Cdouble}, # age
            Ref{Cdouble}, # mass_csp
            Ref{Cdouble}, # lbol_csp
            Ref{Cdouble}, # sfr
            Ref{Cdouble}, # mdust
            Ref{Cdouble}, # mformed
            Ref{Cdouble}, # emlines
        ),
        n_age,
        nline,
        age,
        mass_csp,
        lbol_csp,
        sfr,
        mdust,
        mformed,
        emlines,
    )
    Dict(
        :age => age,
        :mass_csp => mass_csp,
        :lbol_csp => lbol_csp,
        :sfr => sfr,
        :mdust => mdust,
        :mformed => mformed,
        :emlines => emlines,
    )
end


"""
    _get_filter_data()
"""
function _get_filter_data()
    nb = _get_nbands()
    wave_eff = zeros(Cdouble, nb)
    mag_vega = zeros(Cdouble, nb)
    mag_sun = zeros(Cdouble, nb)
    ccall(
        (:__driver_MOD_get_filter_data, libfp),
        Cvoid,
        (Ref{Cint}, Ref{Cdouble}, Ref{Cdouble}, Ref{Cdouble}),
        nb,
        wave_eff,
        mag_vega,
        mag_sun,
    )
    Dict(
        :wave_eff => wave_eff,
        :mag_vega => mag_vega,
        :mag_sun => mag_sun,
    )
end


"""
    _write_isoc(outfile)
"""
function _write_isoc(outfile::Vector{UInt8})
    ccall(
        (:__driver_MOD_write_isoc, libfp),
        Cvoid,
        (Ref{UInt8},),
        outfile,
    )
end
