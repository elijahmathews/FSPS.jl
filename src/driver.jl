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


#
# SKIP set_ssp_params
#


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
            Ref{Cint},
            Ref{Cint},
            Ref{Cdouble},
            Ref{Cdouble},
            Ref{Cdouble},
            Ref{Cdouble},
            Ref{Cdouble},
            Ref{Cdouble},
            Ref{Cdouble},
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
