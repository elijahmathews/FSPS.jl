const libfp = joinpath(@__DIR__, "..", "deps", "fsps")

"""
    is_setup()

Check if FSPS is has been initially set up or not. Returns `true` if FSPS
has been set up or `false` if FSPS has not been set up.
"""
function is_setup()
    unsafe_load(cglobal((:__driver_MOD_is_setup, libfp), Bool))
end


"""
    setup!(; <keyword arguments>)

Set up FSPS.

# Arguments
- `compute_vega_mags::Bool=false`: compute spectral flux density in Vega magnitudes rather than AB magnitudes.
- `vactoair::Bool=false`: compute wavelengths in air rather than in vacuum.

# Examples
```julia-repl
julia> using FSPS

julia> FSPS.get_setup_vars()
Dict{Symbol, Bool} with 2 entries:
  :compute_vega_mags => 0
  :vactoair          => 0

julia> FSPS.setup!(; compute_vega_mags=true)

julia> FSPS.get_setup_vars()
Dict{Symbol, Bool} with 2 entries:
  :compute_vega_mags => 1
  :vactoair          => 0
```
"""
function setup!(; compute_vega_mags::Bool=false, vactoair::Bool=false)
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
    set_ssp_params!(; <keyword arguments>)

Set the parameters that affect the simple stellar population computation.

# Arguments

- `imf_type::Integer=2`: Set the type of stellar [initial mass function](https://en.wikipedia.org/wiki/Initial_mass_function) used.
  - `imf_type=0`: Use the [Salpeter (1955)](https://doi.org/10.1086/145971) initial mass function.
  - `imf_type=1`: Use the [Chabrier (2003)](https://doi.org/10.1086/376392) initial mass function.
  - `imf_type=2`: Use the [Kroupa (2001)](https://doi.org/10.1046/j.1365-8711.2001.04022.x) initial mass function.
  - `imf_type=3`: Use the [van Dokkum (2008)](https://doi.org/10.1086/525014) initial mass function.
  - `imf_type=4`: Use the [Davé (2008)](https://doi.org/10.1111/j.1365-2966.2008.12866.x) initial mass function.
- `imf_upper_limit::Real=120`: The upper limit of the initial mass function in units of solar masses. Note that if this is above the maximum mass in the isochrones then those stars will not contribute to the spectrum but will affect the overall initial mass function normalization.
- `imf_lower_limit::Real=0.08`: The lower limit of the initial mass function in units of solar masses. Note that if this is below the minimum mass in the isochrones then those stars will not contribute to the spectrum but will affect the overall initial mass function normalization.
- `imf1::Real=1.3`: Logarithmic slope of the initial mass function over the range `0.08 < M/Msol < 0.5`. Only used if `imf_type=2`.
- `imf2::Real=2.3`: Logarithmic slope of the initial mass function over the range `0.5 < M/Msol < 1`. Only used if `imf_type=2`.
- `imf3::Real=2.3`: Logarithmic slope of the initial mass function over the range `1 < M/Msol < imf_upper_limit`. Only used if `imf_type=2`.
- `vdmc::Real=0.08`: Initial mass function parameter defined by [van Dokkum (2008)](https://doi.org/10.1086/525014). Only used if `imf_type=3`.
- `mdave::Real=0.5`: Initial mass function parameter defined by [Davé (2008)](https://doi.org/10.1111/j.1365-2966.2008.12866.x). Only used if `imf_type=4`.
- `dell::Real=0.0`: Shift in bolometric luminosity in the isochrones for the thermally-pulsating asymptotic giant branch stars. Note that the meaning of this parameter and `delt` reflect the calibrations presented by [Conroy et al. (2009)](https://doi.org/10.1088/0004-637X/699/1/486). Only used if `tpagb_norm_type=1`.
- `delt::Real=0.0`: Shift in effective temperature in the isochrones for the thermally-pulsating asymptotic giant branch stars. Note that the meaning of this parameter and `dell` reflect the calibrations presented by [Conroy et al. (2009)](https://doi.org/10.1088/0004-637X/699/1/486). Only used if `tpagb_norm_type=1`.
- `sbss::Real=0.0`: Specific frequency of blue straggler stars as parameterized by [Conroy et al. (2009)](https://doi.org/10.1088/0004-637X/699/1/486).
- `fbhb::Real=0.0`: Fraction of horizontal branch stars that are blue as parameterized by [Conroy et al. (2009)](https://doi.org/10.1088/0004-637X/699/1/486). The blue horizontal branch stars are uniformly spread in effective temperature to 10⁴ K.
- `pagb::Real=1.0`: Weight given to the post-asymptotic giant branch phase. A value of `pagb=0.0` turns off post-asymptotic giant branch stars; a value of `pagb=1.0` implies that the [Vassiliadis & Wood (1994)](https://doi.org/10.1086/191962) tracks are implemented as-is.
- `add_stellar_remnants::Integer=true`: Add stellar remnants in the stellar mass computation.
- `tpagb_norm_type::Integer=2`: Specify the thermally-pulsating asymptotic giant branch normalization scheme.
  - `tpagb_norm_type=0`: Use the Padova ([Marigo & Girardi 2007](https://doi.org/10.1051/0004-6361:20066772)) isochrones.
  - `tpagb_norm_type=1`: Use the [Conroy & Gunn (2010)](https://doi.org/10.1088/0004-637X/712/2/833) normalization.
  - `tpagb_norm_type=2`: Use the [Villaume et al. (2015)](https://doi.org/10.1088/0004-637X/806/1/82) normalization.
- `add_agb_dust_model::Bool=true`: Turn on the asymptotic giant branch circumstellar dust model presented by [Villaume (2014)](https://doi.org/2014AAS...22315115V). The asymptotic giant branch emission is scaled by the parameter `agb_dust`.
- `agb_dust::Real=1.0`: Scale the circumstellar asymptotic giant branch dust emission.
- `redgb::Real=1.0`: Modify weight given to the red giant branch. Only available with the BaSTI isochrone set.
- `agb::Real=1.0`: Modify weight given to the thermally-pulsating asymptotic giant branch. Only has an effect with FSPS version 3.1 or later.
- `masscut::Real=150.0`: Truncate the initial mass function above this value (in solar masses).
- `fcstar::Real=1.0`: Fraction of stars that the Padova isochrones identify as carbon stars that FSPS assigns to a carbon star spectrum. Set to `0.0` if all carbon stars should be regular M-type stars.
- `evtype::Integer=-1`: Compute simple stellar populations for only the given evolutionary type. All phases used when `evtype=-1`.
- `use_wr_spectra::Bool=true`: Turn on the Wolf-Rayet spectral library. If `false`, the default library will be used instead.
- `logt_wmb_hot::Real=0.0`: Use the WMBasic hot star library ([Eldridge et al. 2017](https://doi.org/10.1017/pasa.2017.51) above this value of effective temperature or 25000 K (whichever is larger).
- `smooth_lsf::Bool=false`: Apply smoothing of the simple stellar populations by a wavelength dependent line spread function. See [`set_lsf`](@ref) for more details. Only takes effect if `smooth_velocity=true`.
"""
function set_ssp_params!(; imf_type::Integer=2, imf_upper_limit::Real=120.0, imf_lower_limit::Real=0.08,
                         imf1::Real=1.3, imf2::Real=2.3, imf3::Real=2.3, vdmc::Real=0.08, mdave::Real=0.5,
                         dell::Real=0.0, delt::Real=0.0, sbss::Real=0.0, fbhb::Real=0.0, pagb::Real=1.0,
                         add_stellar_remnants::Bool=true, tpagb_norm_type::Integer=2,
                         add_agb_dust_model::Bool=true, agb_dust::Real=1.0, redgb::Real=1.0, agb::Real=1.0,
                         masscut::Real=150.0, fcstar::Real=1.0, evtype::Integer=-1, use_wr_spectra::Bool=true,
                         logt_wmb_hot::Real=0.0, smooth_lsf::Bool=false)
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


"""
    _set_csp_params!(...)
"""
function _set_csp_params!(; smooth_velocity::Integer=true, redshift_colors::Integer=false,
                          compute_light_ages::Integer=false, nebemlineinspec::Integer=true,
                          dust_type::Integer=0, add_dust_emission::Integer=true,
                          add_neb_emission::Integer=false, add_neb_continuum::Integer=true,
                          cloudy_dust::Integer=false, add_igm_absorption::Integer=false,
                          zmet::Integer=1, sfh::Integer=0, wgp1::Integer=1, wgp2::Integer=1,
                          wgp3::Integer=1, tau::Real=1.0, constant::Real=0.0, tage::Real=0.0,
                          fburst::Real=0.0, tburst::Real=11.0, dust1::Real=0.0, dust2::Real=0.0,
                          logzsol::Real=0.0, zred::Real=0.0, pmetals::Real=2.0, dust_clumps::Real=-99.0,
                          frac_nodust::Real=0.0, dust_index::Real=-0.7, dust_tesc::Real=7.0,
                          frac_obrun::Real=0.0, uvb::Real=1.0, mwr::Real=3.1, dust1_index::Real=-1.0,
                          sf_start::Real=0.0, sf_trunc::Real=0.0, sf_slope::Real=0.0,
                          duste_gamma::Real=0.01, duste_umin::Real=1.0, duste_qpah::Real=3.5,
                          sigma_smooth::Real=0.0, min_wave_smooth::Real=1e3, max_wave_smooth::Real=1e4,
                          gas_logu::Real=-2.0, gas_logz::Real=0.0, igm_factor::Real=1.0,
                          fagn::Real=0.0, agn_tau::Real=10.0)
    temp_smooth_velocity = Cint[smooth_velocity]
    temp_redshift_colors = Cint[redshift_colors]
    temp_compute_light_ages = Cint[compute_light_ages]
    temp_nebemlineinspec = Cint[nebemlineinspec]
    temp_dust_type = Cint[dust_type]
    temp_add_dust_emission = Cint[add_dust_emission]
    temp_add_neb_emission = Cint[add_neb_emission]
    temp_add_neb_continuum = Cint[add_neb_continuum]
    temp_cloudy_dust = Cint[cloudy_dust]
    temp_add_igm_absorption = Cint[add_igm_absorption]
    temp_zmet = Cint[zmet]
    temp_sfh = Cint[sfh]
    temp_wgp1 = Cint[wgp1]
    temp_wgp2 = Cint[wgp2]
    temp_wgp3 = Cint[wgp3]
    temp_tau = Cdouble[tau]
    temp_constant = Cdouble[constant]
    temp_tage = Cdouble[tage]
    temp_fburst = Cdouble[fburst]
    temp_tburst = Cdouble[tburst]
    temp_dust1 = Cdouble[dust1]
    temp_dust2 = Cdouble[dust2]
    temp_logzsol = Cdouble[logzsol]
    temp_zred = Cdouble[zred]
    temp_pmetals = Cdouble[pmetals]
    temp_dust_clumps = Cdouble[dust_clumps]
    temp_frac_nodust = Cdouble[frac_nodust]
    temp_dust_index = Cdouble[dust_index]
    temp_dust_tesc = Cdouble[dust_tesc]
    temp_frac_obrun = Cdouble[frac_obrun]
    temp_uvb = Cdouble[uvb]
    temp_mwr = Cdouble[mwr]
    temp_dust1_index = Cdouble[dust1_index]
    temp_sf_start = Cdouble[sf_start]
    temp_sf_trunc = Cdouble[sf_trunc]
    temp_sf_slope = Cdouble[sf_slope]
    temp_duste_gamma = Cdouble[duste_gamma]
    temp_duste_umin = Cdouble[duste_umin]
    temp_duste_qpah = Cdouble[duste_qpah]
    temp_sigma_smooth = Cdouble[sigma_smooth]
    temp_min_wave_smooth = Cdouble[min_wave_smooth]
    temp_max_wave_smooth = Cdouble[max_wave_smooth]
    temp_gas_logu = Cdouble[gas_logu]
    temp_gas_logz = Cdouble[gas_logz]
    temp_igm_factor = Cdouble[igm_factor]
    temp_fagn = Cdouble[fagn]
    temp_agn_tau = Cdouble[agn_tau]
    ccall(
        (:__driver_MOD_set_csp_params, libfp),
        Cvoid,
        (
            Ref{Cint},    # smooth_velocity
            Ref{Cint},    # redshift_colors
            Ref{Cint},    # compute_light_ages
            Ref{Cint},    # nebemlineinspec
            Ref{Cint},    # dust_type
            Ref{Cint},    # add_dust_emission
            Ref{Cint},    # add_neb_emission
            Ref{Cint},    # add_neb_continuum
            Ref{Cint},    # cloudy_dust
            Ref{Cint},    # add_igm_absorption
            Ref{Cint},    # zmet
            Ref{Cint},    # sfh
            Ref{Cint},    # wgp1
            Ref{Cint},    # wgp2
            Ref{Cint},    # wgp3
            Ref{Cdouble}, # tau
            Ref{Cdouble}, # constant
            Ref{Cdouble}, # tage
            Ref{Cdouble}, # fburst
            Ref{Cdouble}, # tburst
            Ref{Cdouble}, # dust1
            Ref{Cdouble}, # dust2
            Ref{Cdouble}, # logzsol
            Ref{Cdouble}, # zred
            Ref{Cdouble}, # pmetals
            Ref{Cdouble}, # dust_clumps
            Ref{Cdouble}, # frac_nodust
            Ref{Cdouble}, # dust_index
            Ref{Cdouble}, # dust_tesc
            Ref{Cdouble}, # frac_obrun
            Ref{Cdouble}, # uvb
            Ref{Cdouble}, # mwr
            Ref{Cdouble}, # dust1_index
            Ref{Cdouble}, # sf_start
            Ref{Cdouble}, # sf_trunc
            Ref{Cdouble}, # sf_slope
            Ref{Cdouble}, # duste_gamma
            Ref{Cdouble}, # duste_umin
            Ref{Cdouble}, # duste_qpah
            Ref{Cdouble}, # sigma_smooth
            Ref{Cdouble}, # min_wave_smooth
            Ref{Cdouble}, # max_wave_smooth
            Ref{Cdouble}, # gas_logu
            Ref{Cdouble}, # gas_logz
            Ref{Cdouble}, # igm_factor
            Ref{Cdouble}, # fagn
            Ref{Cdouble}, # agn_tau
        ),
        temp_smooth_velocity,
        temp_redshift_colors,
        temp_compute_light_ages,
        temp_nebemlineinspec,
        temp_dust_type,
        temp_add_dust_emission,
        temp_add_neb_emission,
        temp_add_neb_continuum,
        temp_cloudy_dust,
        temp_add_igm_absorption,
        temp_zmet,
        temp_sfh,
        temp_wgp1,
        temp_wgp2,
        temp_wgp3,
        temp_tau,
        temp_constant,
        temp_tage,
        temp_fburst,
        temp_tburst,
        temp_dust1,
        temp_dust2,
        temp_logzsol,
        temp_zred,
        temp_pmetals,
        temp_dust_clumps,
        temp_frac_nodust,
        temp_dust_index,
        temp_dust_tesc,
        temp_frac_obrun,
        temp_uvb,
        temp_mwr,
        temp_dust1_index,
        temp_sf_start,
        temp_sf_trunc,
        temp_sf_slope,
        temp_duste_gamma,
        temp_duste_umin,
        temp_duste_qpah,
        temp_sigma_smooth,
        temp_min_wave_smooth,
        temp_max_wave_smooth,
        temp_gas_logu,
        temp_gas_logz,
        temp_igm_factor,
        temp_fagn,
        temp_agn_tau,
    )
end


"""
    _ssps!()
"""
function _ssps!()
    ccall(
        (:__driver_MOD_ssps, libfp),
        Cvoid,
        (Ref{Cvoid},),
        nothing::Cvoid,
    )
end


"""
    _ssp!(zi)
"""
function _ssp!(zi::Integer)
    temp_zi = Cint[zi]
    ccall(
        (:__driver_MOD_ssp, libfp),
        Cvoid,
        (Ref{Cint},),
        temp_zi,
    )
end


"""
    compute_zdep!(zcontinuous::Integer)

Compute the full composite stellar population (and simple stellar populations if not
already cached) after interpolation in metallicity specified by `zcontinuous`.

# Arguments
- `zcontinuous::Integer`: Flag specifying how interpolation in metallicity of the simple stellar populations is performed before computing composite stellar population models.
  - `zcontinuous=0`: No interpolation, use only metallicity index specified by `zmet`.
  - `zcontinuous=1`: The simple stellar populations are interpolated to the value of `logzsol` before the spectra and magnitudes are computed, and the value of `zmet` is ignored.
  - `zcontinuous=2`: The simple stellar populations are convolved with a metallicity distribution function specified by the `logzsol` and `pmetals` parameters. The value of `zmet` is ignored.
  - `zcontinuous=3`: Use all available simple stellar population metallicities when computing the composite model, for use exclusively with tabular star formation histories where the metallicity evolution as a function of age is given (see [`set_tabular_sfh`](@ref)). The values of `zmet` and `logzsol` are ignored. Furthermore, `add_neb_emission` must be set to `false`.
"""
function compute_zdep!(zcontinuous::Integer)
    ns = _get_nspec()
    n_age = _get_ntfull()
    temp_zcontinuous = Cint[zcontinuous]
    ccall(
        (:__driver_MOD_compute_zdep, libfp),
        Cvoid,
        (Ref{Cint}, Ref{Cint}, Ref{Cint},),
        ns,
        n_age,
        temp_zcontinuous,
    )
end


"""
    _get_spec_vector()

Retrieve the spectrum `ocompsp(1)%spec` from FSPS, returns a `Vector{Cdouble}`.
Used in the case where `tage > 0`.
"""
function _get_spec_vector()
    ns = _get_nspec()
    n_age = 1
    spec_out = zeros(Cdouble, ns)
    ccall(
        (:__driver_MOD_get_spec, libfp),
        Cvoid,
        (Ref{Cint}, Ref{Cint}, Ref{Cdouble}),
        ns,
        n_age,
        spec_out,
    )
    return spec_out
end


"""
    _get_spec_matrix()

Retrieve the spectra `ocompsp(i)%spec` from FSPS for all `i`, returns a `Matrix{Cdouble}`.
Used in the case where `tage == 0`.
"""
function _get_spec_matrix()
    ns = _get_nspec()
    n_age = _get_ntfull()
    spec_out = zeros(Cdouble, ns, n_age)
    ccall(
        (:__driver_MOD_get_spec, libfp),
        Cvoid,
        (Ref{Cint}, Ref{Cint}, Ref{Cdouble}),
        ns,
        n_age,
        spec_out,
    )
    return spec_out
end


"""
    _get_mags(zred, mag_compute)
"""
function _get_mags(zred::Real, mag_compute::Vector{I}) where {I <: Integer}
    ns = _get_nspec()
    n_age = _get_ntfull()
    n_bands = _get_nbands()
    temp_zred = Cdouble[zred]
    temp_mag_compute = Vector{Cint}(mag_compute)
    mags = zeros(Cdouble, n_age, n_bands)
    ccall(
        (:__driver_MOD_get_mags, libfp),
        Cvoid,
        (Ref{Cint}, Ref{Cint}, Ref{Cint}, Ref{Cdouble}, Ref{Cint}, Ref{Cdouble}),
        ns,
        n_age,
        n_bands,
        temp_zred,
        temp_mag_compute,
        mags,
    )
    return mags
end


"""
    _interp_ssp(zpos, tpos)
"""
function _interp_ssp(zpos::Real, tpos::Real)
    ns = _get_nspec()
    temp_zpos = Cdouble[zpos]
    temp_tpos = Cdouble[tpos]
    spec = zeros(Cdouble, ns)
    mass = Cdouble[0.0]
    lbol = Cdouble[0.0]
    ccall(
        (:__driver_MOD_interp_ssp, libfp),
        Cvoid,
        (Ref{Cint}, Ref{Cdouble}, Ref{Cdouble}, Ref{Cdouble}, Ref{Cdouble}, Ref{Cdouble}),
        ns,
        temp_zpos,
        temp_tpos,
        spec,
        mass,
        lbol,
    )
    Dict(
        :spec => spec,
        :mass => mass,
        :lbol => lbol,
    )
end


"""
    _smooth_spectrum(wave, spec, sigma_broad, minw, maxw)
"""
function _smooth_spectrum!(ns::Integer, wave::Vector{R}, spec::Vector{R}, sigma_broad::Real,
                           minw::Real, maxw::Real) where {R <: Real}
    temp_ns = Cint[ns]
    temp_wave = Vector{Cdouble}(wave)
    temp_spec = Vector{Cdouble}(spec)
    temp_sigma_broad = Cdouble[sigma_broad]
    temp_minw = Cdouble[minw]
    temp_maxw = Cdouble[maxw]
    ccall(
        (:__driver_MOD_smooth_spectrum, libfp),
        Cvoid,
        (Ref{Cint}, Ref{Cdouble}, Ref{Cdouble}, Ref{Cdouble}, Ref{Cdouble}, Ref{Cdouble}),
        temp_ns,
        temp_wave,
        temp_spec,
        temp_sigma_broad,
        temp_minw,
        temp_maxw,
    )
    spec = temp_spec
end


"""
    _stellar_spectrum!(...)
"""
function _stellar_spectrum!(mact::Real, logt::Real, lbol::Real, logg::Real, phase::Real,
                            ffco::Real, lmdot::Real, wght::Real, spec_out::Vector{R}) where {R <: Real}
    ns = _get_nspec()
    temp_mact = Cdouble[mact]
    temp_logt = Cdouble[logt]
    temp_lbol = Cdouble[lbol]
    temp_logg = Cdouble[logg]
    temp_phase = Cdouble[phase]
    temp_ffco = Cdouble[ffco]
    temp_lmdot = Cdouble[lmdot]
    temp_wght = Cdouble[wght]
    temp_spec_out = Vector{Cdouble}(spec_out)
    ccall(
        (:__driver_MOD_stellar_spectrum, libfp),
        Cvoid,
        (
         Ref{Cint},    # ns
         Ref{Cdouble}, # mact
         Ref{Cdouble}, # logt
         Ref{Cdouble}, # lbol
         Ref{Cdouble}, # logg
         Ref{Cdouble}, # phase
         Ref{Cdouble}, # ffco
         Ref{Cdouble}, # lmdot
         Ref{Cdouble}, # wght
         Ref{Cdouble}, # spec_out
        ),
        ns,
        temp_mact,
        temp_logt,
        temp_lbol,
        temp_logg,
        temp_phase,
        temp_ffco,
        temp_lmdot,
        temp_wght,
        temp_spec_out,
    )
    spec_out = temp_spec_out
end


"""
    _get_ssp_spec!(ssp_spec_out, ssp_mass_out, ssp_lbol_out)
"""
function _get_ssp_spec!(ssp_spec_out::Matrix{R}, ssp_mass_out::Vector{R}, ssp_lbol_out::Vector{R}) where {R <: Real}
    ns = _get_nspec()
    n_age = _get_ntfull()
    n_z = _get_nz()
    temp_ssp_spec_out = Matrix{Cdouble}(ssp_spec_out)
    temp_ssp_mass_out = Vector{Cdouble}(ssp_mass_out)
    temp_ssp_lbol_out = Vector{Cdouble}(ssp_lbol_out)
    ccall(
        (:__driver_MOD_get_ssp_spec, libfp),
        Cvoid,
        (Ref{Cint}, Ref{Cint}, Ref{Cint}, Ref{Cdouble}, Ref{Cdouble}, Ref{Cdouble}),
        ns,
        n_age,
        n_z,
        temp_ssp_spec_out,
        temp_ssp_mass_out,
        temp_ssp_lbol_out,
    )
    ssp_spec_out = temp_ssp_spec_out
    ssp_mass_out = temp_ssp_mass_out
    ssp_lbol_out = temp_ssp_lbol_out
end


"""
    _set_sfh_tab!(ntab, age, sfr, met)
"""
function _get_sfh_tab!(ntab::Integer, age::Vector{R}, sfr::Vector{R}, met::Vector{R}) where {R <: Real}
    temp_ntab = Cint[ntab]
    temp_age = Vector{Cdouble}(age)
    temp_sfr = Vector{Cdouble}(sfr)
    temp_met = Vector{Cdouble}(met)
    ccall(
        (:__driver_MOD_set_sfh_tab, libfp),
        Cvoid,
        (Ref{Cint}, Ref{Cdouble}, Ref{Cdouble}, Ref{Cdouble}),
        temp_ntab,
        temp_age,
        temp_sfr,
        temp_met,
    )
end


"""
    _set_ssp_lsf!(nsv, sigma, wlo, whi)
"""
function _set_ssp_lsf!(nsv::Integer, sigma::Vector{R}, wlo::Real, whi::Real) where {R <: Real}
    temp_nsv = Cint[nsv]
    temp_sigma = Vector{Cdouble}(sigma)
    temp_wlo = Cdouble[wlo]
    temp_whi = Cdouble[whi]
    ccall(
        (:__driver_MOD_set_ssp_lsf, libfp),
        Cvoid,
        (Ref{Cint}, Ref{Cdouble}, Ref{Cdouble}, Ref{Cdouble}),
        temp_nsv,
        temp_sigma,
        temp_wlo,
        temp_whi,
    )
end


"""
    get_setup_vars()

Retrieve the `compute_vega_mags` and `vactoair` values FSPS was set up with.

# Examples
```julia-repl
julia> using FSPS

julia> FSPS.get_setup_vars()
Dict{Symbol, Bool} with 2 entries:
  :compute_vega_mags => 0
  :vactoair          => 0

julia> FSPS.setup!(; compute_vega_mags=true)

julia> FSPS.get_setup_vars()
Dict{Symbol, Bool} with 2 entries:
  :compute_vega_mags => 1
  :vactoair          => 0
```
"""
function get_setup_vars()
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
