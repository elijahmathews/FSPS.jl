Base.@kwdef struct StellarPopulation
    add_agb_dust_model::Bool = true
    add_dust_emission::Bool = true
    add_igm_absorption::Bool = false
    add_neb_emission::Bool = false
    add_neb_continuum::Bool = true
    add_stellar_remnants::Bool = true
    redshift_colors::Bool = false
    compute_light_ages::Bool = false
    nebemlineinspec::Bool = true
    smooth_velocity::Bool = true
    smooth_lsf::Bool = false
    cloudy_dust::Bool = false
    agb_dust::Real = 1.0
    tpagb_norm_type::Integer = 2
    dell::Real = 0.0
    delt::Real = 0.0
    redgb::Real = 1.0
    agb::Real = 1.0
    fcstar::Real = 1.0
    fbhb::Real = 0.0
    sbss::Real = 0.0
    pagb::Real = 1.0
    zred::Real = 0.0
    zmet::Integer = 1
    logzsol::Real = 0.0
    pmetals::Real = 2.0
    imf_type::Integer = 2
    imf_upper_limit::Real = 120.0
    imf_lower_limit::Real = 0.08
    imf1::Real = 1.3
    imf2::Real = 2.3
    imf3::Real = 2.3
    vdmc::Real = 0.08
    mdave::Real = 0.5
    evtype::Integer = -1
    use_wr_spectra::Bool = true
    logt_wmb_hot::Real = 0.0
    masscut::Real = 150.0
    sigma_smooth::Real = 0.0
    min_wave_smooth::Real = 1.0e3
    max_wave_smooth::Real = 1.0e4
    gas_logu::Real = -2.0
    gas_logz::Real = 0.0
    igm_factor::Real = 1.0
    sfh::Integer = 0
    tau::Real = 1.0
    constant::Real = 0.0
    sf_start::Real = 0.0
    sf_trunc::Real = 0.0
    tage::Real = 0.0
    dust_tesc::Real = 7.0
    fburst::Real = 0.0
    tburst::Real = 11.0
    sf_slope::Real = 0.0
    dust_type::Integer = 0
    dust1::Real = 0.0
    dust2::Real = 0.0
    dust_clumps::Real = -99.0
    frac_nodust::Real = 0.0
    frac_obrun::Real = 0.0
    dust_index::Real = -0.7
    dust1_index::Real = -1.0
    mwr::Real = 3.1
    uvb::Real = 1.0
    wgp1::Integer = 1
    wgp2::Integer = 1
    wgp3::Integer = 1
    duste_gamma::Real = 0.01
    duste_umin::Real = 1.0
    duste_qpad::Real = 3.5
    fagn::Real = 0.0
    agn_tau::Real = 10.0
end
