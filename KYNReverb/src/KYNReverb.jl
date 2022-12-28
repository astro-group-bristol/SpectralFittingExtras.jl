module KYNReverb

using LibXSPEC_KYNReverb_jll

using SpectralFitting
using SpectralFitting:
    @xspecmodel,
    @wrap_xspec_model_ccall,
    AbstractSpectralModel,
    Additive,
    FitParam,
    parameter_type,
    FreeParameters

@xspecmodel (:KYNrefrev, libXSPEC_kynrefrev) struct XS_KYNRefrev{T,F} <:
                                                    AbstractSpectralModel{Additive}
    "Normalisation."
    K::FitParam{T}
    "Black hole spin."
    a::FitParam{T}
    "Observer inclination angle (degrees)."
    θₒ::FitParam{T}
    "Inner disc radius."
    r_inner::FitParam{T}
    "Switch for inner edge: 0: integrate from inner edge. 1: If inner edge below marginally stable, integrate above marginally stable only. 2: Integrate from inner edge in units of marginally stable."
    ms::FitParam{T}
    "Outer disc radius."
    r_outer::FitParam{T}
    "Lower azimuth of non-zero disc emissivity (degrees)."
    ϕ::FitParam{T}
    dϕ::FitParam{T}
    "Black hole mass in units of 10⁸ solar masses."
    M::FitParam{T}
    "Coronal source height."
    h::FitParam{T}
    "Power law index of the primary flux."
    Γ::FitParam{T}
    "Intrinsic (if negative) or observed (positive) primary isotropic flux in the X-ray energy range 2-10 keV in units of L/Ledd."
    L::FitParam{T}
    "Ratio of primary to reflected flux. 0: Primary source hidden (only reflection). 1: Self-consistent model, with positive towards observer, negative towards the disc."
    NpNr::FitParam{T}
    "Density profile normalisation in 10^15 cm⁻³ if positive. Ionisation profile if negative."
    density::FitParam{T}
    "Density profile in `abs(density) × r^density_profile`."
    density_profile::FitParam{T}
    "Fe abundance in solar units."
    Fe_abundance::FitParam{T}
    "Fraction of thermalised flux from overal incident flux. 0: Only reverberation of reflection radiation. Negative: only thermal radiation. Positive: Both."
    thermalised_flux::FitParam{T}
    "Accretion rate in units of Ledd (positive) or solar mass per Julian year (negative)."
    accretion_rate::FitParam{T}
    "Spectral hardenning factor."
    f_column::FitParam{T}
    "Position of cloud centre in α impact parameter."
    α_cloud::FitParam{T}
    "Position of cloud centre in β impact parameter."
    β_cloud::FitParam{T}
    "Radius of the obscuring cloud."
    r_cloud::FitParam{T}
    "Overall Doppler shift."
    z::FitParam{T}
    "Limb darkening. 0: isotropic emission. 1: Laor. 2: Haardt."
    limb::FitParam{T}
    "Reflection table to use. 1: Reflion. 2. Reflionx."
    table::FitParam{T}
    "Switch for how to compute reflection spectral. 1: Use computed ionisation parameter ξ, for interpolation in reflion use proper total incident intensity with shifted cut-offs. 2: Use ξ, correspondent to the computed normalisation of the incident flux (i.e. do no shift the cut-offs when computing total incident intensity)."
    switch::FitParam{T}
    "Which FITS table to use. Current ntable=80 is correct for this model."
    n_table::FitParam{T}
    "Number of grid points in radius. If negative, then dependent on the height."
    n_radius::FitParam{T}
    "Divisions of radius. 0: equidistant, 1: expontential, >1: mixed with constant log step in inner regius and constant linear step in outer region."
    r_division::FitParam{T}
    "Number of grid points to use in azimuth."
    n_ϕ::FitParam{T}
    "Length of the time bin (in gravitational units)."
    dt::FitParam{T}
    "Number of time subbins per one time bin"
    nt::FitParam{T}
    t1_f1_e1::FitParam{T}
    t2_f2_e2::FitParam{T}
    "Lower value of the energy refence band for lag or amplitude energy dependence."
    Eref_lower::FitParam{T}
    "Upper value of the energy refence band for lag or amplitude energy dependence."
    Eref_upper::FitParam{T}
    "Lag shift for lag-energy dependence."
    dt_Af::FitParam{T}
    "Multiplicative factor for the amplitude-energy dependence."
    k_qf::FitParam{T}
    """
    - defines output in the XSPEC (photar array)
            - 0: spectrum for time interval defined by par32 and par33
        - _the following values correspond to energy dependent Fourier transform 
        at the frequency band defined by par32 and par33:_
            - -1: real part of FT of the relative reflection
            - -2: imaginary part of FT of the relative reflection
            - -3: amplitude of FT of the relative reflection
            - -4: phase of FT of the relative reflection
            - -5: amplitude for the relative reflection divided by amplitude in the 
                reference energy band defined by par34 and par35  (integration in 
                frequencies is done in real and imaginary parts first and then 
                the amplitudes are computed)
            - -6: lag for the relative reflection with respect to reference energy 
                band defined by par34 and par35 (integration in frequencies is 
                done in real and imaginary parts first and then the lags are 
                computed with frequency at half of the wrapping frequency or 
                middle of the frequency band)
            - -7: amplitude  for the relative reflection divided by amplitude in 
                the reference energy band defined by par34 and par35 (integration 
                in frequencies here is done in amplitudes directly)
            - -8: lag for the relative reflection with respect to reference energy 
                band defined by par34 and par35 (integration in frequencies here 
                is done in lags directly)
            - 1: real part of FT including primary radiation
            - 2: imaginary part of FT including primary radiation
            - 3: amplitude of FT including primary radiation
            - 4: phase of FT including primary radiation
            - 5: amplitude including the primary radiation divided by amplitude in 
                the reference energy band defined by par34 and par35 (integration 
                in frequencies is done in real and imaginary parts first and then 
                the amplitudes are computed)
            - 6: lag diluted by primary radiation with respect to reference energy 
                band defined by par34 and par35 (integration in frequencies is 
                done in real and imaginary parts first and then the lags are 
                computed with frequency at half of the wrapping frequency or 
                middle of the frequency band)
            - 7: amplitude including the primary radiation divided by amplitude in 
                the reference energy band defined by par34 and par35 (integration 
                in frequencies here is done in amplitudes directly)
            - 8: lag diluted by primary radiation with respect to reference energy 
                band defined by par34 and par35 (integration in frequencies here 
                is done in lags directly)
        - _the following values correspond to frequency dependent Fourier 
        transform for the energy band of interest defined by par32 and par33:_
            - -11: real part of FT of the relative reflection
            - -12: imaginary part of FT of the relative reflection
            - -13: amplitude of FT of the relative reflection
            - -14: phase of FT of the relative reflection
            - -15: amplitude  for the relative reflection divided by amplitude in 
                the reference energy band defined by par34 and par35 (rebinning 
                here is done in real and imaginary parts first and then the 
                amplitudes are computed)
            - -16: lag for the relative reflection with respect to reference energy 
                band defined by par34 and par35 (rebinning here is done in real 
                and imaginary parts first and then the lags are computed)
            - -17: amplitude  for the relative reflection divided by amplitude in 
                the reference energy band defined by par34 and par35 (rebinning 
                here is done in amplitudes directly)
            - -18: lag for the relative reflection with respect to reference energy 
                band defined by par34 and par35 (rebinning here is done in lags 
                directly)
            - 11: real part of FT including primary radiation
            - 12: imaginary part of FT including primary radiation
            - 13: amplitude of FT including primary radiation
            - 14: phase of FT including primary radiation
            - 15: amplitude including the primary radiation divided by amplitude in 
                the reference energy band defined by par34 and par35 (rebinning 
                here is done in real and imaginary parts first and then the 
                amplitudes are computed)
            - 16: lag diluted by primary radiation with respect to reference energy 
                band defined by par34 and par35 (rebinning here is done in real 
                and imaginary parts first and then the lags are computed)
            - 17: amplitude including the primary radiation divided by amplitude in 
                the reference energy band defined by par34 and par35 (rebinning 
                here is done in amplitudes directly)
            - 18: lag diluted by primary radiation with respect to reference energy 
                band defined by par34 and par35 (rebinning here is done in lags 
                directly)
    """
    output_type::FitParam{T}
    "Number of threads."
    n_threads::FitParam{T}
    function XS_KYNRefrev(;
        K = FitParam(1.0),
        a = FitParam(1.0, lower_limit = 0.0, upper_limit = 1.0),
        θₒ = FitParam(30.0, lower_limit = 0.0, upper_limit = 90.0),
        r_inner = FitParam(1.0),
        ms = FitParam(1.0),
        r_outer = FitParam(1000.0),
        ϕ = FitParam(0.0),
        dϕ = FitParam(360.0),
        M = FitParam(0.1, lower_limit = 1e-8, upper_limit = 1e3),
        h = FitParam(3.0),
        Γ = FitParam(2.0),
        L = FitParam(0.001),
        NpNr = FitParam(1.0),
        density = FitParam(1.0),
        density_profile = FitParam(0.0),
        Fe_abundance = FitParam(1.0),
        thermalised_flux = FitParam(0.0),
        accretion_rate = FitParam(0.1),
        f_column = FitParam(2.4),
        α_cloud = FitParam(-6.0),
        β_cloud = FitParam(0.0),
        r_cloud = FitParam(0.0),
        z = FitParam(0.0),
        limb = FitParam(0.0),
        table = FitParam(2.0),
        switch = FitParam(1.0),
        n_table = FitParam(80.0),
        n_radius = FitParam(-4488.0),
        r_division = FitParam(-1.0),
        n_ϕ = FitParam(180.0),
        dt = FitParam(1.0),
        nt = FitParam(1.0),
        t1_f1_e1 = FitParam(0.3),
        t2_f2_e2 = FitParam(0.8),
        Eref_lower = FitParam(1.0),
        Eref_upper = FitParam(3.0),
        dt_Af = FitParam(0.0),
        k_qf = FitParam(1.0),
        output_type = FitParam(16.0),
        n_threads = FitParam(1.0),
    )
        new{parameter_type(K),FreeParameters{(:K, :a, :θₒ)}}(
            K,
            a,
            θₒ,
            r_inner,
            ms,
            r_outer,
            ϕ,
            dϕ,
            M,
            h,
            Γ,
            L,
            NpNr,
            density,
            density_profile,
            Fe_abundance,
            thermalised_flux,
            accretion_rate,
            f_column,
            α_cloud,
            β_cloud,
            r_cloud,
            z,
            limb,
            table,
            switch,
            n_table,
            n_radius,
            r_division,
            n_ϕ,
            dt,
            nt,
            t1_f1_e1,
            t2_f2_e2,
            Eref_lower,
            Eref_upper,
            dt_Af,
            k_qf,
            output_type,
            n_threads,
        )
    end
end

function __init__()
    # set the KYDIR for the model data directory
    ccall(
        (:FPMSTR, LibXSPEC_jll.libXSFunctions),
        Cint,
        (Cstring, Cstring),
        "KYDIR",
        LibXSPEC_KYNReverb_jll.artifact_dir * "/lib",
    )
end


export XS_KYNRefrev

end # module
