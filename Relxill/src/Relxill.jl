module Relxill

using LibXSPEC_jll
using LibXSPEC_Relxill_jll

using SpectralFitting
using SpectralFitting: AbstractSpectralModel, Additive, FitParam

using XSPECModels: @xspecmodel, @wrap_xspec_model_ccall

@xspecmodel (:lmodrelline, libXSPEC_relxill) struct XS_Relline{T} <:
                                                    AbstractSpectralModel{T,Additive}
    "Normalisation."
    K::T
    "Rest frame emission energy."
    lineE::T
    "Emissivity index 1."
    index1::T
    "Emissivity index 2."
    index2::T
    "Break radius between the emissivity indices."
    r_break::T
    "Black hole spin."
    a::T
    "Observer inclination."
    θ_obs::T
    "Inner disc radius (set to -1 for ISCO)."
    inner_r::T
    "Outer disc radius."
    outer_r::T
    "Redshift factor."
    z::T
    "Limb darkening index."
    limb::T
end
function XS_Relline(;
    K = FitParam(1.0),
    lineE = FitParam(6.4, upper_limit = 10.0),
    index1 = FitParam(3.0, frozen = true),
    index2 = FitParam(3.0, frozen = true),
    r_break = FitParam(15.0, frozen = true),
    a = FitParam(0.998, upper_limit = 0.998, frozen = true),
    θ_obs = FitParam(30.0, upper_limit = 90.0),
    inner_r = FitParam(-1.0, frozen = true),
    outer_r = FitParam(400.0, frozen = true),
    z = FitParam(0.0, frozen = true),
    limb = FitParam(0.0, frozen = true),
)
    XS_Relline(K, lineE, index1, index2, r_break, a, θ_obs, inner_r, outer_r, z, limb)
end

@xspecmodel (:lmodrelxill, libXSPEC_relxill) struct XS_Relxill{T} <:
                                                    AbstractSpectralModel{T,Additive}
    "Normalisation."
    K::T
    "Emissivity index 1."
    index1::T
    "Emissivity index 2."
    index2::T
    "Break radius between the emissivity indices."
    r_break::T
    "Black hole spin."
    a::T
    "Observer inclination."
    θ_obs::T
    "Inner disc radius (set to -1 for ISCO)."
    inner_r::T
    "Outer disc radius."
    outer_r::T
    "Redshift factor."
    z::T
    "Photon index."
    Gamma::T
    "Log of ionisation parameter."
    logxi::T
    "Iron abundance."
    Afe::T
    "High-energy cut off"
    Ecut::T
    "Reflection fraction."
    refl_frac::T
end
function XS_Relxill(;
    K = FitParam(1.0),
    index1 = FitParam(3.0, frozen = true),
    index2 = FitParam(3.0, frozen = true),
    r_break = FitParam(15.0, frozen = true),
    a = FitParam(0.998, upper_limit = 0.998, frozen = true),
    # theta_obs upper limit of 87 wasn't strictly obeyed resulting in
    # *** relxill error : incl 87.001  is not in the required range between 3-87 deg
    # this can also happen for the lower limit
    # *** relxill error : incl 3.000  is not in the required range between 3-87 deg
    θ_obs = FitParam(30.0, lower_limit = 4.0, upper_limit = 86.0),
    inner_r = FitParam(-1.0, frozen = true),
    outer_r = FitParam(400.0, frozen = true),
    z = FitParam(0.0, frozen = true),
    Gamma = FitParam(2.0, lower_limit = 1.0, upper_limit = 3.4),
    logxi = FitParam(1.0, upper_limit = 4.7, frozen = true),
    Afe = FitParam(1.0, frozen = true),
    Ecut = FitParam(300.0, frozen = true),
    refl_frac = FitParam(1.0),
)
    XS_Relxill(
        K,
        index1,
        index2,
        r_break,
        a,
        θ_obs,
        inner_r,
        outer_r,
        z,
        Gamma,
        logxi,
        Afe,
        Ecut,
        refl_frac,
    )
end

SpectralFitting.register_model_data(XS_Relxill, "xillver/xillver-a-Ec5.fits.gz")

export XS_Relline, XS_Relxill

end # module Relxill
