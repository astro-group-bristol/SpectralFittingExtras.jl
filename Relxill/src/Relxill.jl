module Relxill

using LibXSPEC_jll
using LibXSPEC_Relxill_jll

using SpectralFitting
using SpectralFitting:
    @xspecmodel,
    @wrap_xspec_model_ccall,
    AbstractSpectralModel,
    Additive,
    FitParam,
    parameter_type,
    FreeParameters

@xspecmodel (:lmodrelline, libXSPEC_relxill) struct XS_Relline{T,F} <: AbstractSpectralModel{T,Additive}
    "Normalisation."
    K::T
    "Rest frame emission energy."
    lineE::T
    "Photon index 1."
    index1::T
    "Photon index 2."
    index2::T
    "Break radius between the photon indices."
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
    lineE = FitParam(6.4),
    index1 = FitParam(3.0),
    index2 = FitParam(3.0),
    r_break = FitParam(15.0),
    a = FitParam(0.998, upper_limit = 0.998),
    θ_obs = FitParam(30.0),
    inner_r = FitParam(-1.0),
    outer_r = FitParam(400.0),
    z = FitParam(0.0),
    limb = FitParam(0.0),
)
    XS_Relline{typeof(K),FreeParameters{(:K, :a)}}(
        K,
        lineE,
        index1,
        index2,
        r_break,
        a,
        θ_obs,
        inner_r,
        outer_r,
        z,
        limb,
    )
end

export XS_Relline

end # module Relxill
