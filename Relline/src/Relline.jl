module Relline

using LibXSPEC_jll
using LibXSPEC_Relline_jll
using LibXSPEC_Relline_jll: libXSPEC_relline

using SpectralFitting
using SpectralFitting:
    AbstractSpectralModel, Additive, FitParam, FreeParameters

@xspecmodel type = Float32 struct XS_Relline{T,F} <: AbstractSpectralModel{T,Additive}
    K::T
    lineE::T
    index1::T
    index2::T
    r_break::T
    a::T
    θ_obs::T
    inner_r::T
    outer_r::T
    z::T
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

function SpectralFitting._unsafe_ffi_invoke!(
    output::Vector{Float32},
    error_vec,
    input::Vector{Float32},
    params::Vector{Float32},
    ::Type{<:XS_Relline},
)
    ccall(
        (:tdrelline_, libXSPEC_relline),
        Cvoid,
        (Ptr{Float32}, Ref{Int32}, Ptr{Float32}, Ref{Int32}, Ptr{Float32}),
        input,
        length(input) - 1,
        params,
        2,
        output,
    )
end

export XS_Relline

function __init__()
    @warn (
        "This model is outdated (relline v0.4a), and has numerical problems at shallow and steep observer inclination. To use the most up-to-date `relline` model, please install `Relxill` with `relline` v0.5a. This can be done via the Julia package manager."
    )
end

end # module Relline
