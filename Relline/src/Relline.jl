module Relline


using LibXSPEC_jll
using LibXSPEC_Relline_jll

using SpectralFitting
using SpectralFitting:
    AbstractSpectralModel, Additive, FitParam, parameter_type, FreeParameters

struct XS_Relline{T,F} <: AbstractSpectralModel{Additive}
    K::FitParam{T}
    lineE::FitParam{T}
    index1::FitParam{T}
    index2::FitParam{T}
    r_break::FitParam{T}
    a::FitParam{T}
    θ_obs::FitParam{T}
    inner_r::FitParam{T}
    outer_r::FitParam{T}
    z::FitParam{T}
    limb::FitParam{T}
    function XS_Relline(;
        K = FitParam(1.0),
        lineE = FitParam(6.4),
        index1 = FitParam(3.0),
        index2 = FitParam(3.0),
        r_break = FitParam(15.0),
        a = FitParam(0.998),
        θ_obs = FitParam(30.0),
        inner_r = FitParam(-1.0),
        outer_r = FitParam(400.0),
        z = FitParam(0.0),
        limb = FitParam(0.0),
    )

        new{parameter_type(K),FreeParameters{(:K, :a)}}(
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
end

@inline function SpectralFitting.invoke!(flux, energy, ::Type{<:XS_Relline}, params...)
    _relline!(flux, energy, params)
end

# needs 32 bit floats ...
@inline function _relline!(flux::Vector{Float32}, energy::Vector{Float32}, params::Vector{Float32})
    ccall(
        (:tdrelline_, libXSPEC_relline),
        Cvoid,
        (Ptr{Float32}, Ref{Int32}, Ptr{Float32}, Ref{Int32}, Ptr{Float32}),
        energy,
        length(energy) - 1,
        params,
        2,
        flux,
    )
end

@inline function _relline!(flux, energy, params)
    f32_flux = Float32.(flux)
    f32_energy = Float32.(energy)
    f32_params = [Float32(i) for i in params]
    _relline!(f32_flux, f32_energy, f32_params)
    @. flux = f32_flux
end

SpectralFitting.implementation(::Type{<:XS_Relline}) = SpectralFitting.XSPECImplementation()

export XS_Relline

end # module Relline
