module Warmabs

using LibXSPEC_Warmabs_jll
using LibXSPEC_Warmabs_jll: libXSPEC_fphotems

using SpectralFitting
using SpectralFitting:
    AbstractSpectralModel, Multiplicative, FitParam, FreeParameters

@xspecmodel type = Float32 struct XS_WarmAbsorber{T,F} <:
                                  AbstractSpectralModel{T,Multiplicative}
    column::T
    rlogxi::T
    Cabund::T
    Nabund::T
    Oabund::T
    Fabund::T
    Neabund::T
    Naabund::T
    Mgabund::T
    Alabund::T
    Siabund::T
    Pabund::T
    Sabund::T
    Clabund::T
    Arabund::T
    Kabund::T
    Caabund::T
    Scabund::T
    Tiabund::T
    Vabund::T
    Crabund::T
    Mnabund::T
    Feabund::T
    Coabund::T
    Niabund::T
    Cuabund::T
    Znabund::T
    write_outfile::T
    outfile_idx::T
    vturb::T
    redshift::T
end
function XS_WarmAbsorber(;
    column = FitParam(0.0; lower_limit = -3.0, upper_limit = 2.0),
    rlogxi = FitParam(0.0; lower_limit = -4.0, upper_limit = 5.0),
    Cabund = FitParam(1.0; lower_limit = 0.0, upper_limit = 1000.0),
    Nabund = FitParam(1.0; lower_limit = 0.0, upper_limit = 1000.0),
    Oabund = FitParam(1.0; lower_limit = 0.0, upper_limit = 1000.0),
    Fabund = FitParam(1.0; lower_limit = 0.0, upper_limit = 1000.0),
    Neabund = FitParam(1.0; lower_limit = 0.0, upper_limit = 1000.0),
    Naabund = FitParam(1.0; lower_limit = 0.0, upper_limit = 1000.0),
    Mgabund = FitParam(1.0; lower_limit = 0.0, upper_limit = 1000.0),
    Alabund = FitParam(1.0; lower_limit = 0.0, upper_limit = 1000.0),
    Siabund = FitParam(1.0; lower_limit = 0.0, upper_limit = 1000.0),
    Pabund = FitParam(1.0; lower_limit = 0.0, upper_limit = 1000.0),
    Sabund = FitParam(1.0; lower_limit = 0.0, upper_limit = 1000.0),
    Clabund = FitParam(1.0; lower_limit = 0.0, upper_limit = 1000.0),
    Arabund = FitParam(1.0; lower_limit = 0.0, upper_limit = 1000.0),
    Kabund = FitParam(1.0; lower_limit = 0.0, upper_limit = 1000.0),
    Caabund = FitParam(1.0; lower_limit = 0.0, upper_limit = 1000.0),
    Scabund = FitParam(1.0; lower_limit = 0.0, upper_limit = 1000.0),
    Tiabund = FitParam(1.0; lower_limit = 0.0, upper_limit = 1000.0),
    Vabund = FitParam(1.0; lower_limit = 0.0, upper_limit = 1000.0),
    Crabund = FitParam(1.0; lower_limit = 0.0, upper_limit = 1000.0),
    Mnabund = FitParam(1.0; lower_limit = 0.0, upper_limit = 1000.0),
    Feabund = FitParam(1.0; lower_limit = 0.0, upper_limit = 1000.0),
    Coabund = FitParam(1.0; lower_limit = 0.0, upper_limit = 1000.0),
    Niabund = FitParam(1.0; lower_limit = 0.0, upper_limit = 1000.0),
    Cuabund = FitParam(1.0; lower_limit = 0.0, upper_limit = 1000.0),
    Znabund = FitParam(1.0; lower_limit = 0.0, upper_limit = 1000.0),
    write_outfile = FitParam(0.0; lower_limit = 0.0, upper_limit = 1.0),
    outfile_idx = FitParam(0.0; lower_limit = 0.0, upper_limit = 10000.0),
    vturb = FitParam(0.0; lower_limit = 0.0, upper_limit = 10000.0),
    Redshift = FitParam(0.0; lower_limit = 0.0, upper_limit = 10.0),
)
    XS_WarmAbsorber{typeof(column),SpectralFitting.FreeParameters{()}}(
        column,
        rlogxi,
        Cabund,
        Nabund,
        Oabund,
        Fabund,
        Neabund,
        Naabund,
        Mgabund,
        Alabund,
        Siabund,
        Pabund,
        Sabund,
        Clabund,
        Arabund,
        Kabund,
        Caabund,
        Scabund,
        Tiabund,
        Vabund,
        Crabund,
        Mnabund,
        Feabund,
        Coabund,
        Niabund,
        Cuabund,
        Znabund,
        write_outfile,
        outfile_idx,
        vturb,
        Redshift,
    )
end
function SpectralFitting._unsafe_ffi_invoke!(
    output::Vector{Float32},
    error_vec::Vector{Float32},
    input::Vector{Float32},
    params::Vector{Float32},
    ::Type{<:XS_WarmAbsorber},
)
    ccall(
        (:fwarmabs_, libXSPEC_fphotems),
        Cvoid,
        (Ptr{Float32}, Ref{Int32}, Ptr{Float32}, Ref{Int32}, Ptr{Float32}, Ptr{Float32}),
        input,
        length(input) - 1,
        params,
        1,
        output,
        error_vec,
    )
end

const MODEL_FILES = (
    "atdbwarmabs.fits",
    "xo01_detah2.fits",
    "xo01_detah3.fits",
    "xo01_detah4.fits",
    "xo01_detahl.fits",
    "xo01_detail.fits",
    "xo01_detal2.fits",
    "xo01_detal3.fits",
    "xo01_detal4.fits",
)

# register all model data
for model in (XS_WarmAbsorber,)
    SpectralFitting.register_model_data(
        model,
        (
            ("warmabs/$i", joinpath(LibXSPEC_Warmabs_jll.artifact_dir, "data", i)) for
            i in MODEL_FILES
        )...,
    )
end

export XS_WarmAbsorber

end # module Warmabs
