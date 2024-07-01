module Reflionx

using SpectralFitting

struct ReflectionSpectrum{V}
    energy::V
    flux::V
end

function SpectralFitting.MultiLinearInterpolations.restructure(
    spec::ReflectionSpectrum,
    vs::AbstractVector,
)
    @views begin
        start = 1
        stop = start + length(spec.energy) - 1
        # TODO: energy grid probably always the same?
        start = stop + 1

        stop = start + length(spec.flux) - 1
        flux = vs[start:stop]
        ReflectionSpectrum(spec.energy[1:end], flux[1:end-1])
    end
end

struct TableData{N,T}
    params::NTuple{N,Vector{T}}
    grids::Array{ReflectionSpectrum{Vector{T}},N}
end

function Base.show(io::IO, ::MIME"text/plain", @nospecialize(data::TableData))
    print(
        io,
        """
        Reflionx Table Data (N=$(length(data.grids))):
          . Γ      $(extrema(data.params[1])) (N=$(length(data.params[1])))
          . logξ   $(extrema(data.params[2])) (N=$(length(data.params[2])))
          . E_cut  $(extrema(data.params[3])) (N=$(length(data.params[3])))
        """,
    )
end

function _read_reflionx_data()
    data = get_model_data(ReflionxTable)
    _Γ::Vector{Float64} = data[1]["Γ"]
    _logξ::Vector{Float64} = data[1]["logξ"]
    _E_cut::Vector{Float64} = data[1]["E_cut"]
    grid::Array{ReflectionSpectrum{Vector{Float64}},3} = data[1]["grid"]
    TableData((_Γ, _logξ, _E_cut), grid)
end

struct ReflionxTable{D,T} <: AbstractTableModel{T,Additive}
    table::D
    K::T
    "Photon index (e^-Γ)"
    Γ::T
    "Ionization index"
    logξ::T
    "High energy cutoff energy"
    E_cut::T
end

function SpectralFitting.invoke!(output, domain, model::ReflionxTable)
    spec = SpectralFitting.MultiLinearInterpolations.interpolate!(
        model.table.interp,
        model.table.data.params,
        model.table.data.grids,
        (model.Γ, model.logξ, model.E_cut),
    )
    SpectralFitting.rebin_if_different_domains!(output, domain, spec.energy, spec.flux)
end

function Base.copy(m::ReflionxTable)
    # create a new mutli-linear interpolator cache
    interp = SpectralFitting.MultilinearInterpolator{3}(m.table.data.grids)
    typeof(m)(
        (; interp = interp, data = m.table.data),
        (copy(getproperty(m, f)) for f in fieldnames(typeof(m))[2:end])...,
    )
end

function ReflionxTable(
    data::TableData;
    K = FitParam(1.0),
    Γ = FitParam(2.0),
    logξ = FitParam(2.0),
    E_cut = FitParam(125.0),
)
    interp = SpectralFitting.MultilinearInterpolator{3}(data.grids)
    ReflionxTable((; interp = interp, data = data), K, Γ, logξ, E_cut)
end

export ReflionxTable

include("parsing.jl")

end # module Reflionx
