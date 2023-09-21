module SpectralFittingExtras

@enum ModelType additive multiplicative convolutional

function _to_type_symbol(s::ModelType)
    if s == additive
        :Additive
    elseif s == multiplicative
        :Multiplicative
    else
        :Convolutional
    end
end

sanitize_name(name) = filter(!isnumeric, name)

function ModelType(s::AbstractString)
    if s == "add"
        additive
    elseif s == "mul"
        multiplicative
    elseif s == "con"
        convolutional
    else
        error("Unknown model type $(s)")
    end
end

struct ModelInfo
    name::String
    type::ModelType
    symbol::String
end

struct ModelParameter
    name::String
    value::Float64
    lower_limit::Float64
    upper_limit::Float64
end

function translate_model_dat(path::AbstractString, names)
    data = open(path) do file
        String(read(file))
    end
    model_data = split(data, "\n\n")
    map(zip(names, model_data)) do (name, md)
        _translate_model_dat(name, split(md, "\n"))
    end
end

function process_parameter(line)
    # so that split works okay replace
    items = split(replace(line, "\" \"" => "\"\""))
    ModelParameter(
        sanitize_name(items[1]),
        parse(Float64, items[3]),
        parse(Float64, items[4]),
        parse(Float64, items[6]),
    )
end

function process_model_info(line)
    items = split(line)
    ModelInfo(sanitize_name(items[1]), ModelType(items[6]), items[5])
end

function _translate_model_dat(name, lines)
    itt = Iterators.Stateful(filter(i -> length(i) > 0, lines))
    info = process_model_info(popfirst!(itt))
    params = map(process_parameter, itt)
    _build_model(name, info, params)
end

function _format_default(p::ModelParameter)
    "$(p.name) = FitParam($(p.value); lower_limit = $(p.lower_limit), upper_limit = $(p.upper_limit))"
end

function _build_model(name, info::ModelInfo, params::AbstractVector{<:ModelParameter})
    struct_params = join(["    $(p.name)::T" for p in params], "\n")

    default_params = join(map(_format_default, params), ", ")
    param_names = join([p.name for p in params], ", ")

    constructor = """
    function $(name)(; $(default_params))
        $(name){typeof($(params[1].name)),SpectralFitting.FreeParameters{()}}($param_names)
    end
    """

    """
    @xspecmodel :$(info.symbol) struct $(name){T,F} <: AbstractSpectralModel{T,$(_to_type_symbol(info.type))}
    $(struct_params)  
    end
    $constructor
    """
end

export translate_model_dat

end # module
