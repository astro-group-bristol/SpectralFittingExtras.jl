struct _RunData
    path::String
    Γ::Float64
    logξ::Float64
    E_cut::Float64
end

function _RunData(name::String)
    s = basename(name)
    m = match(r"gamma_(?<gamma>[\d\.]+)_logxi_(?<logxi>[\d\.]+)_Ecut_(?<ecut>[\d\.]+)", s)
    _RunData(
        name,
        parse(Float64, m["gamma"]),
        parse(Float64, m["logxi"]),
        parse(Float64, m["ecut"]),
    )
end

function _read_run_results(path::String)
    energies = Float64[]
    fluxes = Float64[]
    for line in eachline(path)
        values = split(line, " "; keepempty = false)
        push!(energies, parse(Float64, values[1]))
        push!(fluxes, parse(Float64, values[end]))
    end
    ReflectionSpectrum(energies, fluxes)
end

"""
    parse_run(directory_root::String; verbose = true)
    
Parse the output of a simulation run into a table model grid.
"""
function parse_run(directory_root::String; verbose = true)
    runs = _RunData[]
    for (root, _, files) in walkdir(directory_root)
        for f in files
            if endswith(f, "_output.dat")
                N_files = length(runs) + 1
                verbose && print("Gathering $N_files output files...\r")
                push!(runs, _RunData(joinpath(root, f)))
            end
        end
    end

    N_files = length(runs)
    verbose && println("Gathered total $N_files output files. Parsing...")

    # work out the axes for each variable
    gamma = sort(collect(Set([i.Γ for i in runs])))
    logxi = sort(collect(Set([i.logξ for i in runs])))
    ecut = sort(collect(Set([i.E_cut for i in runs])))

    # this is a really lazy way of making sure the data is stored in the right order
    function _arranger(Γ, logξ, E_cut)::String
        for r in runs
            if r.Γ == Γ && r.logξ == logξ && r.E_cut == E_cut
                return r.path
            end
        end
        error("No combination: $((Γ, logξ, E_cut))")
    end

    paths = [_arranger(Γ, logξ, E_cut) for Γ in gamma, logξ in logxi, E_cut in ecut]
    # now we read all the files and get the right values 
    TableData((gamma, logxi, ecut), _read_run_results.(paths))
end
