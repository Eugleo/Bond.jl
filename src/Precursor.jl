module Precursor

using DataFrames
using OrderedCollections

include("Mass.jl")
include("Types.jl")

using .Mass, .Types

function try_stepping!(objects, table; from, to, bond)
    if to in keys(table)
        next_parents = table[to]
        if !isnothing(next_parents) && !(from in next_parents)
            push!(table[to], from)
        end
    else
        table[to] = isnothing(from) ? nothing : [from]
        compute_step!(objects, table; to, bond)
    end
end

function compute_step!(objects, table; to, bond)
    index, remaining_mass, segment_budget, cysteines = to
    if remaining_mass <= 0 || segment_budget < 0 || index >= length(objects)
        return table
    end

    maxskips = max(length(objects) - index - 1, 0)
    for skip = 0:((segment_budget > 0 && cysteines > 0) ? maxskips : 0)
        next = (index + 1) + skip
        next_mass = remaining_mass - objects[next].approx_mass + (skip > 0) * bond
        next_cysteines = cysteines + objects[next].cysteines - (skip > 0) * 2

        if next_mass >= 0
            next = (next, next_mass, segment_budget - (skip > 0), next_cysteines)
            try_stepping!(objects, table; from = to, to = next, bond)
        end
    end

    table
end

function show_solutions(table, objects, segments, alkylation_mass, mods)
    result = DataFrame()

    for (segments_left, free_cysteines, modname, path) in
        walk(table, objects, segments, alkylation_mass, mods)
        push!(
            result,
            Dict(
                "segments" => segments - segments_left,
                "total_cysteines" => sum(i -> i.cysteines, objects[path]),
                "free_cysteines" => free_cysteines,
                "selected" => join([(k in path) ? "X" : "_" for k = 1:length(objects)]),
                "modifications" => modname,
            ),
            cols = :union,
        )
    end

    result
end

function compute_modification_steps(masses, bounds)
    [
        (ns, sum(ns)) for
        ns in Base.product([(0:bounds[mod]) .* m for (mod, m) in masses]...)
    ]
end

# TODO: Implementovat modifikace — prostě se bude procházet tabulka z jiné vrstvy než mass = 0
function walk(table, objects, segments, alkylation_mass, modification_masses)
    solutions = []

    mod_bounds = Dict()
    for mod in keys(modification_masses)
        mod_bounds[mod] = sum(o -> modifications(o, mod), objects)
    end

    modifiers = compute_modification_steps(modification_masses, mod_bounds)


    for final_index = 1:length(objects),
        segments_left = 0:segments,
        free_cysteines = 0:sum(o -> o.cysteines, objects),
        (modname, mod) in modifiers

        final_step = (
            final_index,
            alkylation_mass * free_cysteines + mod,
            segments_left,
            free_cysteines,
        )
        if final_step in keys(table)
            if segments_left >= 6
                println(final_step)
            end
            sols = walk!(table, table[final_step], [final_index], [])
            append!(
                solutions,
                [(segments_left, free_cysteines, modname, s) for s in sols],
            )
        end
    end
    solutions
end

function walk!(table, parents, path, solutions)
    if isnothing(parents)
        # We reached the end of the path
        push!(solutions, copy(path))
    else
        for parent in parents
            push!(path, parent[1])
            walk!(table, table[parent], path, solutions)
            pop!(path)
        end

        solutions
    end
end

function subset_sum(digestides, target_mass, segments, sensitivity = 1000)
    table = Dict{Tuple{Int,Int,Int,Int},Union{Vector{Tuple{Int,Int,Int,Int}},Nothing}}()
    skips = 0:(segments > 0 ? max(length(digestides) - 1, 0) : 0)

    h2o = approximate(mass("H2O"), sensitivity)
    bond = -approximate(mass("H2"), sensitivity)

    for skip in skips
        next = skip + 1
        next_mass = target_mass - digestides[next].approx_mass - h2o

        if next_mass >= 0
            to = (next, next_mass, segments - 1, digestides[next].cysteines)
            try_stepping!(digestides, table; from = nothing, to, bond + h2o)
        end
    end

    table # walk(table, objects)
end

# objects40 = map(
#     (m, cs) -> (mass = m, cysteines = cs),
#     rand(Random.seed!(42), 10:50000, 40),
#     rand(Random.seed!(42), 0:2, 40),
# )

end
