module Precursor

using DataFrames
using OrderedCollections

include("Mass.jl")
using .Mass: mass

const H = ceil(Int, mass("H") * 100)
const BOND = ceil(Int, -mass("H2") * 100)
const H2O = ceil(Int, mass("H2O") * 100)

function try_stepping!(objects, table; from, to)
    if to in keys(table)
        next_parents = table[to]
        if !isnothing(next_parents) && !(from in next_parents)
            push!(table[to], from)
        end
    else
        table[to] = isnothing(from) ? nothing : [from]
        compute_step!(objects, table; to)
    end
end

function compute_step!(objects, table; to)
    index, remaining_mass, segment_budget, cysteines = to
    if remaining_mass <= 0 || segment_budget < 0 || index >= length(objects)
        return table
    end

    maxskips = max(length(objects) - index - 1, 0)
    for skip = 0:((segment_budget > 0 && cysteines > 0) ? maxskips : 0)
        next = (index + 1) + skip
        next_mass = remaining_mass - objects[next].mass + (skip > 0) * (H2O + BOND)
        next_cysteines = cysteines + objects[next].cysteines - (skip > 0) * 2

        if next_mass >= 0
            next = (next, next_mass, segment_budget - (skip > 0), next_cysteines)
            try_stepping!(objects, table; from = to, to = next)
        end
    end

    table
end

function show_solutions(table, objects, segments, alkylation_mass)
    result = DataFrame()

    for (segments_left, free_cysteines, path) in walk(table, objects, segments, alkylation_mass)
        push!(
            result,
            Dict(
                "segments" => segments - segments_left,
                "total_cysteines" => sum(i -> i.cysteines, objects[path]),
                "free_cysteines" => free_cysteines,
                "selected" => join([(k in path) ? "X" : "_" for k in 1:length(objects)])
            ),
            cols = :union
        )
    end

    result
end

# TODO: Implementovat modifikace — prostě se bude procházet tabulka z jiné vrstvy než mass = 0
# Bude to fungovat dobře pro nízký-ish počet modfifikací
function walk(table, objects, segments, alkylation_mass)
    solutions = []
    for final_index = 1:length(objects),
        segments_left = 0:segments,
        free_cysteines = 0:sum(o -> o.cysteines, objects)

        final_step = (final_index, alkylation_mass * free_cysteines, segments_left, free_cysteines)
        if final_step in keys(table)
            if segments_left >= 6
                println(final_step)
            end
            sols = walk!(table, table[final_step], [final_index], [])
            append!(solutions, [(segments_left, free_cysteines, s) for s in sols])
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

function subset_sum(objects::Objects, target_mass, segments) where {Objects}
    table = Dict{Tuple{Int,Int,Int,Int},Union{Vector{Tuple{Int,Int,Int,Int}},Nothing}}()
    skips = 0:(segments > 0 ? max(length(objects) - 1, 0) : 0)
    for skip in skips
        next = skip + 1
        next_mass = target_mass - objects[next].mass - H2O

        if next_mass >= 0
            to = (next, next_mass, segments - 1, objects[next].cysteines)
            try_stepping!(objects, table; from = nothing, to)
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
