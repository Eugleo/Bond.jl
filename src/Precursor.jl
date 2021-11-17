module Precursor

using DataFrames

include("Mass.jl")
include("Types.jl")

using .Mass, .Types

function step!(objects, masses, transitions; from, to, bond)
    if to in keys(transitions)
        next_parents = transitions[to]
        # if !isnothing(next_parents) && !(from in next_parents)
            push!(transitions[to], from)
        # end
    else
        transitions[to] = from[1] == 0 ? nothing : [from]
    end
    compute_step!(objects, masses, transitions; from=to, bond)
end

function compute_step!(objects, weights, transitions; bond, from)
    index, target, jumps_left, cysteines = from
    is_first = index == 0
    if target <= 0 || jumps_left < 0 || (!is_first && index >= length(objects))
        return transitions
    end

    can_jump = is_first || (jumps_left > 0  && cysteines > 0)
    max_jump_size = can_jump ? max(length(objects) - index - 1, 0) : 0
    for skipped = 0:max_jump_size
        jumped_with_bond = !is_first && skipped > 0
        next_object = (index + 1) + skipped
        new_target = target - weights[next_object] + jumped_with_bond * bond
        if new_target >= 0
            to = (
                next_object,
                new_target,
                jumps_left - jumped_with_bond,
                cysteines + objects[next_object].cysteines - jumped_with_bond * 2,
            )
            step!(objects, weights, transitions; from, to, bond)
        end
    end

    transitions
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

Step = Tuple{Int,Int,Int,Int}

function approximate(number, mult)
    ceil(Int, mult * number)
end

function subset_sum(objects, target, max_jumps, sensitivity = 1000)
    println("HeyHoa")
    transitions = Dict{Step,Union{Vector{Step},Nothing}}()
    bond = approximate(mass("H2O") - mass("H2"), sensitivity)
    masses = [approximate(d.mass, sensitivity) for d in objects]
    from = (0, target - approximate(mass("H2O"), sensitivity), max_jumps, 0)
    compute_step!(objects, masses, transitions; bond, from)

    transitions # walk(table, objects)
end

# objects40 = map(
#     (m, cs) -> (mass = m, cysteines = cs),
#     rand(Random.seed!(42), 10:50000, 40),
#     rand(Random.seed!(42), 0:2, 40),
# )

end
