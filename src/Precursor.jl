module Precursor

function go!(
    counter,
    objects,
    index,
    remaining_mass,
    remaining_segments,
    start_segment,
    chosen,
    solutions,
)
    #counter[index, max(remaining_mass + 1, 1), remaining_segments + 1] += 1

    if remaining_mass ≈ 0
        push!(solutions, objects[chosen])
    elseif remaining_mass < 0 || index == length(objects) + 1
        solutions
    else
        if remaining_segments > 0 && !start_segment
            for beginning = (index+1):length(objects)
                go!(
                    counter,
                    objects,
                    beginning,
                    remaining_mass,
                    remaining_segments,
                    true,
                    chosen,
                    solutions,
                )
            end
        end

        go!(
            counter,
            objects,
            index + 1,
            remaining_mass - objects[index].mass,
            remaining_segments - start_segment,
            false,
            push!(chosen, index),
            solutions,
        )
        pop!(chosen)
        solutions
    end
end

function perform_step!(objects, table; from, to)
    index, mass, segments = to
    next_step = (index, mass + 1, segments + 1)

    if isassigned(table, index, mass + 1, segments + 1)
        next_parents = table[next_step...]
        if !isnothing(next_parents) && !(from in next_parents)
            push!(table[next_step...], from)
        end
    else
        table[next_step...] = isnothing(from) ? nothing : [from]
        step!(objects, index, mass, segments, table)
    end
end

function step!(objects, index, remaining_mass, segment_budget, table)
    if remaining_mass <= 0 || segment_budget < 0 || index >= length(objects)
        return table
    end

    maxskips = max(length(objects) - index - 1, 0)
    skips = 0:(segment_budget > 0 ? maxskips : 0)
    for skip in skips
        next = index + skip + 1
        next_mass = remaining_mass - objects[next].mass

        if next_mass >= 0
            from = (index, segment_budget)
            to = (next, next_mass, segment_budget - (skip > 0))
            perform_step!(objects, table; from, to)
        end
    end

    table
end

function walk(table, objects)
    solution = []
    indices, _, segments = size(table)
    for I in CartesianIndices((1:indices, 1, 1:segments))
        index, mass, segment = Tuple(I)
        if isassigned(table, index, mass, segment)
            walk!(table, objects, table[I], objects[index].mass, [index], solution)
        end
    end
    solution
end

function walk!(table, objects, variants, remaining_mass, path, solution)
    if isnothing(variants)
        # We reached the end of the path
        push!(solution, objects[path])
    else
        for (index, segment) in variants
            parents = table[index, remaining_mass+1, segment+1]
            previous_mass = remaining_mass + objects[index].mass
            push!(path, index)
            walk!(table, objects, parents, previous_mass, path, solution)
            pop!(path)
        end
        solution
    end
end

function subset_sum(objects::Objects, target_mass, segments) where {Objects}
    dims = (length(objects), target_mass + 1, segments + 1)
    table = Array{Union{Vector{Tuple{Int,Int}},Nothing}}(undef, dims)
    skips = 0:(segments > 0 ? max(length(objects) - 1, 0) : 0)
    for skip in skips
        next = skip + 1
        next_mass = target_mass - objects[next].mass

        if next_mass >= 0
            to = (next, next_mass, segments - 1)
            perform_step!(objects, table; from = nothing, to)
        end
    end

    walk(table, objects)
end

end
