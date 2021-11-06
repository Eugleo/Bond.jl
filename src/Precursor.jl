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

function subset_sum(objects::Objects, target_mass, segments) where {Objects}
    solutions = Objects[]
    counter = [] #zeros(length(objects) + 1, target_mass + 1, segments + 1)

    Threads.@threads for beginning = 1:length(objects)
        chosen = UInt16[]
        sizehint!(chosen, length(objects))
        go!(
            counter,
            objects,
            beginning,
            target_mass,
            segments - 1,
            true,
            chosen,
            solutions,
        )
    end
    solutions
end

function go!2(objects, index, remaining_mass, segment_budget, start_segment, solution)
    if remaining_mass < 0 || index == length(objects) + 1
        solution
    else
        if segment_budget > 0 && !start_segment
            for beginning = (index+1):length(objects)
                if isnothing(solution[beginning, remaining_mass+1, segment_budget+1])
                    solution[beginning, remaining_mass+1, segment_budget+1] =
                        [(index, segment_budget)]
                else
                    push!(
                        solution[beginning, remaining_mass+1, segment_budget+1],
                        (index, segment_budget),
                    )
                end
                go!2(objects, beginning, remaining_mass, segment_budget, true, solution)
            end
        end

        newmass = remaining_mass - objects[index].mass
        newsegments = segment_budget - start_segment

        if newmass >= 0
            if isnothing(solution[index+1, newmass+1, newsegments+1])
                solution[index+1, newmass+1, newsegments+1] = [(index, segment_budget)]
            else
                push!(
                    solution[index+1, newmass+1, newsegments+1],
                    (index, segment_budget),
                )
            end
        end

        go!2(objects, index + 1, newmass, newsegments, false, solution)

        solution
    end
end


function tried_action(table, objects, index, mass, segments)
    next_mass = mass - objects[index].mass
    if next_mass >= 0
        start_segments = segments - 1
        if start_segments > 0
            start_step = table[index+1, next_mass+1, start_segments+1]
            did_start = !isnothing(start_step) && (index, segments) in start_step
        else
            did_start = false
        end

        elongate_step = table[index+1, next_mass+1, segments+1]
        did_elongate = !isnothing(elongate_step) && (index, segments) in elongate_step

        (starting = did_start, elongating = did_elongate)
    else
        (starting = false, elongating = false)
    end
end

function start_new_segment!(objects, beginning, remaining_mass, segment_budget, table)
    go!3(objects, beginning, remaining_mass, segment_budget, true, table)
end

function elongate_segment!(objects, index, remaining_mass, segment_budget, table)
    go!3(objects, index, remaining_mass, segment_budget, false, table)
end

function add_as_source!(table, here_i, here_seg, next_i, next_mass, next_seg)
    next_step = (next_i, next_mass + 1, next_seg + 1)
    if isnothing(table[next_step...])
        table[next_step...] = []
    end

    println("Linking ", here_i * 5, " as a parent to ", next_i * 5, " while ", next_mass, " remains")

    here = (here_i, here_seg)
    if !(here in table[next_step...])
        push!(table[next_step...], here)
    end
end

function go!3(objects, index, remaining_mass, segment_budget, start_segment, solution)
    if remaining_mass <= 0 || (start_segment && segment_budget <= 0)
        return solution
    end

    # Tady bychom měli taky něco addnout
    if !start_segment && segment_budget > 0
        for beginning = (index+1):length(objects)
            add_as_source!(
                solution,
                index,
                segment_budget,
                beginning,
                remaining_mass,
                segment_budget,
            )

            start_new_segment!(
                objects,
                beginning,
                remaining_mass,
                segment_budget,
                solution,
            )
        end
    end

    next_mass = remaining_mass - objects[index].mass

    if next_mass >= 0 && index + 1 <= length(objects)
        next_segment_budget = segment_budget - start_segment

        tried = tried_action(solution, objects, index, remaining_mass, segment_budget)

        add_as_source!(
            solution,
            index,
            segment_budget,
            index + 1,
            next_mass,
            next_segment_budget,
        )

        if (start_segment && !tried.starting) || (!start_segment && !tried.elongating)
            elongate_segment!(
                objects,
                index + 1,
                next_mass,
                next_segment_budget,
                solution,
            )
        end
    end

    solution
end


function subset_sum2(objects::Objects, target_mass, segments, fun) where {Objects}
    table = Array{Union{Vector{Tuple{Int,Int}},Nothing}}(
        undef,
        (length(objects), target_mass + 1, segments + 1),
    )
    fill!(table, nothing)

    for beginning in eachindex(objects)
        fun(objects, beginning, target_mass, segments, true, table)
    end

    # walk(table, objects)
    table
end

function walk(table, objects)
    solution = []
    successes = table[:, 1, :]
    for I in CartesianIndices(successes)
        if !isnothing(successes[I])
            walk!(table, objects, successes[I], 0, [I[1]], solution)
        end
    end
    solution
end

function walk!(table, objects, starting_points, mass, path, solution)
    if isnothing(starting_points)
        # We reached the end of the path
        push!(solution, objects[path])
    else
        for (index, segment) in Set(starting_points)
            push!(path, index)
            walk!(
                table,
                objects,
                table[index, mass+objects[index].mass+1, segment+1],
                mass,
                path,
                solution,
            )
            pop!(path)
        end
        solution
    end
end

end
