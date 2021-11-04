module Precursor

function go!(
    objects,
    index,
    remaining_mass,
    remaining_segments,
    start_segment,
    chosen,
    solutions,
)
    if remaining_mass ≈ 0
        push!(solutions, objects[chosen])
    elseif remaining_mass < 0 || index == length(objects) + 1
        solutions
    else
        if remaining_segments > 0 && !start_segment
            for beginning = (index+1):length(objects)
                go!(
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

function subset_sum(objects::Objects, target_mass, segments) where Objects
    solutions = Objects[]
    chosen = UInt16[]
    sizehint!(chosen, length(objects))
    for beginning = 1:length(objects)
        go!(objects, beginning, target_mass, segments - 1, true, chosen, solutions)
    end
    solutions
end



end