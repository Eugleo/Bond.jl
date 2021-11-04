module Precursor

function match_precursor!(digestides, index, target_mass, segments, solutions)
    if target_mass ≈ 0
        push!(solutions, digestides[segments])
    elseif target_mass < 0 || index == length(digestides) + 1
        solutions
    else
        match_precursor!(digestides, index + 1, target_mass, segments, solutions)
        match_precursor!(
            digestides,
            index + 1,
            target_mass - digestides[index].mass,
            push!(segments, index),
            solutions,
        )
        pop!(segments)
        solutions
    end
end

function subset_sum(objects, target_mass)
    function go!(index, remaining_mass, chosen, solutions)
        if remaining_mass ≈ 0
            push!(solutions, objects[chosen])
        elseif remaining_mass < 0 || index == length(objects) + 1
            solutions
        else
            go!(index + 1, remaining_mass, chosen, solutions)
            go!(
                index + 1,
                remaining_mass - objects[index].mass,
                push!(chosen, index),
                solutions,
            )
            pop!(chosen)
            solutions
        end
    end

    chosen = []
    sizehint!(chosen, length(objects))
    go!(1, target_mass, chosen, [])
end



end