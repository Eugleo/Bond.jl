module Mass

export mass

# Source: https://www.unimod.org/masses.html
ATOMIC_WEIGHT = Dict(
    "H" => 1.007825035,
    "2H" => 2.014101779,
    "Li" => 7.016003,
    "B" => 11.0093055,
    "C" => 12,
    "13C" => 13.00335483,
    "N" => 14.003074,
    "15N" => 15.00010897,
    "O" => 15.99491463,
    "18O" => 17.9991603,
    "F" => 18.99840322,
    "Na" => 22.9897677,
    "Mg" => 23.9850423,
    "Al" => 26.9815386,
    "P" => 30.973762,
    "S" => 31.9720707,
    "Cl" => 34.96885272,
    "K" => 38.9637074,
    "Ca" => 39.9625906,
    "Cr" => 51.9405098,
    "Mn" => 54.9380471,
    "Fe" => 55.9349393,
    "Ni" => 57.9353462,
    "Co" => 58.9331976,
    "Cu" => 62.9295989,
    "Zn" => 63.9291448,
    "As" => 74.9215942,
    "Br" => 78.9183361,
    "Se" => 79.9165196,
    "Mo" => 97.9054073,
    "Ru" => 101.9043485,
    "Pd" => 105.903478,
    "Ag" => 106.905092,
    "Cd" => 113.903357,
    "I" => 126.904473,
    "Pt" => 194.964766,
    "Au" => 196.966543,
    "Hg" => 201.970617,
)

function get_count(count_str)
    isnothing(count_str) ? 1 : parse(Int, count_str)
end

function mass(formula::AbstractString)
    element = nothing
    count_str = nothing
    total = 0

    for ch in formula
        if isuppercase(ch)
            # Beginning of a new element

            if !isnothing(element)
                # Add the mass of the last element
                total += ATOMIC_WEIGHT[element] * get_count(count_str)
            end

            element = string(ch)
            count_str = nothing
        elseif islowercase(ch)
            # Continue building the element name
            element *= ch
        elseif isdigit(ch)
            # Start building the count, or continue with building the current one
            count_str = isnothing(count_str) ? string(ch) : count_str * ch
        end
    end
    total + ATOMIC_WEIGHT[element] * get_count(count_str)
end

end