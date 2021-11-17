module Types

using FASTX
include("Mass.jl")
using .Mass: mass

export load_protein, Peptide, digest, trypsin

struct Peptide
    sequence::String
    mass::Float32 # MAYBE: Put Float64 here
    cysteines::Int

    function Peptide(sequence)
        ends = mass("H2") + mass("O")
        mass = sum(amino_acid -> AMINO_ACID_MASS[amino_acid], sequence) + ends
        new(sequence, mass, Base.count(aa -> aa == 'C', sequence))
    end
end

# Source: pyteomics' mass module
const AMINO_ACID_MASS = Dict(
    'G' => 57.02146,
    'A' => 71.03711,
    'S' => 87.03203,
    'P' => 97.05276,
    'V' => 99.06841,
    'T' => 101.04768,
    'C' => 103.00919,
    'L' => 113.08406,
    'I' => 113.08406,
    'J' => 113.08406,
    'N' => 114.04293,
    'D' => 115.02694,
    'Q' => 128.05858,
    'K' => 128.09496,
    'E' => 129.04259,
    'M' => 131.04049,
    'H' => 137.05891,
    'F' => 147.06841,
    'U' => 150.95364,
    'R' => 156.10111,
    'Y' => 163.06333,
    'W' => 186.07931,
    'O' => 237.14773,
)

function load_protein(path)
    open(FASTA.Reader, path) do reader
        for record in reader
            close(reader)
            return Peptide(FASTA.sequence(record))
        end
    end
end

trypsin = Regex("(?<=[KR])(?!P)")

function digest(protein, enzyme)
    [Peptide(seq) for seq in split(protein.sequence, enzyme)]
end

function count(biopolymer, amino_acid)
    Base.count(aa -> aa == amino_acid, biopolymer.sequence)
end

end
