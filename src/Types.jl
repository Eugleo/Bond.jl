module Types

using FASTX

export Protein, load_fasta

struct Protein
    sequence::String
end

function load_fasta(path)
    open(FASTA.Reader, path) do reader
        for record in reader
            close(reader)
            return Protein(FASTA.sequence(record))
        end
    end
end

end
