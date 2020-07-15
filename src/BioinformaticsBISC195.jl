module BioinformaticsBISC195

export normalizeDNA,
       composition,
       gc_content,
       complement,
       reverse_complement,
       parse_fasta

# # uncomment the following line if you intend to use BioSequences types
# using BioSequences

"""
    normalizeDNA(::AbstractString)

Ensures that a sequence only contains valid bases
(or `'N'` for unknown bases).
Returns a String.
"""
function normalizeDNA(seq)
    seq = uppercase(string(seq))
    for base in seq
        # note: `N` indicates an unknown base
        occursin(base, "AGCTN") || error("invalid base $base")
    end
    return seq # change to `return LongDNASeq(seq)` if you want to try to use BioSequences types
end

function composition(seq)
    seq = uppercase(seq)
    
    for b in seq
        !occursin(b, "ATGCN") && error("Invalid base, $b")
    end
    counts = Dict(b=>0 for b in collect("ATGCN"))
    for base in seq
        counts[base] += 1
    end
    return counts
end

function gc_content(seq)
    c = composition(seq)
    return (c['G'] + c['C']) / length(seq)
end

function complement(base::Char)
    comp = Dict('A'=>'T',
                'T'=>'A',
                'G'=>'C',
                'C'=>'G',
                'N'=>'N')
    return comp[uppercase(base)]
end

function complement(seq::AbstractString)
    return map(complement, seq)
end

reverse_complement(seq) = reverse(complement(seq))

mutable struct FastaRecord
    header
    sequence
end

# these are "accessor" functions - they're not strictly necessary
function header(fr::FastaRecord)
    return fr.header
end

function sequence(fr::FastaRecord)
    return fr.sequence
end

## Note: Simple functions like those above can be written with shortened syntax:
# header(fr::FastaRecord) = fr.header
# sequence(fr::FastaRecord) = fr.sequence

# sequence! updates the sequence field adds it to the end of the `sequence` field
function sequence!(fr::FastaRecord, seq::AbstractString)
    fr.sequence = seq
end

function parse_fasta(path)
    records = FastaRecord[]
    for line in eachline(path)
        if startswith(line, '>')
            # if the line is a header, we push! a new record with an empty sequence to the `records` vector
            header = line[2:end]
            push!(records, FastaRecord(header, ""))
        else
            # otherwise, we add the line onto the end of the sequence
            record = records[end]
            newseq = sequence(record) * line
            sequence!(record, newseq)
        end
    end
    return records
end

end # module BioinformaticsBISC195