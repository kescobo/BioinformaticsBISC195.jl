module BioinformaticsBISC195

export FastaRecord,
       header,
       sequence,
       sequence!,
       normalizeDNA,
       composition,
       gc_content,
       complement,
       reverse_complement,
       parse_fasta,
       uniquekmers,
       kmer_distance

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
        occursin(base, "AGCTNYRWMKSHVDB") || error("invalid base $base")
    end
    return seq # change to `return LongDNASeq(seq)` if you want to try to use BioSequences types
end

function composition(seq)
    seq = normalizeDNA(seq)
    counts = Dict(b=>0 for b in collect("ATGCN"))
    for base in seq
        haskey(counts, base) ? counts[base] += 1 : counts['N'] +=1
    end
    return counts
end

function gc_content(seq)
    c = composition(seq)
    return (c['G'] + c['C']) / sum(c[b] for b in "ATGC")
end

function complement(base::Char)
    base = uppercase(base)
    comp = Dict('A'=>'T',
                'T'=>'A',
                'G'=>'C',
                'C'=>'G')
    haskey(comp, base) ? comp[base] : "blah"
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

Base.length(fr::FastaRecord) = length(sequence(fr))
gc_content(fr::FastaRecord) = gc_content(sequence(fr))

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

function uniquekmers(seq::AbstractString, k::Int)
    kmers = String[]
    i = 1
    e = lastindex(seq)-k+1
    while i < e
        kmer = seq[i:i+k-1]
        if !all(c-> occursin(c, "ATGC"), kmer)
            i += k
        else
            push!(kmers, kmer)
            i +=1
        end
    end
    return Set(kmers)
end

uniquekmers(record::FastaRecord, k::Int) = uniquekmers(sequence(record), k)

# function kmer_distance(k1, k2)
#     u = union(k1, k2)
#     i = intersect(k1,k2)

#     return 1 - length(i) / length(u)
# end

end # module BioinformaticsBISC195