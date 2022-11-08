using SharedArrays,Distributed,Printf, LinearAlgebra, Statistics
using NLopt
import DCAUtils: read_fasta_alignment, remove_duplicate_sequences, compute_weights, add_pseudocount, compute_weighted_frequencies
using LoopVectorization
using Base.Threads: @threads, nthreads
using FastaIO
using GZip





function n2l(n)
    alphabetN = [ 1,   2,   3,   4,   5,   6,   7,   8,   9,   10,  11,  12,  13, 14,  15,  16,  17,  18,  19,  20, 21]
    alphabetL=  ["A", "C", "D", "E", "F", "G", "H", "I", "K", "L", "M", "N", "P", "Q", "R", "S", "T", "V", "W", "Y", "-"]
    dd =  Dict(alphabetN[i] => alphabetL[i] for i =1:21)
    return dd[n]
end




function ReadFasta(filename::AbstractString,max_gap_fraction::Real, theta::Any, remove_dups::Bool)

    Z = read_fasta_alignment(filename, max_gap_fraction)
    if remove_dups
        Z, _ = remove_duplicate_sequences(Z)
    end

    N, M = size(Z)
    q = round(Int,maximum(Z))

    q > 32 && error("parameter q=$q is too big (max 31 is allowed)")
    W , Meff = compute_weights(Z,q,theta)
    rmul!(W, 1.0/Meff)
    Zint=round.(Int,Z)
    return W, Zint,N,M,q
end


function write_fasta_data(filename::String, data)
    println("Writing $(length(data)) entries to file $filename")

    # note: writefasta can take a writable IO descriptor as
    #       argument instead of a filename; also, it has
    #       an optional mode argument which defaults to "w"
    #       (use "a" for appending to an exisitng file)
    writefasta(filename, data)

#     println("This is the content of the file:")
#     gzopen(filename) do f
#         println(read(f, String))
#     end
end




function compute_APC(J::Array{Float64,4},N,q)
    FN = fill(0.0, N,N)
    for i=1:N-1
        for j=i+1:N
            FN[i,j] = norm(J[1:q-1,1:q-1,i,j],2)
            FN[j,i] =FN[i,j]
        end
    end
    FN=correct_APC(FN)
    return FN
end




function compute_ranking(S::Matrix{Float64}, min_separation::Int = 5)
    N = size(S, 1)
    R = Array{Tuple{Int,Int,Float64}}(undef, div((N-min_separation)*(N-min_separation+1), 2))
    counter = 0
    for i = 1:N-min_separation, j = i+min_separation:N
        counter += 1
        R[counter] = (i, j, S[j,i])
    end

    sort!(R, by=x->x[3], rev=true)
    return R
end

function expandPi(Pi_true,N)
    # N, _ = size(Pi_true)
    # @show N
    c = zeros(N,21)
    c[:,1:20] .= transpose(reshape(Pi_true, (20,N)))
    c[:,21].= (1 .-sum(c, dims=2))[:,1]
    c[:,21]
    return c
end
# c = expandPi(Pi_true)

function expandP(Pi_true, Pij_true, N)
    c = expandPi(Pi_true,N)
    q=21
    # N,q = size(c)
    ce = reshape(repeat(c, N), (N,N,21))
    a = permutedims(reshape(Pij_true, (20,N,20,N)),(2,1,4,3))
    b = zeros((N,21,N,21))
    b[:,1:20, :, 1:20] .= a
    # b[:,21,:,21]
    b[:,21,:,21]
    # a = reshape(Pij_true, (53,20,53,20))
    #b[:,21,:,:] .= (ce .- sum(b, dims=(2))[:,1,:,:])
    b[:,21,:,:] .= (permutedims(ce,(2,1,3)) .- sum(b, dims=(2))[:,1,:,:])
    b[:,:,:,21] .= (permutedims(ce,(1,3,2)) .- sum(b, dims=(4))[:,:,:,1])
    b[:,21,:,21] .=0
    b[:,21,:,21] .= (1 .- sum(b, dims=(2,4))[:,1,:,1])
    b[3,:,5,:]
    return c, b
end

function corrCIJ(Pi_true, Pij_true,N)
    pie, pije = expandP(Pi_true, Pij_true,N)
    N,q = size(pie)
    a = reshape(permutedims(pije,(2,1,4,3)), (1113, 1113)) - vec(transpose(pie)) * transpose(vec( transpose(pie)))
    b = permutedims(reshape(a,(21,N,21,N)), (2,1,4,3))

    c = vec([])
    for i = 1:(N-1)
        for j =i+1:N
            if i!=j
                c = vcat(c, vec(b[i,:,j,:]))
            end
        end
    end

    return c
end




function parsemut(toparse)
    smut = split(toparse, "")
    ina = smut[1]
    outa = smut[end]
    site = parse(Int,join(smut[2:end-1]))
    return ina, outa, site
end




function pos_clean_dict(npzpath)

    y = npzread(npzpath)
    d = Dict()
    count=0
    for i =1:length(y)
        if y[i]
            count+=1
            d[i] = count
        else
            d[i] = -1
        end
    end
    return d
end

function processCSV(csv, npzpath)
    posdict = pos_clean_dict(npzpath)
    todelete = []
    for mut_id = 1:size(csv)[1]
        mut = csv[mut_id, :].mutant
        screenscore = csv[mut_id, :].screenscore
        ina, outa, site = parsemut(mut)
        if posdict[site] == -1
            push!(todelete, mut_id)
        elseif ismissing(screenscore)
            push!(todelete, mut_id)
        else
            csv[mut_id, :].mutant = ina*string(posdict[site])*outa
        end
    end
    delete!(csv, todelete)
    return csv
end


function profile(Pi_true, N, pseudocount)
    pie = expandPi(Pi_true,N)
    pip = (1-pseudocount).*pie .+(pseudocount/21)
    fields = log.(pip)
    proba = exp.(fields)
    proba .= proba./repeat(sum(proba, dims=2), 1, 21)
    return fields, proba
end
