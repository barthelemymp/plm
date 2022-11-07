using SharedArrays,Distributed,Printf, LinearAlgebra, Statistics
using NLopt
import DCAUtils: read_fasta_alignment, remove_duplicate_sequences, compute_weights, add_pseudocount, compute_weighted_frequencies
using LoopVectorization
using Base.Threads: @threads, nthreads
using Distributions
struct PlmAlg
    method::Symbol
    verbose::Bool
    epsconv::Float64
    maxit::Int
end

struct PlmOut
    pslike::Union{Vector{Float64},Float64}
    Jtensor::Array{Float64,4}
    htensor::Array{Float64,2}
    score::Array{Tuple{Int, Int, Float64},1}
end

struct PlmVar
    N::Int
    M::Int
    q::Int
    q2::Int
    lambdaJ::Float64
    lambdaH::Float64
    Z::SharedArray{Int,2}
    W::SharedArray{Float64,1}
    IdxZ::SharedArray{Int,2} #partial index computation for speed up energy calculation

    function PlmVar(N,M,q,q2,lambdaJ, lambdaH, Z,W)
        sZ = SharedArray{Int}(size(Z))
        sZ[:] = Z
        sW = SharedArray{Float64}(size(W))
        sW[:] = W

        IdxZ = Array{Int}(undef, N, M)
        q2=q*q
        for i in 1:M
            for j in 1:N
                IdxZ[j,i] = (j-1) * q2 + q * (Z[j,i] - 1)
            end
        end
        sIdxZ = SharedArray{Int}(size(IdxZ))
        sIdxZ[:] = IdxZ
        new(N,M,q,q2,lambdaJ, lambdaH, sZ, sW, sIdxZ)
    end
end


struct DecVar{N}
    fracdec::Float64
    fracmax::Float64
    blockdecimate::Bool
    dmask::SharedArray{Bool,N}
    function DecVar{N}(fracdec, fracmax, blockdecimate, dmask) where N
        sdmask = SharedArray{Bool}(size(dmask))
        sdmask[:] = dmask
        new(fracdec, fracmax, blockdecimate, sdmask)
    end
end





function optimfunwrapper(x::Vector, g::Vector, site, var)
    g === nothing && (g = zeros(Float64, length(x)))
    return PLsiteAndGrad!(x, g, site,  var)
end

function optimfunwrapper(x::Vector, g::Vector, var)
    g === nothing && (g = zeros(Float64, length(x)))
    return PLsiteAndGradSym!(x, g, var)
end

function inflate_matrix(J::Array{Float64,3},N)
    q,q,NN = size(J)

    @assert (N*(N-1))>>1 == NN

    Jt = zeros(q,q,N,N)
    ctr = 0
    for i in 1:N-1
        for j in i+1:N
            ctr += 1
            Jt[:,:,i,j] = J[:,:,ctr]
            Jt[:,:,j,i] = J[:,:,ctr]'
        end
    end
    return Jt
end

function correct_APC(S::Matrix)
    N = size(S, 1)
    Si = sum(S, dims=1)
    Sj = sum(S, dims=2)
    Sa = sum(S) * (1 - 1/N)

    S -= (Sj * Si) / Sa
    return S
end


function sumexp(vec::Array{Float64,1})
    mysum = 0.0
    @inbounds @simd for i=1:length(vec)
        mysum += exp(vec[i])
    end
    return mysum
end


function plmdca_asym2(Z::Array{T,2},W::Vector{Float64};
                decimation::Bool=false,
                fracmax::Real=0.3,
                fracdec::Real=0.1,
                remove_dups::Bool=true,
                min_separation::Int=1,
                theta=:auto,
                lambdaJ::Real=0.01,
                lambdaH::Real=0.01,
                epsconv::Real=1.0e-6,
                maxit::Int=5000,
                verbose::Bool=true,
                method::Symbol=:LD_LBFGS) where T <: Integer

    all(x -> x > 0, W) || throw(DomainError("vector W should normalized and with all positive elements"))
    isapprox(sum(W), 1) || throw(DomainError("sum(W) ≠ 1. Consider normalizing the vector W"))
    N, M = size(Z)
    M = length(W)
    q = Int(maximum(Z))

    plmalg = PlmAlg(method, verbose, epsconv, maxit)
    plmvar = PlmVar(N, M, q, q * q, lambdaJ, lambdaH, Z, W)
    Jmat, pslike = if decimation  == false
        MinimizePLAsym(plmalg, plmvar)
    else
        decvar = DecVar{2}(fracdec, fracmax, ones(Bool, (N - 1) * q * q, N))
        DecimateAsym!(plmvar, plmalg, decvar)
    end
    #score, FN, Jtensor, htensor =  ComputeScore(Jmat, plmvar, min_separation)
    return Jmat#PlmOut(sdata(pslike), Jtensor, htensor, score)

end
plmdca2(Z,W;kwds...) = plmdca_asym(Z, W;kwds...)

function plmdca_asym2(filename::String;
                theta::Union{Symbol,Real}=:auto,
                max_gap_fraction::Real=0.9,
                remove_dups::Bool=true,
                kwds...)
    time = @elapsed W, Z, N, M, q = ReadFasta(filename, max_gap_fraction, theta, remove_dups)
    println("preprocessing took $time seconds")
    plmdca_asym2(Z, W; kwds...)
end

plmdca2(filename::String; kwds...) = plmdca_asym(filename; kwds...)

function MinimizePLAsym(alg::PlmAlg, var::PlmVar)

    LL = (var.N - 1) * var.q2 + var.q
    x0 = zeros(Float64, LL)
    vecps = SharedArray{Float64}(var.N)
    Jmat = zeros(LL,var.N) |> SharedArray
    Threads.@threads for site = 1:var.N # 1:12
        opt = Opt(alg.method, length(x0))
        ftol_abs!(opt, alg.epsconv)
        xtol_rel!(opt, alg.epsconv)
        xtol_abs!(opt, alg.epsconv)
        ftol_rel!(opt, alg.epsconv)
        maxeval!(opt, alg.maxit)
        min_objective!(opt, (x, g) -> optimfunwrapper(x, g, site, var))
        elapstime = @elapsed  (minf, minx, ret) = optimize(opt, x0)
        alg.verbose && @printf("site = %d\t pl = %.4f\t time = %.4f\t", site, minf, elapstime)
        alg.verbose && println("exit status = $ret")
        vecps[site] = minf
        Jmat[:,site] .= minx
    end
    # Jmat = @distributed hcat for site = 1:var.N # 1:12
    #     opt = Opt(alg.method, length(x0))
    #     ftol_abs!(opt, alg.epsconv)
    #     xtol_rel!(opt, alg.epsconv)
    #     xtol_abs!(opt, alg.epsconv)
    #     ftol_rel!(opt, alg.epsconv)
    #     maxeval!(opt, alg.maxit)
    #     min_objective!(opt, (x, g) -> optimfunwrapper(x, g, site, var))
    #     elapstime = @elapsed  (minf, minx, ret) = optimize(opt, x0)
    #     alg.verbose && @printf("site = %d\t pl = %.4f\t time = %.4f\t", site, minf, elapstime)
    #     alg.verbose && println("exit status = $ret")
    #     vecps[site] = minf
    #     minx
    # end
    return sdata(Jmat), vecps
end

function ComputeUL(alg::PlmAlg, var::PlmVar, site::Int, LL::Int)

    N  = var.N
    q2 = var.q2
    lb = -Inf * ones(Float64, LL)
    ub =  Inf * ones(Float64, LL)
    tiny::Float64 = 1.0e-6
    offset::Int = 0

    for i = 1:site - 1
        for s = 1:q2
            lb[offset + s] = -tiny
            ub[offset + s] =  tiny
        end
    end
    offset += q2

    for i = site + 1:N
        for s = 1:q2
            lb[offset + s] = -tiny
            ub[offset + s] =  tiny
        end
        offset += q2
    end
    return lb, ub
end

function PLsiteAndGrad!(x::Vector{Float64}, grad::Vector{Float64}, site::Int, plmvar::PlmVar)


    LL = length(x)
    q2 = plmvar.q2
    q = plmvar.q
    N = plmvar.N
    M = plmvar.M
    Z = sdata(plmvar.Z)
    W = sdata(plmvar.W)
    IdxZ = sdata(plmvar.IdxZ)

    for i = 1:LL - q
	    grad[i] = 2.0 * plmvar.lambdaJ  * x[i]
	end
	for i = (LL - q + 1):LL
	    grad[i] = 4.0 * plmvar.lambdaH * x[i]
	end

	pseudolike = 0.0
	vecene = zeros(Float64, q)
	lnorm = 0.0
	expvecenesumnorm = zeros(Float64, q)

    @inbounds for m = 1:M
        izm = view(IdxZ, :, m)
        zsm = Z[site,m]
		fillvecene!(vecene, x, site, izm, q, N)
	    lnorm = logsumexp(vecene)
	    expvecenesumnorm .= @. exp(vecene - lnorm)
	    pseudolike -= W[m] * (vecene[ zsm ] - lnorm)
        @avx for i = 1:site - 1
	        for s = 1:q
	            grad[ izm[i] + s ] += W[m] * expvecenesumnorm[s]
	        end
	        grad[ izm[i] + zsm ] -= W[m]
	    end
        @avx for i = site + 1:N
	        for s = 1:q
	            grad[ izm[i] - q2 + s ] += W[m] *  expvecenesumnorm[s]
	        end
	        grad[ izm[i] - q2 + zsm ] -= W[m]
	    end
        @avx for s = 1:q
	        grad[ (N - 1) * q2 + s ] += W[m] * expvecenesumnorm[s]
	    end
		grad[ (N - 1) * q2 + zsm ] -= W[m]
	end
	pseudolike += L2norm_asym(x, plmvar)
end


# Energy filling

function fillvecene!(vecene::Vector{Float64}, x::Vector{Float64}, site::Int, IdxSeq::AbstractArray{Int,1}, q::Int, N::Int)
#     println(typeof(vecene))
#     println(typeof(x))
#     println(typeof(site))
#     println(typeof(IdxSeq))
#     println(typeof(q))
	q2 = q * q
    @inbounds for l = 1:q
        scra::Float64 = 0.0
        @avx for i = 1:site - 1 # Begin sum_i \neq site J
            scra += x[IdxSeq[i] + l]
        end
        # skipping sum over residue site
        @avx for i = site + 1:N
            scra +=  x[IdxSeq[i] - q2 + l]
        end # End sum_i \neq site J
        scra +=  x[(N - 1) * q2 + l] # sum H
        vecene[l] = scra
    end

end


function logsumexp(X::Vector{Float64})
    u = maximum(X)
    isfinite(u) || return float(u)
    return u + log(sum(x -> exp(x - u), X))
end


function L2norm_asym(vec::Array{Float64,1}, plmvar::PlmVar)
    q = plmvar.q
    N = plmvar.N
    lambdaJ = plmvar.lambdaJ
    lambdaH = plmvar.lambdaH

    LL = length(vec)

    mysum1 = 0.0
    @inbounds @avx for i = 1:(LL - q)
        mysum1 += vec[i] * vec[i]
    end
    mysum1 *= lambdaJ

    mysum2 = 0.0
    @inbounds @avx for i = (LL - q + 1):LL
        mysum2 += vec[i] * vec[i]
    end
    mysum2 *= 2lambdaH

    return mysum1 + mysum2
end



function ComputeScore(Jmat::Array{Float64,2}, var::PlmVar, min_separation::Int)

    q = var.q
    N = var.N
    JJ = reshape(Jmat[1:end - q,:], q, q, N - 1, N)
    Jtemp1 = zeros(q, q, Int(N * (N - 1) / 2))
    Jtemp2 = zeros(q, q, Int(N * (N - 1) / 2))
    l = 1
    for i = 1:(N - 1)
        for j = (i + 1):N
            Jtemp1[:,:,l] = JJ[:,:,j - 1,i] # J_ij as estimated from from g_i.
            Jtemp2[:,:,l] = JJ[:,:,i,j]' # J_ij as estimated from from g_j.
            l = l + 1
        end
    end

    hplm = fill(0.0, q, N)
    for i in 1:N
        hplm[:,i] .= Jmat[end - q + 1:end,i]
    end

    Jtensor1 = inflate_matrix(Jtemp1, N)
    Jtensor2 = inflate_matrix(Jtemp2, N)
    Jplm = (Jtensor1 + Jtensor2) / 2 # for the energy I do not want to gauge

    ctr = 0
    for i in 1:N - 1
        for j in i + 1:N
            ctr += 1
            Jtensor1[:,:,i,j] = Jtemp1[:,:,ctr] - repeat(mean(Jtemp1[:,:,ctr], dims=1), q, 1) - repeat(mean(Jtemp1[:,:,ctr], dims=2), 1, q) .+ mean(Jtemp1[:,:,ctr])
            Jtensor1[:,:,j,i] = Jtensor1[:,:,i,j]'
            Jtensor2[:,:,i,j] = Jtemp2[:,:,ctr] - repeat(mean(Jtemp2[:,:,ctr], dims=1), q, 1) - repeat(mean(Jtemp2[:,:,ctr], dims=2), 1, q) .+ mean(Jtemp2[:,:,ctr])
            Jtensor2[:,:,j,i] = Jtensor2[:,:,i,j]'
        end
    end # zerosumgauge the different tensors

    Jtensor = (Jtensor1 + Jtensor2) / 2

    FN = compute_APC(Jtensor, N, q)
    score = compute_ranking(FN, min_separation)
    return score, FN, Jplm, hplm
end




function fillvecene2!(vecene::Vector{Float64}, x::Vector{Float64}, site::Int, IdxSeq::AbstractArray{Int,1}, q::Int, N::Int)
#     printf(typeof(vecene))
#     printf(typeof(x))
#     printf(typeof(site))
#     printf(typeof(IdxSeq))
#     printf(typeof(q))
	q2 = q * q
    @inbounds for l = 1:q
        scra::Float64 = 0.0
        for i = 1:site - 1 # Begin sum_i \neq site J
            scra += x[IdxSeq[i] + l]
        end
        # skipping sum over residue site
        for i = site + 1:N
            scra +=  x[IdxSeq[i] - q2 + l]
        end # End sum_i \neq site J
        scra +=  x[(N - 1) * q2 + l] # sum H
        vecene[l] = scra
    end
end

function get_proba(Jmat, site, izm, plmvar)
    x = Jmat[:,site]
    LL = length(x)

    q = plmvar.q
    q2 = q*q
    N = plmvar.N
    M = plmvar.M
#     Z = sdata(plmvar.Z)
#     W = sdata(plmvar.W)
#     IdxZ = sdata(plmvar.IdxZ)
    vecene = zeros(Float64, q)
    lnorm = 0.0
    expvecenesumnorm = zeros(Float64, q)

#     @inbounds for m = 1:M
#         izm = view(IdxZ, :, m)
#         zsm = Z[site,m]
    fillvecene2!(vecene, x, site, izm, q, N)
    lnorm = logsumexp(vecene)
    expvecenesumnorm .= @. exp(vecene - lnorm)
    return expvecenesumnorm
end


function gibbsstep(Jmat, plmVar)
    q = plmVar.q
    @threads for m =1:plmVar.M #Threads
        izm = view(plmVar.IdxZ, :, m)
        site =rand(1:plmVar.N)
        p = get_proba(Jmat, site, izm,plmVar)
        d = Distributions.Categorical(p)
        l = rand(d)

#         println("$m $site $l")
        plmVar.Z[site, m] = l
        plmVar.IdxZ[site, m] = (site-1) * q*q + q * (l - 1)
    end
end
