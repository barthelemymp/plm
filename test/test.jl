
using JLD: save, load
#cd("C:\\Users\\bartm\\Documents\\These\\plm\\test")
cd("/Data/barth/plm/test")
push!(LOAD_PATH, joinpath(pwd(), "../src"))
using Revise
using plm
using Statistics
using CSV
using ArgParse
using Plots
function parse_commandline()
    s = ArgParseSettings()
    s = ArgParseSettings()
    @add_arg_table s begin
        "--lambdaH"
            help = "Lambda H"
            default = 0.000001
        "--lambdaJ"
            help = "Lambda J"
            default =  0.0001
        "--mult"
            help = "Multiplier"
            arg_type = Int
            default =  1
        "--rinit"
            help = "initialize ramdomly the chain or fron the datasequence"
            action = :store_true
        "fasta"
            help = "path to  fasta"
            required = true
    end
    return parse_args(s)
end

function main()
    parsed_args = parse_commandline()
    println("Parsed args:")
    for (arg,val) in parsed_args
        println("  $arg  =>  $val")
    end

    lambdaH = parse(Float64,parsed_args["lambdaH"])
    lambdaJ =  parse(Float64,parsed_args["lambdaJ"])
    @show lambdaH typeof(lambdaH)
    msadir =  parsed_args["fasta"]
    plmo = plmdca_asym2(joinpath(pwd(), msadir), theta = :auto, lambdaH=lambdaH, lambdaJ=lambdaJ)

    filename = msadir
    theta=:auto
    max_gap_fraction=0.9
    remove_dups=true
    time = @elapsed W, Z, N, M, q = ReadFasta(filename, max_gap_fraction, theta, remove_dups)

    N, M = size(Z)
    M = length(W)
    q = Int(maximum(Z))

    plmvar = PlmVar(N, M, q, q * q, lambdaJ, lambdaH, Z, W)

    mult= parsed_args["mult"]
    plmvarSample = PlmVar(N, plmvar.M*mult, q, q * q, lambdaJ, lambdaH, reshape(repeat(plmvar.Z, mult),(N,:)),repeat(plmvar.W, mult))
    if  parsed_args["rinit"]
        plmvarSample.Z.=21
    end

    # Pi_s, Pij_s, _, _ = compute_weighted_frequencies(convert(Array{Int8,2}, plmvarSample.Z), plmvarSample.q, 0)
    Pi_true, Pij_true, _, _ = compute_weighted_frequencies(convert(Array{Int8,2}, plmvar.Z), plmvarSample.q, :auto)
    # println(Statistics.cor(vec(Pij_s), vec(Pij_true)))
    corr_list = []
    x_list = []

    for s =1:1500
        gibbsstep(plmo, plmvarSample)
        if s%50==0
            Pi_s, Pij_s, _, _ = compute_weighted_frequencies(convert(Array{Int8,2}, plmvarSample.Z), plmvarSample.q, 0)
            # Pi_true, Pij_true, _, _ = compute_weighted_frequencies(convert(Array{Int8,2}, plmvar.Z), plmvarSample.q, :auto)
            corrij_s = corrCIJ(Pi_s, Pij_s, N)
            corrij_true = corrCIJ(Pi_true, Pij_true, N)
            println(Statistics.cor(vec(corrij_s), vec(corrij_true)))
            push!(corr_list,Statistics.cor(vec(corrij_s), vec(corrij_true)))
            push!(x_list, s)
        end
    end



    plmvarSample = PlmVar(N, M*mult, q, q * q, lambdaJ, lambdaH, transpose(repeat(transpose(Z),mult)), repeat(W,mult))
    corr_list2 = []
    x_list2 = []

    for s =1:1500
        gibbsstep(plmo, plmvarSample)
        if s%50==0
            Pi_s, Pij_s, _, _ = compute_weighted_frequencies(convert(Array{Int8,2}, plmvarSample.Z), plmvarSample.q, 0)
            # Pi_true, Pij_true, _, _ = compute_weighted_frequencies(convert(Array{Int8,2}, plmvar.Z), plmvarSample.q, :auto)
            corrij_s = corrCIJ(Pi_s, Pij_s, N)
            corrij_true = corrCIJ(Pi_true, Pij_true, N)
            println(Statistics.cor(vec(corrij_s), vec(corrij_true)))
            push!(corr_list2,Statistics.cor(vec(corrij_s), vec(corrij_true)))
            push!(x_list2, s)
        end
    end
    plt = plot(title="corr for lh $(lambdaH) lj $(lambdaJ)",margins = 5Plots.mm)
    plot!(plt,x_list, corr_list, label="Gap Init")
    plot!(plt,x_list2, corr_list2, label="Data Init")
    xlabel!(plt, "gibsteps")
    ylabel!(plt, "2pt correlation")
    savefig(plt, "../../corr14_lh$(lambdaH)_lj$(lambdaJ).png")













end
main()


#
# lambdaH = 0.000001
# lambdaJ = 0.0001
# msadir = "../data/PF00014.fasta"
# plmo = plmdca_asym2(joinpath(pwd(), msadir), theta = :auto, lambdaH=lambdaH, lambdaJ=lambdaJ)
#
#
# filename = msadir
# theta=:auto
# max_gap_fraction=0.9
# remove_dups=true
# time = @elapsed W, Z, N, M, q = ReadFasta(filename, max_gap_fraction, theta, remove_dups)
#
#
# N, M = size(Z)
# M = length(W)
# q = Int(maximum(Z))
# lambdaJ = 0.01
# lambdaH=0.01
# plmvar = PlmVar(N, M, q, q * q, lambdaJ, lambdaH, Z, W)
#
# mult=4
# plmvarSample = PlmVar(N, plmvar.M*mult, q, q * q, lambdaJ, lambdaH, reshape(repeat(plmvar.Z, mult),(N,:)),repeat(plmvar.W, mult))
# plmvarSample.Z.=21
# Pi_s, Pij_s, _, _ = compute_weighted_frequencies(convert(Array{Int8,2}, plmvarSample.Z), plmvarSample.q, 0)
# Pi_true, Pij_true, _, _ = compute_weighted_frequencies(convert(Array{Int8,2}, plmvar.Z), plmvarSample.q, :auto)
# println(Statistics.cor(vec(Pij_s), vec(Pij_true)))
# for s =1:500
#     gibbsstep(plmo, plmvarSample)
#
#     if s%100==0
#         Pi_s, Pij_s, _, _ = compute_weighted_frequencies(convert(Array{Int8,2}, plmvarSample.Z), plmvarSample.q, 0)
#         # Pi_true, Pij_true, _, _ = compute_weighted_frequencies(convert(Array{Int8,2}, plmvar.Z), plmvarSample.q, :auto)
#         corrij_s = corrCIJ(Pi_s, Pij_s)
#         corrij_true = corrCIJ(Pi_true, Pij_true)
#         println(Statistics.cor(vec(corrij_s), vec(corrij_true)))
#     end
# end
