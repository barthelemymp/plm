
using JLD: save, load
#cd("C:\\Users\\bartm\\Documents\\These\\plm\\test")
cd("/Data/barth/plm/test")
push!(LOAD_PATH, joinpath(pwd(), "../src"))
using Revise
using plm
using Statistics
using CSV
using ArDCA
using ArgParse
using Plots
using DataFrames
using CSV
using NPZ
using StatsBase
using Plots
global indi = 0
plt = plot(title="Mutation spearman correlation")
xticks = (1:9, ["AMIE","B3VI55T","BF520","BRCA1","BRCA1BRCT","CALM1","DLG4","HG"])

spplm_list = []
spplm_full_list = []
spplmprof_list = []
spardca_list = []
xs = []

theta=:auto
max_gap_fraction=0.9
remove_dups=true


startposlist = [1     ,1         ,25    ,1      ,1625       ,1      ,300   ,1]
# multss = []
for famname in ["BRCA1","B3VI55T","BF520","AMIE","BRCA1BRCT","CALM1","DLG4","HG"]#"BLAT"
    global indi+=1
    startpos = startposlist[indi]
    @show indi, famname
    push!(xs,indi)
    datadir = "/Data/barth/mutdata/"
    alipath = datadir* famname*"/"*"ali"*famname*"_clean.fasta"
    wtpath =  datadir * famname*"/"*famname*"_clean.fasta"
    mutationpath =  datadir * famname*"/"*famname*"_mutants_exp.fasta"
    fitnesspath =  datadir * famname*"/"*famname*"_mutants_exp.fit"
    sequences, mutations = read_mutation(mutationpath)
    mutations = mutations[2:]
    fitness = read_fit(fitnesspath)
    fitness = fitness[2:]
    time = @elapsed W, Z, N, M, q = ReadFasta(alipath, max_gap_fraction, theta, remove_dups)
# if generative
    # lambdaJ = 0.0001
    # lambdaH = 0.000001

# if mds
    lambdaJ = 0.01
    lambdaH = 0.0001
    plmvar = PlmVar(N, M, q, q * q, lambdaJ, lambdaH, Z, W)
    Pi_true, Pij_true, _, _ = compute_weighted_frequencies(convert(Array{Int8,2}, plmvar.Z), q, :auto)
    pitrue, pijtrue = expandP(Pi_true, Pij_true,N)
    plmo = plmdca_asym2(joinpath(pwd(), alipath), theta = :auto,verbose=false, lambdaJ=lambdaJ,lambdaH=lambdaH)

    thetawt = 0.0
    Wwt, Zwt, N, Mwt, q = ReadFasta(wtpath, max_gap_fraction, thetawt, remove_dups)
    q=21
    plmVar_wt = PlmVar(N, Mwt, q, q * q, lambdaJ, lambdaH, Zwt, Wwt)
    plmscore, dmsplmscores, dmsexpscores = DMS_score_plmsite(plmo, plmVar_wt, mutations, fitness)
    spplm = corspearman(dmsplmscores, dmsexpscores)
    push!(spplm_list,spplm)
    dmsplmscores, dmsexpscores = DMS_score_plm(plmo, plmVar_wt, mutations, fitness)
    spplm_full = corspearman(dmsplmscores, dmsexpscores)
    push!(spplm_full_list,spplm_full)
    _, prof = profile(Pi_true, N, 0.05)
    _, dmsplmscores, dmsexpscores = DMS_score_plmXprofile(plmo, prof, plmVar_wt, mutations, fitness)
    spplmprof = corspearman(dmsplmscores, dmsexpscores)
    push!(spplmprof_list,spplmprof)
    arnet,arvar=ardca(alipath, verbose=false, lambdaJ=0.01,lambdaH=0.0001; permorder=:NATURAL)
    ardms = dms_single_site(arnet, arvar, 1)[1]
    dmsplmscores, dmsexpscores = DMS_score_ardca(ardms, mutations, fitness)
    spardca = -1* corspearman(dmsplmscores, dmsexpscores)
    push!(spardca_list,spardca)
    @show famname, spplm, spplm_full, spplmprof, spardca
end

# scatter!(plt, xs,spplm_list, label="cond_plm")
# scatter!(plt, xs,spplm_full_list, label="plm")
# scatter!(plt, xs,spplmprof_list, label="plm profile rebalanced")
# scatter!(plt, xs,spardca_list, label="ardca")
# xticks!(plt, xticks)
#
# savefig(plt, "../../mut.png")
