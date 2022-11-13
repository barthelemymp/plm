
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
startposlist = [1     ,1         ,25    ,1      ,1625       ,1      ,300   ,1]
for famname in ["AMIE","B3VI55T","BF520","BRCA1","BRCA1BRCT","CALM1","DLG4","HG"]#"BLAT"
    global indi+=1
    startpos = startposlist[indi]
    @show indi, famname
    push!(xs,indi)
    datadir = "/Data/barth/mutdata/"
    alipath = datadir* famname*"/"*"ali"*famname*"_clean.fasta"
    # "../plm/data/CALM/aliCALM1_clean.fasta"
    wtpath =  datadir * famname*"/"*famname*"_clean.fasta"
    # "../plm/data/CALM/CALM1_clean.fasta"
    mutationpath =  datadir * famname*"/"*"exp"*famname*"_clean.csv"
    # "../plm/data/CALM/Roth2017.csv"
    maskpath = datadir * famname*"/"* famname*"mask.npz.npy"
    # "../plm/data/calm1mask.npz.npy"
    csv = DataFrame(CSV.File(mutationpath, delim=";"))
    processCSV(csv, maskpath, startpos)

    theta=:auto
    max_gap_fraction=0.9
    remove_dups=true
    time = @elapsed W, Z, N, M, q = ReadFasta(alipath, max_gap_fraction, theta, remove_dups)

    lambdaJ =0.01
    lambdaH=0.0001
    plmvar = PlmVar(N, M, q, q * q, lambdaJ, lambdaH, Z, W)
    Pi_true, Pij_true, _, _ = compute_weighted_frequencies(convert(Array{Int8,2}, plmvar.Z), q, :auto)

    plmo = plmdca_asym2(joinpath(pwd(), alipath), theta = :auto, lambdaJ=0.002,lambdaH=0.0001)


    thetawt = 0.0
    Wwt, Zwt, N, Mwt, q = ReadFasta(wtpath, max_gap_fraction, thetawt, remove_dups)

    q=21

    plmVar_wt = PlmVar(N, Mwt, q, q * q, lambdaJ, lambdaH, Zwt, Wwt)

    plmscore, dmsplmscores, dmsexpscores = DMS_score_plmsite(plmo, plmVar_wt, csv)
    spplm = corspearman(dmsplmscores, dmsexpscores)
    push!(spplm_list,spplm)

    dmsplmscores, dmsexpscores = DMS_score_plm(plmo, plmVar_wt, csv)
    spplm_full = corspearman(dmsplmscores, dmsexpscores)
    push!(spplm_full_list,spplm_full)


    _, prof = profile(Pi_true, N, 0.05)
    _, dmsplmscores, dmsexpscores = DMS_score_plmXprofile(plmo, prof, plmVar_wt, csv)
    spplmprof = corspearman(dmsplmscores, dmsexpscores)
    push!(spplmprof_list,spplmprof)
    arnet,arvar=ardca(alipath, verbose=false, lambdaJ=0.01,lambdaH=0.0001; permorder=:NATURAL)
    ardms = dms_single_site(arnet, arvar, 1)[1]
    dmsplmscores, dmsexpscores = DMS_score_ardca(ardms, csv)
    spardca =-1* corspearman(dmsplmscores, dmsexpscores)
    push!(spardca_list,spardca)
    @show spplm, spplmprof, spardca
end

scatter!(plt, xs,spplm_list, label="cond_plm")
scatter!(plt, xs,spplm_full_list, label="plm")
scatter!(plt, xs,spplmprof_list, label="plm profile rebalanced")
scatter!(plt, xs,spardca_list, label="ardca")
xticks!(plt, xticks)

savefig(plt, "../../mut.png")
