using CSV
using NPZ
using StatsBase
using Plots
ind = 0
plt = plot(title="Mutation spearman correlation")
xticks = (1:9, ["AMIE","B3VI55T","BF520","BLAT","BRCA1","BRCA1BRCT","CALM1","DLG4","HG"])

spplm_list = []
spplmprof_list = []
spardca_list = []
xs = []
for famname in ["AMIE","B3VI55T","BF520","BLAT","BRCA1","BRCA1BRCT","CALM1","DLG4","HG"]
    ind+=1
    push!(xs,ind)
    datadir = "/Data/barth/mutdata/"
    alipath = datadir* famname*"/"*"ali"*famname*"_clean.fasta"
    # "../plm/data/CALM/aliCALM1_clean.fasta"
    wtpath =  datadir * famname*"/"*famname*"_clean.fasta"
    # "../plm/data/CALM/CALM1_clean.fasta"
    mutationpath =  datadir * famname*"/"*"exp"*famname*".csv"
    # "../plm/data/CALM/Roth2017.csv"
    maskpath = datadir * famname*"/"* famname*"mask.npz"
    # "../plm/data/calm1mask.npz.npy"
    csv = DataFrame(CSV.File(mutationpath))
    processCSV(csv, maskpath)

    theta=:auto
    max_gap_fraction=0.9
    remove_dups=true
    time = @elapsed W, Z, N, M, q = ReadFasta(alipath, max_gap_fraction, theta, remove_dups)

    lambdaJ =0.002
    lambdaH=0.0001
    plmvar = PlmVar(N, M, q, q * q, lambdaJ, lambdaH, Z, W)

    plmo = plmdca_asym2(joinpath(pwd(), alipath), theta = :auto, lambdaJ=0.002,lambdaH=0.0001)

    plmscore, dmsplmscores, dmsexpscores = DMS_score_plmsite(plmo, plmVar_wt, CSV)
    _, dmsplmscores, dmsexpscores = DMS_score_plmsite(Jmat, plmVar_wt, CSV)
    spplm = corspearman(dmsplmscores, dmsexpscores)
    push!(spplm_list,spplm)
    _, prof = profile(Pi_true, N, 0.05)
    _, dmsplmscores, dmsexpscores = DMS_score_plmXprofile(plmo, prof, plmVar_wt, CSV)
    spplmprof = corspearman(dmsplmscores, dmsexpscores)
    push!(spplmprof_list,spplmprof)
    arnet,arvar=ardca(alipath, verbose=false, lambdaJ=0.002,lambdaH=0.0001; permorder=:NATURAL)
    ardms = dms_single_site(arnet, arvar, 1)[1]
    dmsplmscores, dmsexpscores = DMS_score_ardca(ardms, csv)
    spardca = corspearman(dmsplmscores, dmsexpscores)
    push!(spardca_list,spardca)
end
scatter!(plt, xs,spplm_list, label="plm")
scatter!(plt, xs,spplmprof_list, label="plm profile rebalanced")
scatter!(plt, xs,spardca_list, label="ardca")
xticks!(plt, xticks)

savefig(plt, "../../mut.png")
