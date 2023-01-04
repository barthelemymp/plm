
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
using Base.Threads: @threads, nthreads


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
xticks = (1:9, ["BRCA1","BRCA1BRCT","CALM1","DLG4","AMIE","B3VI55T","BF520","HG"])
nchains = 50000
spplm_list = []
spplm_full_list = []
spplmprof_list = []
spardca_list = []
xs = []

startposlist = dict("AMIE"=>1,"B3VI55T"=>1,"BF520"=>25,"BRCA1"=>1,"BRCA1BRCT"=>1625,"CALM1"=>1,"DLG4"=>300,"HG" =>1) 
multss = []
for famname in ["BRCA1BRCT","CALM1","VIM","DLG4","AMIE","B3VI55T","BF520","HG", "BRCA1"]#"BLAT"
    global indi+=1
    # startpos = startposlist[famname]#startposlist[indi]
    @show indi, famname
    push!(xs,indi)
    datadir = "/Data/barth/mutdata/"
    alipath = datadir* famname*"/"*"ali"*famname*"_clean.fasta"
    # "../plm/data/CALM/aliCALM1_clean.fasta"
    # wtpath =  datadir * famname*"/"*famname*"_clean.fasta"
    # # "../plm/data/CALM/CALM1_clean.fasta"
    # mutationpath =  datadir * famname*"/"*"exp"*famname*"_clean.csv"
    # # "../plm/data/CALM/Roth2017.csv"
    # maskpath = datadir * famname*"/"* famname*"mask.npz.npy"
    # # "../plm/data/calm1mask.npz.npy"
    # csv = DataFrame(CSV.File(mutationpath, delim=";"))
    # processCSV(csv, maskpath, startpos)

    theta=:auto
    max_gap_fraction=0.9
    remove_dups=true
    time = @elapsed W, Z, N, M, q = ReadFasta(alipath, max_gap_fraction, theta, remove_dups)
#     Wwt, Zwt, N, Mwt, q = ReadFasta(alipath, max_gap_fraction, theta, remove_dups)
# end
    lambdaJ = 0.0001
    lambdaH = 0.000001
    plmvar = PlmVar(N, M, q, q * q, lambdaJ, lambdaH, Z, W)
    Pi_true, Pij_true, _, _ = compute_weighted_frequencies(convert(Array{Int8,2}, plmvar.Z), q, :auto)


# mypottsplm = PottsModel(AlignmentTest.sequence_length)
# plmo = plmdca_asym(joinpath(pwd(), "tmp/$(familyName)_$i.fasta"), theta = :auto)
# plm_couplings = deflate_matrix(plmo.Jtensor)
# plm_fields = plmo.htensor
# mypottsplm.fields .= plm_fields
# mypottsplm.couplings .= plm_couplings
# energies = energy(mypottsplm, AlignmentTest)
# spearmanerror = eval_spearman(mypottsplm, AlignmentTest)


    
    plmo = plmdca_asym2(joinpath(pwd(), alipath), theta = :auto,verbose=false, lambdaJ=lambdaJ,lambdaH=lambdaH)
    plmo2 = symetrize_matrix(plmo, q)
    
    #Pi_true, Pij_true, _, _ = compute_weighted_frequencies(convert(Array{Int8,2}, plmvar.Z), plmvarSample.q, :auto)
    pitrue, pijtrue = expandP(Pi_true, Pij_true,N)



    # mult = 2 #multss[i]
    # plmvarSample = PlmVar(N, M*mult, q, q * q, lambdaJ, lambdaH, transpose(repeat(transpose(Z),mult)), repeat(W,mult))

    plmvarSample = PlmVar(N, nchains, q, q * q, lambdaJ, lambdaH, ones(size(Z)[1], nchains), ones(nchains))
    corr_list = []
    x_list = []
    for s =0:2000
        gibbsstep(plmo, plmvarSample)
        if s%100==0
            @show "ëvalin"
            Pi_s, Pij_s, _, _ = compute_weighted_frequencies(convert(Array{Int8,2}, plmvarSample.Z), plmvarSample.q, 0)
            corrij_s = corrCIJ(Pi_s, Pij_s, N)
            corrij_true = corrCIJ(Pi_true, Pij_true, N)
            # println(Statistics.cor(vec(corrij_s), vec(corrij_true)))
            push!(corr_list,Statistics.cor(vec(corrij_s), vec(corrij_true)))
            push!(x_list, s)
        end
    end


    Pi_s, Pij_s, _, _ = compute_weighted_frequencies(convert(Array{Int8,2}, plmvarSample.Z), plmvarSample.q, 0)
    corrij_s = corrCIJ(Pi_s, Pij_s, N)
    corrij_true = corrCIJ(Pi_true, Pij_true, N)

    plt2 = plot(title="$famname corr for lh $(lambdaH) lj $(lambdaJ)",margins = 5Plots.mm)
    scatter!(plt2, vec(corrij_true), vec(corrij_s))
    xlabel!(plt2, "2pt corr True")
    ylabel!(plt2, "2pt corr Sample")
    savefig(plt2, "../../$(famname)corrScatter_lh$(lambdaH)_lj$(lambdaJ)_asym.png")


    # plt = plot(title="$famname corr for lh $(lambdaH) lj $(lambdaJ)",margins = 5Plots.mm)
    # ylims!((0.0,1.0))
    # xlims!((0,2000))
    # plot!(plt,x_list, corr_list, label="Gap Init")
    # plot!(plt,x_list2, corr_list2, label="Data Init")
    # xlabel!(plt, "gibsteps")
    # ylabel!(plt, "2pt correlation")
    # savefig(plt, "../../$(famname)corr_lh$(lambdaH)_lj$(lambdaJ).png")

    # plt2 = plot(title="$famname corr for lh $(lambdaH) lj $(lambdaJ)",margins = 5Plots.mm)
    # scatter!(plt2, vex(corrij_true), vec(corrij_s))
    # xlabel!(plt2, "2pt corr True")
    # ylabel!(plt2, "2pt corr Sample")
    # savefig(plt, "../../$(famname)corrScatter_lh$(lambdaH)_lj$(lambdaJ).png")

    # plt3 = plot(title="$famname pij for lh $(lambdaH) lj $(lambdaJ)",margins = 5Plots.mm)
    # scatter!(plt2, vex(pitrue), vec(pis))
    # xlabel!(plt2, "2pt corr True")
    # ylabel!(plt2, "2pt corr Sample")
    # savefig(plt, "../../$(famname)f1Scatter_lh$(lambdaH)_lj$(lambdaJ).png")

    # plt4 = plot(title="$famname pi for lh $(lambdaH) lj $(lambdaJ)",margins = 5Plots.mm)
    # scatter!(plt2, vex(pijtrue), vec(pijs))
    # xlabel!(plt2, "2pt corr True")
    # ylabel!(plt2, "2pt corr Sample")
    # savefig(plt, "../../$(famname)f2Scatter_lh$(lambdaH)_lj$(lambdaJ).png")







    plmvarSample = PlmVar(N, nchains, q, q * q, lambdaJ, lambdaH, ones(size(Z)[1], nchains), ones(nchains))
    corr_list2 = []
    x_list2 = []
    
    for s =0:2000
        gibbsstep(plmo2, plmvarSample)
        if s%100==0
            Pi_s, Pij_s, _, _ = compute_weighted_frequencies(convert(Array{Int8,2}, plmvarSample.Z), plmvarSample.q, 0)
            corrij_s = corrCIJ(Pi_s, Pij_s, N)
            corrij_true = corrCIJ(Pi_true, Pij_true, N)
            # println(Statistics.cor(vec(corrij_s), vec(corrij_true)))
            push!(corr_list2,Statistics.cor(vec(corrij_s), vec(corrij_true)))
            push!(x_list2, s)
        end
    end
    pis, pijs = expandP(Pi_s, Pij_s,N)
    plt = plot(title="$famname corr for lh $(lambdaH) lj $(lambdaJ)",margins = 5Plots.mm)
    ylims!((0.0,1.0))
    xlims!((0,2000))
    plot!(plt,x_list, corr_list, label="asym")
    plot!(plt,x_list2, corr_list2, label="sym")
    xlabel!(plt, "gibsteps")
    ylabel!(plt, "2pt correlation")
    savefig(plt, "../../$(famname)corr_lh$(lambdaH)_lj$(lambdaJ)_asymvssym.png")

    plt2 = plot(title="$famname corr for lh $(lambdaH) lj $(lambdaJ)",margins = 5Plots.mm)
    scatter!(plt2, vec(corrij_true), vec(corrij_s))
    xlabel!(plt2, "2pt corr True")
    ylabel!(plt2, "2pt corr Sample")
    savefig(plt2, "../../$(famname)corrScatter_lh$(lambdaH)_lj$(lambdaJ)_sym.png")

    # plt3 = plot(title="$famname pij for lh $(lambdaH) lj $(lambdaJ)",margins = 5Plots.mm)
    # scatter!(plt2, vex(pitrue), vec(pis))
    # xlabel!(plt2, "2pt corr True")
    # ylabel!(plt2, "2pt corr Sample")
    # savefig(plt, "../../$(famname)f1Scatter_lh$(lambdaH)_lj$(lambdaJ)_sym.png")

    # plt4 = plot(title="$famname pi for lh $(lambdaH) lj $(lambdaJ)",margins = 5Plots.mm)
    # scatter!(plt2, vex(pijtrue), vec(pijs))
    # xlabel!(plt2, "2pt corr True")
    # ylabel!(plt2, "2pt corr Sample")
    # savefig(plt, "../../$(famname)f2Scatter_lh$(lambdaH)_lj$(lambdaJ)_sym.png")




    # thetawt = 0.0
    # Wwt, Zwt, N, Mwt, q = ReadFasta(wtpath, max_gap_fraction, thetawt, remove_dups)
    #
    # q=21
    #
    # plmVar_wt = PlmVar(N, Mwt, q, q * q, lambdaJ, lambdaH, Zwt, Wwt)
    #
    # plmscore, dmsplmscores, dmsexpscores = DMS_score_plmsite(plmo, plmVar_wt, csv)
    # spplm = corspearman(dmsplmscores, dmsexpscores)
    # push!(spplm_list,spplm)
    #
    # dmsplmscores, dmsexpscores = DMS_score_plm(plmo, plmVar_wt, csv)
    # spplm_full = corspearman(dmsplmscores, dmsexpscores)
    # push!(spplm_full_list,spplm_full)
    #
    #
    # _, prof = profile(Pi_true, N, 0.05)
    # _, dmsplmscores, dmsexpscores = DMS_score_plmXprofile(plmo, prof, plmVar_wt, csv)
    # spplmprof = corspearman(dmsplmscores, dmsexpscores)
    # push!(spplmprof_list,spplmprof)
    # arnet,arvar=ardca(alipath, verbose=false, lambdaJ=0.01,lambdaH=0.0001; permorder=:NATURAL)
    # ardms = dms_single_site(arnet, arvar, 1)[1]
    # dmsplmscores, dmsexpscores = DMS_score_ardca(ardms, csv)
    # spardca =-1* corspearman(dmsplmscores, dmsexpscores)
    # push!(spardca_list,spardca)
    # @show spplm, spplmprof, spardca
end

scatter!(plt, xs,spplm_list, label="cond_plm")
scatter!(plt, xs,spplm_full_list, label="plm")
scatter!(plt, xs,spplmprof_list, label="plm profile rebalanced")
scatter!(plt, xs,spardca_list, label="ardca")
xticks!(plt, xticks)

savefig(plt, "../../mut.png")

