using CSV
using NPZ



alipath = "../plm/data/CALM/aliCALM1_clean.fasta"
wtpath = "../plm/data/CALM/CALM1_clean.fasta"
mutationpath = "../plm/data/CALM/Roth2017.csv"
maskpath = "../plm/data/calm1mask.npz.npy"
csv = DataFrame(CSV.File(mutationpath))
processCSV(csv, maskpath)

theta=:auto
max_gap_fraction=0.9
remove_dups=true
time = @elapsed W, Z, N, M, q = ReadFasta(alipath, max_gap_fraction, theta, remove_dups)

lambdaJ =0.0001
lambdaH=0.000001
plmvar = PlmVar(N, M, q, q * q, lambdaJ, lambdaH, Z, W)

plmo = plmdca_asym2(joinpath(pwd(), alipath), theta = :auto, lambdaH=0.000001, lambdaJ=0.0001)

plmscore, dmsplmscores, dmsexpscores = DMS_score_plmsite(plmo, plmVar_wt, CSV)

_, dmsplmscores, dmsexpscores = DMS_score_plmsite(Jmat, plmVar_wt, CSV)
spplm = corspearman(dmsplmscores, dmsexpscores)


_, prof = profile(Pi_true, N, 0.05)
_, dmsplmscores, dmsexpscores = DMS_score_plmXprofile(plmo, prof, plmVar_wt, CSV)
spplmprof = corspearman(dmsplmscores, dmsexpscores)


arnet,arvar=ardca(alipath, verbose=false, lambdaJ=0.002,lambdaH=0.0001; permorder=:NATURAL)
ardms = dms_single_site(arnet, arvar, 1)[1]
dmsplmscores, dmsexpscores = DMS_score_ardca(ardms, csv)
spardca = corspearman(dmsplmscores, dmsexpscores)
