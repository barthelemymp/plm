from Bio import SeqIO
import numpy as np

def makeMatchmap(seq):
    N = len(seq[0])
    maskMatch = [True]*N
    for i in range(len(seq)):
        assert len(seq[i]) == N
        for po in range(N):
            if seq[i][po] == ".":
                maskMatch[po] = False
    return maskMatch

def write_fasta(dictionary, filename):
    """
    Takes a dictionary and writes it to a fasta file
    Must specify the filename when caling the function
    """
    with open(filename, "w") as outfile:
        for key, value in dictionary.items():
            outfile.write(">"+key + "\n")
            outfile.write(value)
            outfile.write("\n")

def processfasta(fastaPath, file_out, wtPath,wtPath_out, masksavepath):
    seqs = []
    names = []
    for records in SeqIO.parse(fastaPath,"fasta"):
        names.append(records.id )
        seqs.append([x for x in list(records.seq)])

    maskMatch = makeMatchmap(seqs)
    np.save(masksavepath,maskMatch)

    md = {}
    for records in SeqIO.parse(fastaPath,"fasta"):
        names.append(records.id )
        se = "".join([list(records.seq)[i] for i in range(len(list(records.seq))) if maskMatch[i]])
        md[records.id ] = se

    write_fasta(md, file_out)

    mdwt = {}
    for records in SeqIO.parse(wtPath,"fasta"):
        names.append(records.id )
        se = "".join([list(records.seq)[i] for i in range(len(list(records.seq))) if maskMatch[i]])
        mdwt[records.id ] = se

    write_fasta(mdwt, file_out)





for famname in ["AMIE","B3VI55T","BF520","BLAT","BRCA1","BRCA1BRCT","CALM1","DLG4","HG"]: #"AMIE",
    print(famname)
    path_dir = "/Data/barth/mutdata/"
    fastaPath = path_dir + famname+"/"+ "ali"+famname+".fasta"
    wtPath = path_dir + famname+"/"+ famname+".fasta"
    wtPath_out = path_dir + famname+"/"+ famname+"_clean.fasta"
    file_out = path_dir + famname+"/"+ "ali"+famname+"_clean.fasta"
    masksavepath = path_dir + famname+"/"+ famname+"mask.npz"
    processfasta(fastaPath, file_out, masksavepath)
    processfasta(wtPath, wtPath_out, masksavepath)
