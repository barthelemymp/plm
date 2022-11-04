module plm
#
# abstract type Optimizer end
# abstract type SequenceModel end
# abstract type Gauge end
# abstract type FocusingSchedule end
# abstract type Distance end
# abstract type AnnealingSchedule end
include("utils.jl")
include("plmdca.jl")



export write_fasta_data,
       expandP,
       corrCIJ,
       PlmVar,
       ReadFasta,
       plmdca_asym2,
       get_proba,
       gibbsstep,
       compute_weighted_frequencies
end
