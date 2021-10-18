library(devtools)
setwd("~/surfdrive/Shared/Boxtel_General/Scripts/Git_submission/Freek_MutationalPatterns/MutationalPatterns/")
load_all()

ref_genome <- "BSgenome.Hsapiens.UCSC.hg19"
library(ref_genome, character.only = TRUE)
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
library(tidyverse)

setwd("~/surfdrive/Shared/Projects/Freek/mutpatterns/environmental_agents/")

#Read data and transform into correct format
snvs_tb = read_tsv("denovo_subclone_subs_final.txt")
gr = makeGRangesFromDataFrame(snvs_tb, start.field = "Pos", end.field = "Pos", keep.extra.columns = TRUE)
seqlevels(gr) = paste0("chr", c(1:22, "X", "Y"))
seqlengths(gr) = seqlengths(BSgenome.Hsapiens.UCSC.hg19)[paste0("chr", c(1:22, "X", "Y"))]
grl = split(gr, gr$Sample)

# Take samples with enough mutations. Then take the first three as an example.
grl = grl[elementNROWS(grl) > 1000]
grl = grl[1:3]
#names(grl) = paste0("DBADE_", 1:3)
# DBADE concentration: 0.109 uM
#  dibenz[a,h]anthracene diol-epoxide

# Create figure
lesion_fig = plot_lesion_segregation(grl, subsample = 0.33)
ggsave("lesion.pdf", lesion_fig)
ggsave("lesion.png", lesion_fig)

# Calculate data
rl20_tb = calculate_lesion_segregation(grl, names(grl), test = "rl20", 
                                         ref_genome = "BSgenome.Hsapiens.UCSC.hg19", chromosomes = seqlevels(gr))
write_tsv(rl20_tb, "lesion_rl20.txt")

lesion_tb = calculate_lesion_segregation(grl, names(grl))
write_tsv(lesion_tb, "lesion.txt")
