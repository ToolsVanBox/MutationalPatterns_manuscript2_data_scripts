library(devtools)
library(tidyverse)
ref_genome <- "BSgenome.Hsapiens.UCSC.hg19"
library(ref_genome, character.only = TRUE)
setwd("~/surfdrive/Shared/Boxtel_General/Scripts/Git_submission/Freek_MutationalPatterns/MutationalPatterns/")
load_all()


out_dir = "~/surfdrive/Shared/Projects/Freek/mutpatterns/hmf/"
if (!dir.exists(out_dir)){
    dir.create(out_dir)
}
setwd(out_dir)

signatures = get_known_signatures()
mut_mat = read.table("~/surfdrive/Shared/Projects/Freek/mutpatterns/HMF_mut_mat_noPASS.txt", sep = "\t")
fit_res = fit_to_signatures(mut_mat, signatures)
highest_sbs22_samples = names(sort(prop.table(fit_res$contribution, 2)["SBS22",], decreasing = TRUE)[1:10])

hmf_overview_all = suppressWarnings(read_tsv("hmf_metadata.tsv"))
#hmf_overview = hmf_overview_all[hmf_overview_all$sampleId %in% highest_sbs22_samples,]
hmf_overview = hmf_overview_all[hmf_overview_all$primaryTumorLocation %in% c("Biliary", "Kidney", "Liver"),]
hmf_dir = "~/hpc/pmc_vanboxtel/projects/Axel_GenoEcoli/HMF_analysis/somatics"
vcf_fnames = file.path(hmf_dir, hmf_overview$setName, paste0(hmf_overview$sampleId, ".purple.somatic.vcf"))

#Set sample names
sample_names = hmf_overview$sampleId

#Read data
grl = read_vcfs_as_granges(vcf_fnames, sample_names, ref_genome, type = "all")
grl = purrr::map(as.list(grl), function(gr){gr = gr[gr$FILTER == "PASS"]}) %>% 
    purrr::map(.remove_multi_alts_variants) %>% 
    GRangesList()

snv_grl = purrr::map(as.list(grl), function(gr){gr = gr[.find_substitution(gr) & width(gr$REF) == 1]}) %>% 
    GRangesList()

t2a_grl = purrr::map(as.list(snv_grl), function(gr){gr = gr[mut_type(gr) == "T>A"]}) %>% 
    GRangesList()

rl20_t2a_tb = calculate_lesion_segregation(t2a_grl, names(t2a_grl), test = "rl20", ref_genome = ref_genome, chromosomes = seqlevels(t2a_grl[[1]]))
rl20_tb = calculate_lesion_segregation(snv_grl, names(snv_grl), test = "rl20", ref_genome = ref_genome, chromosomes = seqlevels(snv_grl[[1]]))
plot_lesion_segregation(t2a_grl[1:3])
