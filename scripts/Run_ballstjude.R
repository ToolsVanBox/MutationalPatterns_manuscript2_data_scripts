library(devtools)
setwd("~/surfdrive/Shared/vanBoxtelLab (Groupfolder)/Projects/Freek/MutationalPatterns/")
load_all()
library(tidyverse)

#Helper function to read St.Jude like vcf file.
read_vcflike_sample = function(fname){
    vcf_sample = read_tsv(fname, col_names = F, comment = "#", col_types = "cicccdci", skip = 4)
    colnames(vcf_sample) = c("CHROM", "POS", "ID", "REF", "ALT", "VAF", "placeholder", "flanking")
    vcf_sample = vcf_sample %>%
        dplyr::select(-placeholder, -flanking)
}

#Function to read St.Jude like vcf files.
read_vcflike = function(fnames){
    vcfs = purrr::map(fnames, read_vcflike_sample) %>% 
        bind_rows()
    vcfs$sampleID = str_match(vcfs$ID, ".*-.*-(.*)-.*-.*")[,2]
    vcfs_gr = makeGRangesFromDataFrame(vcfs, keep.extra.columns = T, start.field = "POS", end.field = "POS")
    seqlevelsStyle(vcfs_gr) = "UCSC"
    chromosomes = levels(seqnames(vcfs_gr))
    seqlengths(vcfs_gr) = seqlengths(BSgenome.Hsapiens.UCSC.hg19)[chromosomes]
    return(vcfs_gr)
}
    
ref_genome = "BSgenome.Hsapiens.UCSC.hg19"

dir = "~/surfdrive/Shared/vanBoxtelLab (Groupfolder)/Projects/Freek/mutpatterns/local_mut_patterns/"
if (!dir.exists(dir)){
    dir.create(dir)
}
setwd(dir)

# Read in data
ball_vcf_fnames = list.files("~/surfdrive/Shared/vanBoxtelLab (Groupfolder)/Boxtel_General/Data/Mutation_data/SNVs/Published_data/ALL_AML_NBL_OS_WT_Ma_et_al_2018_Nature/BALL/", full.names = T)
ball_vcf_fnames = ball_vcf_fnames[grep(".*vcf$", ball_vcf_fnames)]
gr_total = read_vcflike(ball_vcf_fnames)
gr = unique(gr_total)
GenomeInfoDb::genome(gr) = 'hg19'

# Calculate cossim
chromosomes = paste0("chr", c(seq(1:22), "X"))
sims = determine_regional_similarity(gr, ref_genome, chromosomes, window_size = 200, stepsize = 50, oligo_correction = FALSE)
saveRDS(sims, "sims.rds")
sims_tri_corrected = determine_regional_similarity(gr, ref_genome, chromosomes, window_size = 200, stepsize = 50, oligo_correction = TRUE)
saveRDS(sims_tri_corrected, "sims_tri_corrected.rds")
sims_ext = determine_regional_similarity(gr, ref_genome, chromosomes, window_size = 200, stepsize = 50, oligo_correction = FALSE, extension = 2)
saveRDS(sims_ext, "sims_ext.rds")
sims_noext = determine_regional_similarity(gr, ref_genome, chromosomes, window_size = 200, stepsize = 50, oligo_correction = FALSE, extension = 0)
saveRDS(sims_noext, "sims_noext.rds")

# Plot results
sims_fig = plot_regional_similarity(sims)
ggsave("sims.pdf", sims_fig)
sims_fig_corrected = plot_regional_similarity(sims_tri_corrected, oligo_correction = TRUE)
ggsave("sims_corrected.pdf", sims_fig_corrected)
sims_fig_ext = plot_regional_similarity(sims_ext)
ggsave("sims_ext.pdf", sims_fig_ext)
sims_fig_noext = plot_regional_similarity(sims_noext)
ggsave("sims_no_ext.pdf", sims_fig_noext)






