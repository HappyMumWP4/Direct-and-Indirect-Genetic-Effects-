# Install the required packages (only first time)
install.packages("bigsnpr")

# Load the required packages
library(bigsnpr)


#########################################
# LD MATRIX  AND SET OF VARIANTS TO USE #
#########################################

# Obtain HapMap3 SNPs and LD correlation matrix downloadable at https://ndownloader.figshare.com/files/25503788
# Use https://figshare.com/ndownloader/files/37802721 for HapMap3+ variants
info <- readRDS("/path/to/map.rds") # retrieved from https://figshare.com/articles/dataset/European_LD_reference/13034123?file=25503788


#########################################
# LOAD AND TRANSFORM SUMMARY STATISTICS #
#########################################

# Read in the summary statistic file
sumstats_all <- bigreadr::fread2("/path/to/summary_statistics.txt") 

# Select the needed variables
# LDpred 2 requires the headers with the following exact naming
# Modify the name of sumstats_all$columns in order to match the required info
sumstats <- data.frame(chr= sumstats_all$CHROM,
                       pos= sumstats_all$POS,
                       rsid= sumstats_all$ID,
                       a0= sumstats_all$ALT,
                       a1= sumstats_all$REF, #effect allele for which BETA/OR are calculated
                       # n_eff= sumstats_all$NEFFDIV2,
                       n_case= sumstats_all$NCAS,
                       n_control= sumstats_all$NCON,
                       beta_se= sumstats_all$SE,
                       p= sumstats_all$PVAL,
                       beta= sumstats_all$BETA
                       # OR= sumstats_all$OR
                       # INFO= sumstats_all$INFO
)


# If OR is provided instead of beta, transform the OR into log(OR)
sumstats$beta <- log(sumstats$OR)

# In case of binary trait with separate N for cases and controls
sumstats$n_eff <- 4 / (1 / sumstats$n_case + 1 / sumstats$n_control)
sumstats$n_case <- sumstats$n_control <- NULL

# Filter out HapMap3 SNPs
sumstats <- sumstats[sumstats$rsid %in% info$rsid,] #ne ha perse un sacco


##########################################################
# MATCH VARIANTS BETWEEN GENOTYPE AND SUMMARY STATISTICS #
##########################################################

# Extract the SNP information from the genotype
map <- data.frame(chr= info$chr,
                  rsid= info$rsid,
                  pos= info$pos,
                  a1= info$a1,
                  a0= info$a0)

# Perform SNP matching
df_beta <- snp_match(sumstats, map) # Error: not enough variants have been matched --> possibly due to different genome builds. NB: this function is case sensitive
df_beta <- snp_match(sumstats, map, join_by_pos = FALSE)  # use rsid instead of pos as a solution to previous error
# If the error persists, try:
map_mod<- snp_modifyBuild(info_snp=map,
                          liftOver='/path/to/liftOver', # retrieved from: https://privefl.github.io/bigsnpr/reference/snp_modifyBuild.html
                          from = "hg18",
                          to = "hg19",
                          check_reverse = TRUE)
df_beta <- snp_match(sumstats, map_mod)
# Alternatively, according to the documentation, try:
df_beta <- snp_match(sumstats, map, join_by_pos = FALSE, match.min.prop=0.05)


###################################################
# IF YOU WANT TO PERFORM QC OF SUMMARY STATISTICS #
###################################################

# info_snp <- tidyr::drop_na(tibble::as_tibble(info))
info_snp<- info[info$rsid %in% df_beta$rsid,]

# Better to use af of GWAS and INFO scores as well (then can use 0.7 instead of 0.5 in L35)
sd_ldref <- with(info_snp, sqrt(2 * af_UKBB * (1 - af_UKBB)))

# IF BINARY TRAIT
sd_ss <- with(df_beta, 2 / sqrt(n_eff * beta_se^2 + beta^2))

# IF CONTINUOUS TRAIT
sd_y = with(df_beta, sqrt(quantile(0.5 * (n_eff * beta_se^2 + beta^2), 0.01)))
sd_ss = with(df_beta, sd_y / sqrt(n_eff * beta_se^2 + beta^2))

# Estimate SNPs to remove
is_bad <- sd_ss < (0.5 * sd_ldref) | sd_ss > (sd_ldref + 0.1) | sd_ss < 0.05 | sd_ldref < 0.05

library(ggplot2)
qplot(sd_ldref, sd_ss, color = is_bad, alpha = I(0.5)) +
  theme_bigstatsr() +
  coord_equal() +
  scale_color_viridis_d(direction = -1) +
  geom_abline(linetype = 2, color = "red") +
  labs(x = "SD from allele frequencies of the LD reference",
       y = "SD from summary statistics",
       color = "Removed?")

# qplot(sd_af, sd_ss, color = is_bad, alpha = I(0.5)) +
#   theme_bigstatsr() +
#   coord_equal() +
#   scale_color_viridis_d(direction = -1) +
#   geom_abline(linetype = 2, color = "red") +
#   labs(x = "Standard deviations derived from allele frequencies",
#        y = "Standard deviations derived from the summary statistics",
#        color = "Removed?")

df_beta_no_qc <- df_beta
df_beta <- df_beta[!is_bad, ]


######################################
# COMPUTE LDpred2 SCORES GENOME-WIDE #
######################################

# Create a correlation matrix
NCORES <- nb_cores()
tmp <- tempfile(tmpdir = "/path/to/new/folder/tmp-data")

for (chr in 1:22) {
  
  cat(chr, ".. ", sep = "")
  
  ## indices in 'df_beta'
  ind.chr <- which(df_beta$chr == chr)
  ## indices in 'map_ldref'
  ind.chr2 <- df_beta$`_NUM_ID_`[ind.chr]
  ## indices in 'corr_chr'
  ind.chr3 <- match(ind.chr2, which(map$chr == chr))
  
  corr_chr <- readRDS(paste0("/path/to/folder/LD_chr", chr, ".rds"))[ind.chr3, ind.chr3] # "LD_chr" folder includes LD blocks retrieved from: https://figshare.com/articles/dataset/European_LD_reference_with_blocks_/19213299/1?file=34133082
  
  if (chr == 1) {
    corr <- as_SFBM(corr_chr, tmp, compact = TRUE)
  } else {
    corr$add_columns(corr_chr, nrow(corr))
  }
}

# Select LD of the corresponding SNPs
ld <- info[info$rsid %in% df_beta$rsid,]

# Estimate h2 from LD Score regression
ldsc <- with(df_beta, snp_ldsc(ld$ld, length(ld$ld), chi2 = (beta / beta_se)^2,
                               sample_size = n_eff, blocks = NULL))
ldsc_h2_est <- ldsc[["h2"]]



############################################################
# ESTIMATE ADJUSTED BETAS WITH THE LDpred2 AUTOMATIC MODEL #
############################################################

coef_shrink <- 0.95  # Reduce this up to 0.4 if you have some (large) mismatch with the LD ref

set.seed(42)  # to get the same result every time

# Can take minutes to hours
multi_auto <- snp_ldpred2_auto(corr,
                               df_beta,
                               h2_init = ldsc_h2_est,
                               vec_p_init = seq_log(1e-4, 0.2, length.out = 30), # 0.2 as tutorial, 0.9 as in the paper
                               ncores = NCORES,
                               # use_MLE = FALSE,  # uncomment if you have convergence issues or when power is low (need v1.11.9)
                               allow_jump_sign = FALSE, shrink_corr = coef_shrink)

# Verify whether the chains “converged” by looking at the path of the chains
library(ggplot2)
auto <- multi_auto[[1]]  # first chain
plot_grid(
  qplot(y = auto$path_p_est) + 
    theme_bigstatsr() + 
    geom_hline(yintercept = auto$p_est, col = "blue") +
    scale_y_log10() +
    labs(y = "p"),
  qplot(y = auto$path_h2_est) + 
    theme_bigstatsr() + 
    geom_hline(yintercept = auto$h2_est, col = "blue") +
    labs(y = "h2"),
  ncol = 1, align = "hv"
)

# In the LDpred2 paper, we proposed an automatic way of filtering bad chains by comparing the scale of the resulting predictions.
# Here we recommend an equivalent and simpler alternative:
range <- sapply(multi_auto, function(auto) diff(range(auto$corr_est)))
keep <- which(range > (0.95 * quantile(range, 0.95, na.rm = TRUE)))

# To get the final effects you should only use chains that pass this filtering
beta_auto_means <- rowMeans(sapply(multi_auto[keep], function(auto) auto$beta_est))

# Save adjusted betas
adjusted_betas <- data.frame(CHR= df_beta$chr,
                             POS= df_beta$pos,
                             RSID = df_beta$rsid,
                             A0= df_beta$a0,
                             A1= df_beta$a1,
                             BETA_ADJ = beta_auto_means)
adjusted_betas_reverse <- data.frame(CHR= df_beta$chr,
                                     POS= df_beta$pos,
                                     RSID = df_beta$rsid,
                                     A1= df_beta$a1,
                                     A2= df_beta$a0,
                                     BETA_ADJ = beta_auto_means)

# For PLINK format, we use "adjusted beta reverse" because PLINK assumes A1 as effect allele and A2 as other allele
write.table(adjusted_betas_reverse, "/path/to/adjusted_weights.txt", sep="\t", dec =".", quote= FALSE, row.names=FALSE, col.names= FALSE)

