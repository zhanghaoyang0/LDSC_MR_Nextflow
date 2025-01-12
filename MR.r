libs <- c("TwoSampleMR", "gsmr", "dplyr", "ggplot2", "ggsci")
lapply(libs, require, character.only = TRUE)

args <- commandArgs(trailingOnly = TRUE)
path_gwas1 <- args[1]
path_gwas2 <- args[2]
path_iv <- args[3]
path_ld <- args[4]
path_out <- args[5]

# path_gwas1 = 'data/trait1.txt.gz'
# path_gwas2 = 'data/trait2.txt.gz'
# path_iv = 'clump/trait1.iv'
# path_ld = 'clump/trait1.ld'
# path_out = "result/trait1_TO_trait2"

# read gwas
print('loading data...')
df1_raw <- read.table(path_gwas1, header = 1, sep = "\t")
df2_raw <- read.table(path_gwas2, header = 1, sep = "\t")
snp <- unlist(read.table(path_iv))
print('done !')

# subset iv
print('running mr ...')
df1 <- df1_raw %>%
    filter(SNP %in% snp) %>%
    select(-CHR, -POS)
df2 <- df2_raw %>%
    filter(SNP %in% snp) %>%
    select(-CHR, -POS)

# harmonise
df1 <- df1 %>%
    mutate(id.exposure = "exposure", exposure = "exposure") %>%
    rename("pval.exposure" = "P", "effect_allele.exposure" = "A1", "other_allele.exposure" = "A2", "samplesize.exposure" = "N", "beta.exposure" = "BETA", "se.exposure" = "SE", "eaf.exposure" = "FRQ")
df2 <- df2 %>%
    mutate(id.outcome = "outcome", outcome = "outcome") %>%
    rename("pval.outcome" = "P", "effect_allele.outcome" = "A1", "other_allele.outcome" = "A2", "samplesize.outcome" = "N", "beta.outcome" = "BETA", "se.outcome" = "SE", "eaf.outcome" = "FRQ")
dat <- harmonise_data(df1, df2) %>% filter(mr_keep == T)

# ivw, egger, weighted median, weighted mode
method_list <- c("mr_wald_ratio", "mr_egger_regression", "mr_weighted_median", "mr_ivw", "mr_weighted_mode")
mr_out <- mr(dat, method_list = method_list)
mr_out <- mr_out %>% select(method, nsnp, b, se, pval)
print('done!')

# gsmr
print('running gsmr ...')
get_gsmr_para <- function(pop = "eas") {
    n_ref <<- ifelse(pop == "eas", 481, 489) # Sample size of the 1000g eas (nrow of fam)
    gwas_thresh <<- 5e-8 # GWAS threshold to select SNPs as the instruments for the GSMR analysis
    single_snp_heidi_thresh <<- 0.01 # p-value threshold for single-SNP-based HEIDI-outlier analysis | default is 0.01
    multi_snp_heidi_thresh <<- 0.01 # p-value threshold for multi-SNP-based HEIDI-outlier analysis | default is 0.01
    nsnps_thresh <<- 5 # the minimum number of instruments required for the GSMR analysis | default is 10
    heidi_outlier_flag <<- T # flag for HEIDI-outlier analysis
    ld_r2_thresh <<- 0.05 # LD r2 threshold to remove SNPs in high LD
    ld_fdr_thresh <<- 0.05 # FDR threshold to remove the chance correlations between the SNP instruments
}
get_gsmr_para("eas")

# read ld matrix, for gsmr
ld <- read.table(path_ld)
colnames(ld) <- rownames(ld) <- read.table(gsub("ld", "snplist", path_ld))[, 1]
ldrho <- ld[rownames(ld) %in% dat$SNP, colnames(ld) %in% dat$SNP]
snp_coeff_id <- rownames(ldrho)

gsmr_out <- gsmr(dat$beta.exposure, dat$se.exposure, dat$pval.exposure, dat$beta.outcome, dat$se.outcome, dat$pval.outcome,
    ldrho, snp_coeff_id, n_ref, heidi_outlier_flag,
    gwas_thresh = gwas_thresh, single_snp_heidi_thresh, multi_snp_heidi_thresh, nsnps_thresh, ld_r2_thresh, ld_fdr_thresh
)
gsmr_out <- c("GSMR", length(gsmr_out$used_index), gsmr_out$bxy, gsmr_out$bxy_se, gsmr_out$bxy_pval)
print('done!')

# merge result
res <- rbind(mr_out, gsmr_out)
res
write.csv(res, paste0(path_out, "_MR.csv"), row.names=F)

# pleiotropy
pleio <- mr_pleiotropy_test(dat)
# heterogeneity
hetero <- mr_heterogeneity(dat)

# weak instruments bias
get_f <- function(dat) {
    n <- dat$samplesize.exposure[1]
    k <- nrow(dat)
    r <- get_r_from_bsen(dat$beta.exposure, dat$se.exposure, dat$samplesize.exposure)
    f <- (n - k - 1) * (sum(r^2)) / (1 - sum(r^2)) / k
    return(f)
}
get_f(dat)
# 74.38804

print('ploting ...')
res <- res %>% mutate(b = as.numeric(b))
res$a <- ifelse(res$method == "MR Egger", pleio$egger_intercept, 0) # I only add intercept for MR Egger since other methods have not provided intercepts
p <- ggplot(data = dat, aes(x = beta.exposure, y = beta.outcome)) +
    geom_errorbar(aes(ymin = beta.outcome - se.outcome, ymax = beta.outcome + se.outcome), colour = "grey", width = 0) +
    geom_errorbarh(aes(xmin = beta.exposure - se.exposure, xmax = beta.exposure + se.exposure), colour = "grey", height = 0) +
    geom_point(size = 3, shape = 16, fill = "blue", alpha = 0.6) +
    geom_abline(data = res, aes(intercept = a, slope = b, colour = method), show.legend = TRUE) +
    labs(colour = "Method", x = "GWAS effect on exposure", y = "GWAS effect on outcome") +
    theme_bw() +
    theme(
        legend.position = c(0.22, 0.83),
        legend.title = element_text(size = 11),
        legend.text = element_text(size = 10),
        panel.background = element_blank(),
        axis.text = element_text(size = 12, margin = margin(t = 10, r = 10, b = 10, l = 10)),
        axis.title = element_text(size = 14, margin = margin(t = 10, r = 10, b = 10, l = 10))
    ) +
    guides(colour = guide_legend(title.position = "top", ncol = 1)) +
    scale_fill_nejm()

png(paste0(path_out, "_MR_scatter.png"))
print(p)
dev.off()
print('done!')


write.csv(pleio, paste0(path_out, "_MR_pleiotropy.csv"), row.names=F)
write.csv(hetero, paste0(path_out, "_MR_heterogenity.csv"), row.names=F)
print('all analysis finish!')