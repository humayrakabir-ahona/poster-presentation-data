#Ferroptosis_Heatmap_Rscript_by_humayra.kabir.ahona
# ============================================================
# Figure (REAL DATA VERSION)
# "Ferroptosis-associated transcriptional landscape across
#  (glucocorticoid-resistant and glucocorticoid-sensitive
#  multiple myeloma cells) cancer phenotyoes."
#
# Source data: GSE182638_norm_counts_TPM_GRCh38_p13_NCBI.tsv
#   (user-uploaded, NCBI GEO-style TPM matrix, NCBI GeneID rows,
#    18 samples: GSM5534024-GSM5534041)
#
# IMPORTANT SCOPE / HONESTY NOTES:
# 1. GSE182638 does NOT include an explicit ferroptosis-sensitivity
#    annotation for MM1S/MM1R. This figure shows transcriptional
#    differences ALONG A FERROPTOSIS-RELEVANT GENE SET between the
#    two glucocorticoid-phenotype groups. It is hypothesis-generating
#    for ferroptosis vulnerability -- it does NOT demonstrate or prove
#    ferroptosis sensitivity/resistance itself.
# 2. The series matrix supplied has no embedded sample-characteristics
#    header (GSM ID columns only, no metadata rows), and this sandbox
#    cannot reach NCBI GEO's sample-characteristics page (blocked by
#    reCAPTCHA / robots.txt from this environment). Sample identity
#    (MM1S vs MM1R) was therefore recovered directly from the
#    expression data using NR3C1 (the glucocorticoid receptor, GR):
#    MM1S is GR-positive and glucocorticoid-sensitive; MM1R is a
#    GR-negative subclone and is glucocorticoid-resistant (this is
#    the defining, literature-established distinction between the two
#    lines, not an arbitrary marker choice). NR3C1 TPM in this dataset
#    splits cleanly at the largest gap between consecutive sorted
#    values (~1-5 TPM vs ~8-25 TPM; a ~2x jump at the gap point vs
#    <15% steps elsewhere), giving an unambiguous 9-vs-9 split that
#    matches the expected two-cell-line design.
#    ** If you have the original GEO sample-characteristics table,
#    swap in the true labels below -- they should supersede this
#    inferred split. **
# 3. Ferroptosis gene set is fixed a priori from established
#    literature/FerrDb core regulators (Dixon 2012; Yang 2014;
#    Doll 2019; Zhou & Bao 2020 FerrDb), spanning the cystine/GSH,
#    lipid-peroxidation, and iron-handling axes -- NOT chosen post hoc
#    from what looked interesting in this dataset.
# ============================================================

suppressPackageStartupMessages({
  library(pheatmap)
  library(RColorBrewer)
})

tpm_path <- "/home/claude/GSE182638_TPM.tsv"
tpm <- read.delim(tpm_path, header = TRUE, check.names = FALSE, row.names = 1)

## ---- 1. A priori ferroptosis gene set (Entrez GeneID -> symbol) ----
gene_set <- data.frame(
  entrez = c(23657, 2879, 2182, 84883, 7037, 2495, 2512, 8031, 10162),
  symbol = c("SLC7A11","GPX4","ACSL4","AIFM2","TFRC","FTH1","FTL","NCOA4","LPCAT3"),
  axis   = c("Cystine/GSH","Lipid defense","Lipid remodel","Lipid defense",
             "Iron uptake","Iron storage","Iron storage","Ferritinophagy","Lipid remodel"),
  stringsAsFactors = FALSE
)

stopifnot(all(as.character(gene_set$entrez) %in% rownames(tpm)))

expr <- as.matrix(tpm[as.character(gene_set$entrez), , drop = FALSE])
rownames(expr) <- gene_set$symbol

## ---- 2. Recover MM1S (GC-sensitive) vs MM1R (GC-resistant) via NR3C1/GR ----
NR3C1_ID <- "2908"
stopifnot(NR3C1_ID %in% rownames(tpm))
gr <- as.numeric(tpm[NR3C1_ID, ])
names(gr) <- colnames(tpm)

## simple, transparent split: rank samples by NR3C1 (GR) expression and
## split into the lower and upper half (9 vs 9), matching the expected
## two-cell-line design. The split point (4.70 -> 8.40 TPM) is a >2-fold
## jump -- well outside the <15% step-to-step variation seen elsewhere
## in the sorted values -- so this is not an arbitrary bisection, it
## lands on the same boundary the data itself shows most clearly.
n <- length(gr)
rnk <- rank(gr, ties.method = "first")
phenotype <- ifelse(rnk <= n / 2, "Resistant (MM1R)", "Sensitive (MM1S)")
phenotype <- factor(phenotype, levels = c("Sensitive (MM1S)", "Resistant (MM1R)"))
names(phenotype) <- colnames(tpm)

cat("Group sizes:\n"); print(table(phenotype))
cat("\nGR (NR3C1) TPM by group:\n")
print(round(sort(gr), 2))

## ---- 3. Transform + scale for visualization ----
log_expr <- log2(expr + 1)
z <- t(scale(t(log_expr)))  # per-gene z-score across all 18 samples

## order columns by phenotype for a cleaner block layout, but still let
## pheatmap cluster within/across as it wants (cluster_cols = TRUE keeps
## the dendrogram informative; column order below only sets the input
## order pheatmap starts from)
ord <- order(phenotype)
z <- z[, ord]
phenotype_ord <- phenotype[ord]

annotation_col <- data.frame(Phenotype = phenotype_ord, row.names = colnames(z))
ann_colors <- list(
  Phenotype = c("Sensitive (MM1S)" = "#2C7FB8", "Resistant (MM1R)" = "#D95F02"),
  Axis = c("Cystine/GSH" = "#00BFC4", "Ferritinophagy" = "#F8766D",
           "Iron storage" = "#C77CFF", "Iron uptake" = "#7CAE00",
           "Lipid defense" = "#FF61C3", "Lipid remodel" = "#00E5EE")
)
annotation_row <- data.frame(Axis = gene_set$axis, row.names = gene_set$symbol)

heat_colors <- colorRampPalette(rev(brewer.pal(11, "RdBu")))(100)

## ---- 4. Render: square, poster-ready ----
out_path <- "/home/claude/figure2_real_ferroptosis_heatmap.png"
side_px <- 3000
png(out_path, width = side_px, height = side_px, res = 300, type = "cairo")

pheatmap(
  z,
  color = heat_colors,
  annotation_col = annotation_col,
  annotation_row = annotation_row,
  annotation_colors = ann_colors,
  cluster_rows = TRUE,
  cluster_cols = TRUE,
  clustering_method = "ward.D2",
  show_colnames = FALSE,
  fontsize_row = 13,
  fontsize_col = 10,
  fontsize = 9.5,
  cellwidth = 22,
  cellheight = 60,
  treeheight_row = 35,
  treeheight_col = 35,
  border_color = "white",
  legend = TRUE,
  main = "Ferroptosis-associated transcriptional landscape across\nglucocorticoid-resistant and glucocorticoid-sensitive MM cells",
  filename = NA
)

dev.off()
cat("\nSaved:", out_path, "\n")
