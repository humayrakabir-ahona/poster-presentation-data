###############################################################################
# Volcano plot: MM1R vs MM1S, untreated
# Differential expression — log2FC vs -log10(adjusted p)
#
# Requires: ggplot2, ggrepel, dplyr
#   install.packages(c("ggplot2", "ggrepel", "dplyr"))
###############################################################################

library(ggplot2)
library(ggrepel)
library(dplyr)

## ---------------------------------------------------------------------------
## 1. DATA
## ---------------------------------------------------------------------------
# Replace this block with your real DE results. You need a data.frame with:
#   gene       : gene symbol (character)
#   log2FC     : log2 fold change, MM1R vs MM1S (numeric)
#   padj       : BH-adjusted p-value (numeric)
#   highlight  : optional logical column flagging a special gene set to draw
#                as open circles (e.g. a curated panel, housekeeping genes,
#                genes tested by qPCR) — set to FALSE for all rows if unused
#
# e.g. df <- read.csv("your_DE_results.csv")

set.seed(1)
n_genes <- 4000

# Simulated example data only — calibrated to resemble the reported shape
# (long tails, several highly significant genes on each side, a small set
# of near-zero "highlight" genes). DELETE this block once you plug in real
# results above.
log2FC <- c(rnorm(n_genes, mean = 0, sd = 0.6))
padj   <- 10^(-abs(rnorm(n_genes, mean = 3, sd = 3)) * (0.3 + abs(log2FC)))
padj   <- pmin(padj, 1)

df <- data.frame(
  gene   = paste0("GENE", seq_len(n_genes)),
  log2FC = log2FC,
  padj   = padj
)

# Inject a handful of strongly significant genes at specific positions so the
# labeled set below has something to point to — remove once using real data.
named_hits <- tribble(
  ~gene,       ~log2FC, ~padj,
  "NLRP11",     4.6,     10^-121,
  "ADGRE5",     4.3,     10^-118,
  "NES",        4.7,     10^-113,
  "TUBB4A",     3.0,     10^-75,
  "PTPRCAP",    5.6,     10^-71,
  "BTK",        2.7,     10^-69,
  "CD52",       2.9,     10^-46,
  "LINC01518",  5.6,     10^-43,
  "GNG7",       1.9,     10^-42,
  "LHX8",      -2.3,     10^-96,
  "GNE",       -2.0,     10^-65,
  "MBP",       -5.7,     10^-46,
  "NR3C1",     -1.7,     10^-46,
  "VAV2",      -1.8,     10^-42,
  "JAKMIP1",   -1.4,     10^-42
)
df <- df %>%
  filter(!gene %in% named_hits$gene) %>%
  bind_rows(named_hits)

# A small "highlight" panel plotted as open circles near the origin
# (e.g. a set of control/housekeeping genes tested independently).
highlight_genes <- data.frame(
  gene   = paste0("CTRL", 1:8),
  log2FC = c(0.05, -0.05, 0.1, -0.1, 0.15, -0.15, 0.2, -0.2),
  padj   = runif(8, 0.05, 0.9)
)
highlight_genes$highlight <- TRUE
df$highlight <- FALSE
df <- bind_rows(df, highlight_genes)

genes_to_label <- named_hits$gene

## ---------------------------------------------------------------------------
## 2. THRESHOLDS & COLOR ASSIGNMENT
## ---------------------------------------------------------------------------
fc_cutoff <- 1        # |log2FC| threshold (dashed vertical lines)
p_cutoff  <- 0.05      # adjusted-p threshold (dashed horizontal line)

df <- df %>%
  mutate(
    neglog10p = -log10(pmax(padj, 1e-300)),
    sig_class = case_when(
      highlight                                  ~ "highlight",
      padj < p_cutoff & log2FC >=  fc_cutoff      ~ "up",
      padj < p_cutoff & log2FC <= -fc_cutoff      ~ "down",
      TRUE                                        ~ "ns"
    )
  )

color_map <- c(
  up        = "#E07B1A",  # orange
  down      = "#3E7CB1",  # blue
  ns        = "grey75",
  highlight = "black"
)

## ---------------------------------------------------------------------------
## 3. PLOT
## ---------------------------------------------------------------------------
label_df <- df %>% filter(gene %in% genes_to_label)

p <- ggplot(df %>% filter(sig_class != "highlight"),
            aes(x = log2FC, y = neglog10p, color = sig_class)) +
  geom_point(size = 2, alpha = 0.85) +
  geom_point(
    data = df %>% filter(sig_class == "highlight"),
    aes(x = log2FC, y = neglog10p),
    shape = 21, size = 3, stroke = 1.1, color = "black", fill = NA
  ) +
  geom_vline(xintercept = c(-fc_cutoff, fc_cutoff),
             linetype = "dashed", color = "grey40", linewidth = 0.5) +
  geom_hline(yintercept = -log10(p_cutoff),
             linetype = "dashed", color = "grey40", linewidth = 0.5) +
  geom_text_repel(
    data = label_df,
    aes(label = gene),
    color = "black", size = 5, fontface = "plain",
    max.overlaps = Inf,
    box.padding = 0.6,
    point.padding = 0.3,
    min.segment.length = 0,
    segment.color = "grey30"
  ) +
  scale_color_manual(values = color_map, guide = "none") +
  labs(
    x = NULL,
    y = expression(-log[10]~"adjusted"~italic(p))
  ) +
  theme_minimal(base_size = 18) +
  theme(
    panel.grid.minor   = element_blank(),
    panel.grid.major   = element_line(color = "grey90"),
    axis.line          = element_blank(),
    axis.text          = element_text(color = "black", size = 16),
    axis.title.y       = element_text(size = 20, margin = margin(r = 10)),
    plot.margin        = margin(t = 15, r = 20, b = 15, l = 15)
  )

## ---------------------------------------------------------------------------
## 4. SAVE
## ---------------------------------------------------------------------------
ggsave("volcano_MM1R_vs_MM1S_untreated.png",
       plot = p, width = 13, height = 10, dpi = 300, bg = "white")

print(p)
