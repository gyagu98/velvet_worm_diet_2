# Trial 5 figures v2 - manuscript set (Biology Letters)
# Usage (from the Trial_5 folder):
#   Rscript trial5_figures_v2.R trial5_results ../Trial_4/results
# Arg 1 = Trial 5 results folder (has trial5_all_sequences.tsv)
# Arg 2 = Trial 4 results folder (has all_samples.validated_prey.tsv) - optional,
#         needed for Fig 1 (database-only call vs BLAST verdict).
suppressPackageStartupMessages({
  library(ggplot2); library(dplyr); library(readr); library(tidyr); library(patchwork)
})
cat("Running Trial 5 figures v2\n")
args <- commandArgs(trailingOnly = TRUE)
d5 <- if (length(args) >= 1) args[1] else "trial5_results"
d4 <- if (length(args) >= 2) args[2] else "../Trial_4/results"
out <- file.path(d5, "figures_v2"); dir.create(out, showWarnings = FALSE)

all <- read_tsv(file.path(d5, "trial5_all_sequences.tsv"), show_col_types = FALSE,
                col_types = cols(.default = "c")) %>%
  mutate(across(c(pident, aln_len, bitscore), as.numeric))

lv <- c("PREY", "host (velvet worm)", "no BLAST hit", "lab contamination",
        "implausible (marine)", "non-animal")
lab <- c("PREY" = "Verified prey", "host (velvet worm)" = "Host (Onychophora)",
         "no BLAST hit" = "No BLAST hit", "lab contamination" = "Lab contaminant",
         "implausible (marine)" = "Implausible (marine)", "non-animal" = "Non-animal")
# Okabe-Ito colour-blind safe palette
pal <- c("Verified prey" = "#D55E00", "Host (Onychophora)" = "#0072B2",
         "No BLAST hit" = "#BBBBBB", "Lab contaminant" = "#CC79A7",
         "Implausible (marine)" = "#56B4E9", "Non-animal" = "#009E73")
all <- all %>% mutate(fate = factor(lab[verdict], levels = lab[lv])) %>%
  mutate(fate = droplevels(fate))
th <- theme_classic(base_size = 10) +
  theme(legend.position = "bottom", plot.title = element_text(face = "bold", size = 10),
        axis.text = element_text(colour = "black"))
save <- function(p, name, w, h) {
  ggsave(file.path(out, paste0(name, ".pdf")), p, width = w, height = h)
  ggsave(file.path(out, paste0(name, ".png")), p, width = w, height = h, dpi = 300)
}

# ---- Fig 1: database-only call (Trial 4) vs BLAST verdict (Trial 5) -------------
p1 <- NULL
f4 <- file.path(d4, "all_samples.validated_prey.tsv")
if (file.exists(f4)) {
  t4 <- read_tsv(f4, show_col_types = FALSE, col_types = cols(.default = "c")) %>%
    mutate(db_order = ifelse(grepl("^conflict", order), "Ambiguous", order))
  key <- function(x) sub(" .*", "", x)  # BLAST may shorten IDs at the first space
  j <- t4 %>% mutate(k = key(seq_id)) %>%
    left_join(all %>% mutate(k = key(seq)) %>% select(k, fate), by = "k")
  if (sum(!is.na(j$fate)) == 0) {
    cat("NOTE: Trial 4 IDs did not match Trial 5 IDs - Fig 1 skipped\n")
  } else {
    j <- j %>% mutate(fate = factor(ifelse(is.na(fate), "No BLAST hit", as.character(fate)),
                                    levels = levels(all$fate)))
    ord <- j %>% count(db_order) %>% arrange(n) %>% pull(db_order)
    p1 <- ggplot(j, aes(y = factor(db_order, levels = ord), fill = fate)) +
      geom_bar(colour = "grey20", linewidth = 0.2, width = 0.75) +
      scale_fill_manual(values = pal, name = NULL, drop = TRUE) +
      scale_x_continuous(expand = expansion(mult = c(0, 0.05))) +
      labs(x = "Unique COX1 sequences", y = "Order assigned by MIDORI2 matching",
           title = "Database-only prey calls vs. BLAST verification") + th
    save(p1, "fig1_database_vs_blast", 6.5, 4.5)
  }
} else cat("NOTE: Trial 4 table not found at", f4, "- Fig 1 skipped\n")

# ---- Fig 2: fate of every candidate sequence per worm (with counts) -------------
cnt <- all %>% count(sample, fate)
tot <- all %>% count(sample, name = "N")
p2 <- ggplot(cnt, aes(sample, n, fill = fate)) +
  geom_col(colour = "grey20", linewidth = 0.2, width = 0.7) +
  geom_text(aes(label = ifelse(n >= 2, n, ""), group = fate), position = position_stack(vjust = 0.5), size = 3,
            colour = ifelse(cnt$fate %in% c("Host (Onychophora)", "Verified prey"), "white", "black")) +
  geom_text(data = tot, aes(sample, N, label = paste0("n = ", N)), inherit.aes = FALSE,
            vjust = -0.5, size = 3) +
  scale_fill_manual(values = pal, name = NULL) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.1))) +
  labs(x = "Individual", y = "Candidate COX1 sequences", title = "Fate of candidate sequences") + th
save(p2, "fig2_sequence_fate", 4.5, 4.5)

# ---- Fig 3: identity vs alignment length of best BLAST hit ----------------------
hits <- all %>% filter(!is.na(pident))
p3 <- ggplot(hits, aes(aln_len, pident, colour = fate, shape = sample)) +
  annotate("rect", xmin = 150, xmax = Inf, ymin = 95, ymax = Inf, alpha = 0.08, fill = "black") +
  geom_hline(yintercept = c(90, 95), linetype = c("dashed", "dotted"), colour = "grey40") +
  geom_vline(xintercept = c(80, 150), linetype = c("dashed", "dotted"), colour = "grey40") +
  geom_point(size = 2.4, alpha = 0.85) +
  scale_colour_manual(values = pal, guide = "none") + scale_shape_manual(values = c(16, 17, 15, 18), name = NULL) +
  labs(x = "Alignment length of best hit (bp)", y = "Identity to best hit (%)",
       title = "Best BLAST hit per sequence",
       caption = "Dashed: order-level threshold (80 bp, 90%). Dotted/shaded: family-level (150 bp, 95%).") +
  th + theme(plot.caption = element_text(size = 7, hjust = 0))
save(p3, "fig3_identity_vs_length", 5.5, 4.5)

# ---- Fig 4: verified prey detail -------------------------------------------------
prey <- all %>% filter(verdict == "PREY")
samples <- sort(unique(all$sample))
grid <- expand.grid(sample = samples, stringsAsFactors = FALSE) %>%
  left_join(prey %>% select(sample, class, order, family, top_hit_name, pident, aln_len, report_level),
            by = "sample") %>%
  mutate(txt = ifelse(is.na(order), "No verified prey",
    paste0(order, if_else(report_level == "family" & !is.na(family) & family != "",
                          paste0(": ", family), ""),
           "\nbest hit ", top_hit_name, "\n", round(pident, 1), "% over ", aln_len, " bp")),
    has = !is.na(order))
p4 <- ggplot(grid, aes(x = 1, y = factor(sample, levels = rev(samples)))) +
  geom_tile(aes(fill = has), colour = "white", linewidth = 1.5, width = 1) +
  geom_text(aes(label = txt, colour = has), size = 2.8, lineheight = 0.95) +
  scale_fill_manual(values = c(`TRUE` = "#D55E00", `FALSE` = "grey93"), guide = "none") +
  scale_colour_manual(values = c(`TRUE` = "white", `FALSE` = "grey40"), guide = "none") +
  labs(x = NULL, y = NULL, title = "BLAST-verified prey") + th +
  theme(axis.text.x = element_blank(), axis.ticks = element_blank(), axis.line = element_blank())
save(p4, "fig4_verified_prey", 4, 3)

# ---- Composite main figure -------------------------------------------------------
main <- if (!is.null(p1)) {
  (p1 | (p2 / p4 + plot_layout(heights = c(1.4, 1)))) / p3 +
    plot_layout(heights = c(1.1, 1), guides = "collect") +
    plot_annotation(tag_levels = "a") & theme(legend.position = "bottom")
} else {
  (p2 | p4) / p3 + plot_layout(guides = "collect") +
    plot_annotation(tag_levels = "a") & theme(legend.position = "bottom")
}
save(main, "fig_main_composite", 10, 10)

# ---- Supplementary table ---------------------------------------------------------
write_tsv(all %>% select(sample, seq, fate, top_hit, top_hit_name, pident, aln_len, bitscore,
                         phylum, class, order, family, report_level),
          file.path(out, "table_S1_all_sequences.tsv"))
write_tsv(all %>% count(sample, fate) %>% pivot_wider(names_from = fate, values_from = n, values_fill = 0),
          file.path(out, "table_1_summary.tsv"))
cat("Done - figures in", out, "\n")
