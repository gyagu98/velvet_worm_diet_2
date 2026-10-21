# Trial 5 figures: BLAST-verified prey and what each sequence turned out to be.
# Usage: Rscript trial5_figures.R trial5_results
suppressPackageStartupMessages({library(ggplot2); library(dplyr); library(readr); library(tidyr)})
args <- commandArgs(trailingOnly = TRUE); d <- if (length(args)) args[1] else "trial5_results"
cat("Running Trial 5 figures. Reading:", file.path(d, "trial5_all_sequences.tsv"), "\n")
all <- read_tsv(file.path(d, "trial5_all_sequences.tsv"), show_col_types = FALSE)
out <- file.path(d, "figures"); dir.create(out, showWarnings = FALSE)

lv <- c("PREY", "host (velvet worm)", "no BLAST hit", "lab contamination", "implausible (marine)", "non-animal")
all$verdict <- factor(all$verdict, levels = intersect(lv, unique(all$verdict)))
th <- theme_classic(base_size = 11) + theme(legend.position = "bottom")

# Fig A: fate of every candidate COX1 sequence, per worm
pA <- ggplot(all, aes(sample, fill = verdict)) + geom_bar(colour = "grey20", linewidth = 0.2) +
  scale_fill_grey(start = 0.1, end = 0.9, name = NULL) +
  labs(x = NULL, y = "Candidate COX1 sequences", title = "BLAST verification of candidate prey sequences") + th
ggsave(file.path(out, "fig_A_sequence_fate.pdf"), pA, width = 6, height = 4.5)
ggsave(file.path(out, "fig_A_sequence_fate.png"), pA, width = 6, height = 4.5, dpi = 300)

# Fig B: verified prey orders per worm
prey <- all %>% filter(verdict == "PREY") %>%
  mutate(label = paste0(order, " (", class, ")"))
if (nrow(prey)) {
  pB <- prey %>% count(sample, label) %>%
    complete(sample = unique(all$sample), label, fill = list(n = 0)) %>%
    ggplot(aes(sample, label, fill = n > 0)) + geom_tile(colour = "white") +
    geom_text(aes(label = ifelse(n > 0, n, "")), colour = "white") +
    scale_fill_manual(values = c(`FALSE` = "grey92", `TRUE` = "grey20"), guide = "none") +
    labs(x = NULL, y = NULL, title = "BLAST-verified prey orders") + th
  ggsave(file.path(out, "fig_B_verified_prey.pdf"), pB, width = 5.5, height = 1.5 + 0.4 * n_distinct(prey$label))
  ggsave(file.path(out, "fig_B_verified_prey.png"), pB, width = 5.5, height = 1.5 + 0.4 * n_distinct(prey$label), dpi = 300)
}
cat("Done - figures in", out, "\n")
