# Velvet worm diet - Trial 4 figures (validated, de-duplicated prey sequences)
# Input : Trial_4/all_samples.validated_prey.tsv   (made by trial4_validate.py)
# Run   : setwd(".../folder that contains Trial_4");  source(file.choose())

suppressPackageStartupMessages({ library(tidyverse); library(patchwork) })
message("Running Trial 4 figures v2 (unique sequences, contaminants removed)")

args <- commandArgs(trailingOnly = TRUE)
res <- if (length(args) >= 1) args[1] else "results"
cands <- c(file.path(res, "all_samples.validated_prey.tsv"), "results/all_samples.validated_prey.tsv",
           "Trial_4/results/all_samples.validated_prey.tsv", "all_samples.validated_prey.tsv")
f <- cands[file.exists(cands)][1]
if (is.na(f)) stop("Cannot find all_samples.validated_prey.tsv - run this from the Trial_4 folder")
message("Reading: ", f)
d <- read_tsv(f, show_col_types = FALSE, col_types = cols(.default = "c")) %>%
  mutate(n_contigs_merged = as.integer(n_contigs_merged),
         blast_flag = replace_na(blast_flag, ""))

# drop anything BLAST says is velvet worm; conflicting order calls kept as their own row
d <- d %>% filter(!str_detect(blast_flag, "HOST")) %>%
  mutate(order = if_else(str_starts(order, "conflict"), "Unresolved (conflicting)", order),
         label = if_else(order == "Unresolved (conflicting)", order, paste0(order, " (", class, ")")),
         tier  = factor(tier, levels = c("strong", "standard")))

out <- file.path(dirname(f), "figures"); dir.create(out, showWarnings = FALSE, recursive = TRUE)

counts <- d %>% count(sample, label, tier, name = "n")
tot    <- counts %>% group_by(sample, label) %>% summarise(n = sum(n), .groups = "drop")
write_tsv(tot %>% pivot_wider(names_from = sample, values_from = n, values_fill = 0) %>%
            arrange(label), file.path(out, "table_unique_prey_sequences_by_order.tsv"))

tier_cols <- c(strong = "grey15", standard = "grey65")
theme_set(theme_classic(base_size = 10) +
            theme(plot.title = element_text(face = "bold", hjust = 0.5),
                  legend.position = "bottom"))

bar <- function(df, title) {
  ord <- df %>% group_by(label) %>% summarise(n = sum(n)) %>% arrange(n) %>% pull(label)
  ggplot(df, aes(n, factor(label, levels = ord), fill = tier)) +
    geom_col(width = 0.7) +
    scale_fill_manual(values = tier_cols, drop = FALSE,
                      labels = c(strong = ">=150 bp, >=95% id", standard = ">=80 bp, >=90% id"),
                      name = NULL) +
    scale_x_continuous(expand = expansion(mult = c(0, 0.1)), breaks = scales::breaks_pretty()) +
    labs(x = "Unique prey COX1 sequences", y = NULL, title = title)
}

per <- map(sort(unique(counts$sample)), function(s) {
  p <- bar(filter(counts, sample == s), s)
  ggsave(file.path(out, paste0("fig_", s, ".pdf")), p, width = 5, height = 4)
  ggsave(file.path(out, paste0("fig_", s, ".png")), p, width = 5, height = 4, dpi = 300)
  p
})

pooled <- counts %>% group_by(label, tier) %>% summarise(n = sum(n), .groups = "drop")
p_all <- bar(pooled, "All worms")
ggsave(file.path(out, "fig_all_worms.pdf"), p_all, width = 5.5, height = 5)

pa <- tot %>% complete(sample, label, fill = list(n = 0)) %>% mutate(present = n > 0)
p_pa <- ggplot(pa, aes(sample, label)) +
  geom_tile(aes(fill = present), colour = "white") +
  geom_text(aes(label = if_else(n > 0, as.character(n), "")), size = 3, colour = "white") +
  scale_fill_manual(values = c(`TRUE` = "grey20", `FALSE` = "grey92"), guide = "none") +
  labs(x = NULL, y = NULL, title = "Prey orders detected per worm")
ggsave(file.path(out, "fig_presence_absence.pdf"), p_pa, width = 4.5, height = 5)

panel <- (wrap_plots(per, nrow = 1) / (p_all | p_pa)) +
  plot_layout(guides = "collect") + plot_annotation(tag_levels = "a") &
  theme(legend.position = "bottom")
ggsave(file.path(out, "fig_all_in_one.pdf"), panel, width = 13, height = 10)
ggsave(file.path(out, "fig_all_in_one.png"), panel, width = 13, height = 10, dpi = 300)

message("Done - figures in ", out)
