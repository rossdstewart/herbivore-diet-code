library(readxl)
library(janitor)
library(dplyr)
library(tidyr)
library(tibble)
library(ggplot2)
library(forcats)
library(vegan)

setwd("C:/Users/Admin/Desktop/R/Rose et al")
# -----------------------------
# 1) Settings
# -----------------------------
path_to_file <- "Darius_3_formatted_hit_May15_cleanedDS.xlsx"
sheet_name   <- "Darius_3_formatted_hit_May15_KR"

min_reads_per_sample  <- 100
min_samples_per_animal <- 3
min_reads_per_order   <- 100

# -----------------------------
# 2) Helpers
# -----------------------------
base_name <- function(x) gsub("[\\._].*$", "", x)

clean_animal_name <- function(x) {
  dplyr::recode(
    x,
    "bontebok" = "Bontebok",
    "buffalo" = "Buffalo",
    "zebra" = "Burchell's zebra",
    "capezebra" = "Cape zebra",
    "eland" = "Eland",
    "elephant" = "Elephant",
    "giraffe" = "Giraffe",
    "impala" = "Impala",
    "rhino" = "White Rhino",
    "springbok" = "Springbok",
    "waterbuck" = "Waterbuck",
    "wildebeest" = "Wildebeest",
    .default = x
  )
}

standardise_animal_name <- function(x) {
  dplyr::recode(
    x,
    "Burchell's zebra" = "Burchell's zebra",
    "White Rhino" = "White rhinoceros",
    .default = x
  )
}

normalise_animal_names <- function(x) {
  x %>%
    base_name() %>%
    clean_animal_name() %>%
    standardise_animal_name()
}

# -----------------------------
# 3) Animal size classes
# edit if you want different bins
# -----------------------------
size_df <- tibble(
  animal = c(
    "Bontebok", "Buffalo", "Cape zebra", "Eland", "Elephant",
    "Giraffe", "Impala", "Plains zebra", "Springbok",
    "Waterbuck", "White rhinoceros", "Wildebeest"
  ),
  size_class = c(
    "Medium", "Large", "Large", "Large", "Mega",
    "Mega", "Medium", "Large", "Medium",
    "Large", "Mega", "Large"
  )
)

# -----------------------------
# 4) Load data
# -----------------------------
raw <- read_excel(path_to_file, sheet = sheet_name) %>%
  clean_names()

if ("order" %in% names(raw)) {
  raw <- raw %>% rename(plant_order = order)
}

num_cols <- raw %>%
  select(where(is.numeric)) %>%
  names()

if (length(num_cols) == 0) stop("No numeric diet columns found.")

# -----------------------------
# 5) Drop low-read sample columns
# -----------------------------
col_totals <- raw %>%
  summarise(across(all_of(num_cols), ~ sum(.x, na.rm = TRUE))) %>%
  pivot_longer(everything(), names_to = "sample_col", values_to = "total_reads")

keep_num_cols <- col_totals %>%
  filter(total_reads >= min_reads_per_sample) %>%
  pull(sample_col)

raw <- raw %>%
  select(plant_order, all_of(keep_num_cols))

num_cols <- keep_num_cols

if (length(num_cols) == 0) stop("No sample columns remain after filtering.")

# -----------------------------
# 6) Keep animals with enough retained samples
# -----------------------------
animal_n <- tibble(col = num_cols) %>%
  mutate(animal = normalise_animal_names(col)) %>%
  count(animal, name = "n_samples")

keep_animals <- animal_n %>%
  filter(n_samples >= min_samples_per_animal) %>%
  pull(animal)

if (length(keep_animals) == 0) stop("No animals remain after filtering.")

# -----------------------------
# 7) Aggregate to ORDER x ANIMAL
# -----------------------------
ord_animal <- raw %>%
  filter(!is.na(plant_order), plant_order != "") %>%
  pivot_longer(cols = all_of(num_cols), names_to = "col", values_to = "value") %>%
  mutate(animal = normalise_animal_names(col)) %>%
  filter(animal %in% keep_animals) %>%
  group_by(plant_order, animal) %>%
  summarise(value = sum(value, na.rm = TRUE), .groups = "drop")

# -----------------------------
# 8) Drop low-abundance orders
# -----------------------------
keep_orders <- ord_animal %>%
  group_by(plant_order) %>%
  summarise(total_reads = sum(value), .groups = "drop") %>%
  filter(total_reads >= min_reads_per_order) %>%
  pull(plant_order)

ord_animal <- ord_animal %>%
  filter(plant_order %in% keep_orders)

# -----------------------------
# 9) Join size class and aggregate to ORDER x SIZE
# -----------------------------
ord_size <- ord_animal %>%
  left_join(size_df, by = "animal") %>%
  filter(!is.na(size_class)) %>%
  group_by(plant_order, size_class) %>%
  summarise(value = sum(value, na.rm = TRUE), .groups = "drop")

# -----------------------------
# 10) Relative composition within size class
# -----------------------------
ord_size_rel <- ord_size %>%
  group_by(size_class) %>%
  mutate(rel_abund = value / sum(value)) %>%
  ungroup()

# -----------------------------
# 11) Top orders for plotting
# -----------------------------
top_orders <- ord_size_rel %>%
  group_by(plant_order) %>%
  summarise(total_rel = sum(rel_abund), .groups = "drop") %>%
  slice_max(order_by = total_rel, n = 12, with_ties = FALSE) %>%
  pull(plant_order)

plot_dat <- ord_size_rel %>%
  mutate(
    plant_order = ifelse(plant_order %in% top_orders, plant_order, "Other")
  ) %>%
  group_by(size_class, plant_order) %>%
  summarise(rel_abund = sum(rel_abund), .groups = "drop") %>%
  mutate(
    size_class = factor(size_class, levels = c("Medium", "Large", "Mega"))
  )

# -----------------------------
# 12) Stacked bar plot
# -----------------------------
ggplot(plot_dat, aes(x = size_class, y = rel_abund, fill = plant_order)) +
  geom_col(position = "fill") +
  labs(
    x = "Animal size class",
    y = "Relative diet composition",
    fill = "Plant order",
    title = "Diet composition across animal size classes"
  ) +
  theme_bw(base_size = 12)

# -----------------------------
# 13) Heatmap version
# -----------------------------
heat_dat <- ord_size_rel %>%
  filter(plant_order %in% top_orders) %>%
  mutate(
    size_class = factor(size_class, levels = c("Medium", "Large", "Mega")),
    plant_order = fct_reorder(plant_order, rel_abund, .fun = max)
  )

ggplot(heat_dat, aes(x = size_class, y = plant_order, fill = rel_abund)) +
  geom_tile() +
  labs(
    x = "Animal size class",
    y = "Plant order",
    fill = "Relative abundance",
    title = "Order-level diet composition by animal size class"
  ) +
  theme_bw(base_size = 12)




library(dplyr)
library(tidyr)
library(tibble)
library(purrr)
library(ggplot2)

# ord_animal from your earlier pipeline
# rows = plant_order, columns = animals after filtering

mat_rel <- ord_animal %>%
  group_by(plant_order, animal) %>%
  summarise(value = sum(value, na.rm = TRUE), .groups = "drop") %>%
  pivot_wider(names_from = animal, values_from = value, values_fill = 0) %>%
  column_to_rownames("plant_order") %>%
  as.matrix()

mat_rel <- sweep(mat_rel, 2, colSums(mat_rel), "/")
mat_rel[is.na(mat_rel)] <- 0

ord_animal %>%
  count(plant_order) %>%
  filter(n > 1)


diet_long <- as.data.frame(mat_rel) %>%
  rownames_to_column("plant_order") %>%
  pivot_longer(-plant_order, names_to = "animal", values_to = "rel_abund") %>%
  left_join(size_df, by = "animal")

order_tests <- diet_long %>%
  group_by(plant_order) %>%
  summarise(
    p_value = kruskal.test(rel_abund ~ size_class)$p.value,
    med_medium = median(rel_abund[size_class == "Medium"], na.rm = TRUE),
    med_large  = median(rel_abund[size_class == "Large"], na.rm = TRUE),
    med_mega   = median(rel_abund[size_class == "Mega"], na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(p_adj = p.adjust(p_value, method = "BH")) %>%
  arrange(p_adj)

order_tests

##################################################
library(readr)
library(dplyr)
library(tidyr)
library(stringr)
library(tibble)

# PanTHERIA WR05 file
url_pantheria <- "https://esapubs.org/archive/ecol/E090/184/PanTHERIA_1-0_WR05_Aug2008.txt"

pantheria <- read_tsv(url_pantheria, show_col_types = FALSE)

# keep the columns we need
mass_raw <- pantheria %>%
  transmute(
    binomial = MSW05_Binomial,
    body_mass_g = `5-1_AdultBodyMass_g`
  ) %>%
  filter(!is.na(binomial), !is.na(body_mass_g))

mass_raw

target_species <- tribble(
  ~animal,               ~binomial,
  "Buffalo",             "Syncerus caffer",
  "Elephant",            "Loxodonta africana",
  "Giraffe",             "Giraffa camelopardalis",
  "Impala",              "Aepyceros melampus",
  "Burchell's zebra",    "Equus burchellii",
  "Springbok",           "Antidorcas marsupialis",
  "Waterbuck",           "Kobus ellipsiprymnus",
  "White rhinoceros",    "Ceratotherium simum",
  "Wildebeest",          "Connochaetes taurinus"
)

mass_df <- target_species %>%
  left_join(mass_raw, by = "binomial") %>%
  mutate(body_mass_kg = body_mass_g / 1000)

mass_df

mass_raw %>%
  filter(str_detect(binomial, "Syncerus|Loxodonta|Giraffa|Aepyceros|Equus|Kobus|Ceratotherium|Connochaetes")) %>%
  arrange(binomial)

diet_long_mass <- as.data.frame(mat_rel) %>%
  rownames_to_column("plant_order") %>%
  pivot_longer(-plant_order, names_to = "animal", values_to = "rel_abund") %>%
  left_join(mass_df %>% select(animal, body_mass_kg), by = "animal")

order_lm <- diet_long_mass %>%
  group_by(plant_order) %>%
  group_modify(~{
    fit <- lm(rel_abund ~ log10(body_mass_kg), data = .x)
    s <- summary(fit)
    tibble(
      slope = coef(fit)[2],
      p_value = coef(s)[2, 4],
      r2 = s$r.squared
    )
  }) %>%
  ungroup() %>%
  mutate(p_adj = p.adjust(p_value, method = "BH")) %>%
  arrange(p_adj)

order_lm

order_lm %>%
  slice_max(order_by = abs(slope), n = 12) %>%
  ggplot(aes(x = reorder(plant_order, slope), y = slope)) +
  geom_col() +
  coord_flip() +
  labs(x = NULL, y = "Slope vs log10(body mass)") +
  theme_bw()


library(dplyr)
library(ggplot2)
library(forcats)

library(dplyr)
library(ggplot2)
library(forcats)

plot_dat <- order_lm %>%
  slice_max(order_by = abs(slope), n = 12, with_ties = FALSE) %>%
  mutate(
    direction = ifelse(slope > 0, "Higher in larger animals", "Higher in smaller animals"),
    plant_order = fct_reorder(plant_order, slope),
    p_lab = paste0("p = ", signif(p_adj, 2)),
    sig = case_when(
      p_adj < 0.001 ~ "***",
      p_adj < 0.01  ~ "**",
      p_adj < 0.05  ~ "*",
      TRUE ~ ""
    )
  )

x_pad <- 0.08 * max(abs(plot_dat$slope), na.rm = TRUE)

ggplot(plot_dat, aes(x = plant_order, y = slope, fill = direction)) +
  geom_col(width = 0.8) +
  geom_hline(yintercept = 0, linetype = "dashed", linewidth = 0.6, colour = "grey35") +
  geom_text(
    aes(
      y = ifelse(slope > 0, slope + x_pad, slope - x_pad),
      label = paste0(sig, " ", p_lab)
    ),
    size = 3.5,
    hjust = ifelse(plot_dat$slope > 0, 0, 1)
  ) +
  coord_flip(clip = "off") +
  scale_fill_manual(
    values = c(
      "Higher in larger animals" = "#2C7FB8",
      "Higher in smaller animals" = "#D95F0E"
    )
  ) +
  labs(
    title = "Plant orders associated with herbivore body mass",
    subtitle = "Bars show regression coefficients; labels show BH-adjusted p-values",
    x = NULL,
    y = "Regression coefficient for log10(adult body mass)",
    fill = NULL
  ) +
  theme_bw(base_size = 13) +
  theme(
    panel.grid.major.y = element_blank(),
    panel.grid.minor = element_blank(),
    plot.title = element_text(face = "bold"),
    plot.subtitle = element_text(size = 11),
    legend.position = "top",
    plot.margin = margin(5.5, 40, 5.5, 5.5)
  )
###################################
#
X_mass <- t(mat_rel)

mass_df <- tibble(
  animal = rownames(X_mass)
) %>%
  left_join(mass_df, by = "animal")

dist_bray_mass <- vegdist(X_mass, method = "bray")

set.seed(123)
adon_mass <- adonis2(
  dist_bray_mass ~ log10(body_mass_kg),
  data = mass_df,
  permutations = 9999
)
adon_mass


# overall multivariate support stats from PERMANOVA
r2_lab <- round(adon_mass$R2[1], 3)
f_lab  <- round(adon_mass$F[1], 2)
p_lab  <- signif(adon_mass$`Pr(>F)`[1], 2)

subtitle_txt <- paste0(
  "Overall multivariate signal: PERMANOVA on Bray-Curtis, R² = ", r2_lab,
  ", F = ", f_lab, ", p = ", p_lab
)

plot_dat <- order_lm %>%
  slice_max(order_by = abs(slope), n = 12, with_ties = FALSE) %>%
  mutate(
    direction = ifelse(
      slope > 0,
      "Higher in larger animals",
      "Higher in smaller animals"
    ),
    direction = factor(
      direction,
      levels = c("Higher in larger animals", "Higher in smaller animals")
    ),
    plant_order = fct_reorder(plant_order, slope),
    r2_lab = paste0("R² = ", round(r2, 2))
  )

x_pad <- 0.02 * max(abs(plot_dat$slope), na.rm = TRUE)

ggplot(plot_dat, aes(x = plant_order, y = slope, fill = direction)) +
  geom_col(width = 0.8) +
  geom_hline(yintercept = 0, linetype = "dashed", linewidth = 0.6, colour = "grey35") +
  geom_text(
    aes(
      y = ifelse(slope > 0, slope + x_pad, 0.013),
      label = round(slope,3)
    ),
    size = 3.5,
    hjust = ifelse(plot_dat$slope > 0, 0, 1)
  ) +
  coord_flip(clip = "off") +
  scale_fill_manual(
    values = c(
      "Higher in larger animals" = "#2C7FB8",
      "Higher in smaller animals" = "#D95F0E"
    )
  ) +
  labs(
    title = "Plant orders associated with herbivore body mass",
    subtitle = subtitle_txt,
    x = NULL,
    y = "Regression coefficient for log10(adult body mass)",
    fill = NULL
  ) +
  theme_bw(base_size = 13) +
  theme(
    panel.grid.major.y = element_blank(),
    panel.grid.minor = element_blank(),
    plot.title = element_text(face = "bold"),
    plot.subtitle = element_text(size = 10.5),
    legend.position = "top",
    plot.margin = margin(5.5, 45, 5.5, 5.5)
  )
#
#
#
#
#
selected_order <- "Poales"

plot_one <- diet_long_mass %>%
  filter(plant_order == selected_order)

ggplot(plot_one, aes(x = log10(body_mass_kg), y = rel_abund)) +
  geom_point(size = 3) +
  geom_smooth(method = "lm", se = FALSE, linewidth = 0.9, colour = "#2C7FB8") +
  geom_text(aes(label = animal), nudge_y = 0.01, size = 3.5) +
  labs(
    title = paste(selected_order, "increases with herbivore body mass"),
    subtitle = "Points are species; line shows the fitted linear trend used to estimate the slope",
    x = "log10(adult body mass, kg)",
    y = paste("Relative abundance of", selected_order)
  ) +
  theme_bw(base_size = 13) +
  theme(
    panel.grid.minor = element_blank(),
    plot.title = element_text(face = "bold")
  )
#
# mat_rel must be plant_order x animal
# mass_df must have: animal, body_mass_kg

# 1) animal x order matrix
Y <- t(mat_rel)

# 2) make sure mass table matches
meta <- tibble(animal = rownames(Y)) %>%
  left_join(mass_df, by = "animal") %>%
  mutate(log_mass = log10(body_mass_kg))

# 3) pairwise diet distances
diet_dist <- as.matrix(vegdist(Y, method = "bray"))

# 4) pairwise body-mass differences
mass_diff <- outer(meta$log_mass, meta$log_mass, function(a, b) abs(a - b))

# 5) pull upper triangles into one table
pair_dat <- expand.grid(i = seq_len(nrow(Y)), j = seq_len(nrow(Y))) %>%
  filter(i < j) %>%
  transmute(
    animal_1 = meta$animal[i],
    animal_2 = meta$animal[j],
    body_mass_diff = mass_diff[cbind(i, j)],
    bray_curtis = diet_dist[cbind(i, j)]
  )

# 6) simple observed correlation
cor_test <- cor.test(
  pair_dat$body_mass_diff,
  pair_dat$bray_curtis,
  method = "spearman",
  exact = FALSE
)

cor_test
set.seed(123)

obs_r <- cor(pair_dat$body_mass_diff, pair_dat$bray_curtis, method = "spearman")

perm_r <- replicate(9999, {
  perm_mass <- sample(meta$log_mass)
  perm_diff <- outer(perm_mass, perm_mass, function(a, b) abs(a - b))
  
  perm_pair <- expand.grid(i = seq_len(nrow(Y)), j = seq_len(nrow(Y))) %>%
    filter(i < j) %>%
    transmute(
      body_mass_diff = perm_diff[cbind(i, j)],
      bray_curtis = diet_dist[cbind(i, j)]
    )
  
  cor(perm_pair$body_mass_diff, perm_pair$bray_curtis, method = "spearman")
})

p_perm <- (sum(abs(perm_r) >= abs(obs_r)) + 1) / (length(perm_r) + 1)

support_stats <- tibble(
  spearman_rho = obs_r,
  permutation_p = p_perm
)

support_stats
subtitle_txt <- paste0(
  "Each point is a species pair; rho = ",
  round(obs_r, 2),
  ", permutation p = ",
  signif(p_perm, 2)
)

ggplot(pair_dat, aes(x = body_mass_diff, y = bray_curtis)) +
  geom_point(size = 3, alpha = 0.85, colour = "#2C7FB8") +
  geom_smooth(method = "lm", se = FALSE, linewidth = 0.9, colour = "grey30") +
  labs(
    title = "Diet composition divergence increases with body-size difference",
    subtitle = subtitle_txt,
    x = "Difference in log10(adult body mass)",
    y = "Bray-Curtis diet dissimilarity"
  ) +
  theme_bw(base_size = 13) +
  theme(
    panel.grid.minor = element_blank(),
    plot.title = element_text(face = "bold")
  )

########################
##########################
##############
##########

raw <- read_excel(path_to_file, sheet = sheet_name) %>%
  clean_names()

if ("order" %in% names(raw)) {
  raw <- raw %>% rename(plant_order = order)
}

num_cols <- raw %>%
  select(where(is.numeric)) %>%
  names()

if (length(num_cols) == 0) stop("No numeric sample columns found.")

# -----------------------------
# 5) Drop low-read sample columns
# -----------------------------
col_totals <- raw %>%
  summarise(across(all_of(num_cols), ~ sum(.x, na.rm = TRUE))) %>%
  pivot_longer(everything(), names_to = "sample_col", values_to = "total_reads")

keep_num_cols <- col_totals %>%
  filter(total_reads >= min_reads_per_sample) %>%
  pull(sample_col)

if (length(keep_num_cols) == 0) stop("No samples remain after min_reads_per_sample filtering.")

raw_filt <- raw %>%
  select(plant_order, all_of(keep_num_cols))

# -----------------------------
# 6) Keep only animals with enough retained samples
# -----------------------------
animal_n <- tibble(sample_col = keep_num_cols) %>%
  mutate(animal = normalise_animal_names(sample_col)) %>%
  count(animal, name = "n_samples")

keep_animals <- animal_n %>%
  filter(n_samples >= min_samples_per_animal) %>%
  pull(animal)

if (length(keep_animals) == 0) stop("No animals remain after min_samples_per_animal filtering.")

keep_sample_cols <- tibble(sample_col = keep_num_cols) %>%
  mutate(animal = normalise_animal_names(sample_col)) %>%
  filter(animal %in% keep_animals) %>%
  pull(sample_col)

raw_filt <- raw_filt %>%
  select(plant_order, all_of(keep_sample_cols))

# -----------------------------
# 7) Drop low-abundance plant orders across retained samples
# -----------------------------
order_totals <- raw_filt %>%
  filter(!is.na(plant_order), plant_order != "") %>%
  group_by(plant_order) %>%
  summarise(total_reads = sum(across(all_of(keep_sample_cols)), na.rm = TRUE), .groups = "drop")

keep_orders <- order_totals %>%
  filter(total_reads >= min_reads_per_order) %>%
  pull(plant_order)

raw_filt <- raw_filt %>%
  filter(plant_order %in% keep_orders)

# -----------------------------
# 8) Build sample x order matrix
# -----------------------------
sample_order <- raw_filt %>%
  filter(!is.na(plant_order), plant_order != "") %>%
  pivot_longer(cols = all_of(keep_sample_cols), names_to = "sample_col", values_to = "reads") %>%
  group_by(sample_col, plant_order) %>%
  summarise(reads = sum(reads, na.rm = TRUE), .groups = "drop") %>%
  pivot_wider(names_from = plant_order, values_from = reads, values_fill = 0)

mat_sample <- sample_order %>%
  column_to_rownames("sample_col") %>%
  as.matrix()

# -----------------------------
# 9) Shannon diversity per sample
# -----------------------------
sample_div <- tibble(
  sample_col = rownames(mat_sample),
  shannon = vegan::diversity(mat_sample, index = "shannon")
) %>%
  mutate(
    animal = normalise_animal_names(sample_col)
  ) %>%
  left_join(mass_df, by = "animal") %>%
  filter(!is.na(body_mass_kg))


# -----------------------------
# 10) Linear model
# -----------------------------
fit <- lm(shannon ~ log10(body_mass_kg), data = sample_div)
fit_sum <- summary(fit)

slope_lab <- round(coef(fit)[2], 3)
r2_lab    <- round(fit_sum$r.squared, 3)
p_lab     <- signif(coef(fit_sum)[2, 4], 2)

subtitle_txt <- paste0(
  "LM: slope = ", slope_lab,
  ", R² = ", r2_lab,
  ", p = ", p_lab
)

# -----------------------------
# 11) Plot
# -----------------------------
p <- ggplot(sample_div, aes(x = log10(body_mass_kg), y = shannon, colour = animal)) +
  geom_point(size = 3, alpha = 0.85) +
  geom_smooth(method = "lm", se = TRUE, linewidth = 0.9, colour = "black") +
  labs(
    title = "Plant diversity across anaimal body size",
    subtitle = subtitle_txt,
    x = "log10(adult body mass, kg)",
    y = "Plant diversity per sample (Shannon index)",
    colour = "Animal"
  ) +
  theme_bw(base_size = 13) +
  theme(
    panel.grid.minor = element_blank(),
    plot.title = element_text(face = "bold")
  )

ggsave(
  filename = "plant_diversity_vs_body_size_by_feeding_type.pdf",
  plot = p,
  width = 9,
  height = 7,
  dpi = 300
)


#########################################################################################################
##
##
##

library(readxl)
library(janitor)
library(dplyr)
library(tidyr)
library(tibble)
library(ggplot2)
library(forcats)
library(vegan)

# -----------------------------
# 1) Settings
# -----------------------------
path_to_file <- "Darius_3_formatted_hit_May15_cleanedDS.xlsx"
sheet_name   <- "Darius_3_formatted_hit_May15_KR"

min_reads_per_sample   <- 100
min_samples_per_animal <- 3
min_reads_per_order    <- 100

# -----------------------------
# 2) Helpers
# -----------------------------
base_name <- function(x) gsub("[\\._].*$", "", x)

clean_animal_name <- function(x) {
  dplyr::recode(
    x,
    "bontebok"   = "Bontebok",
    "buffalo"    = "Buffalo",
    "zebra"      = "Burchell's zebra",
    "capezebra"  = "Cape zebra",
    "eland"      = "Eland",
    "elephant"   = "Elephant",
    "giraffe"    = "Giraffe",
    "impala"     = "Impala",
    "rhino"      = "White Rhino",
    "springbok"  = "Springbok",
    "waterbuck"  = "Waterbuck",
    "wildebeest" = "Wildebeest",
    .default = x
  )
}

standardise_animal_name <- function(x) {
  dplyr::recode(
    x,
    "Burchell's zebra" = "Burchell's zebra",
    "White Rhino" = "White rhinoceros",
    .default = x
  )
}

normalise_animal_names <- function(x) {
  x %>%
    base_name() %>%
    clean_animal_name() %>%
    standardise_animal_name()
}

# -----------------------------
# 3) Load data
# -----------------------------
raw <- read_excel(path_to_file, sheet = sheet_name) %>%
  clean_names()

if ("order" %in% names(raw)) {
  raw <- raw %>% rename(plant_order = order)
}

num_cols <- raw %>%
  select(where(is.numeric)) %>%
  names()

if (length(num_cols) == 0) stop("No numeric sample columns found.")

# -----------------------------
# 4) Drop low-read sample columns
# -----------------------------
col_totals <- raw %>%
  summarise(across(all_of(num_cols), ~ sum(.x, na.rm = TRUE))) %>%
  pivot_longer(everything(), names_to = "sample_col", values_to = "total_reads")

keep_num_cols <- col_totals %>%
  filter(total_reads >= min_reads_per_sample) %>%
  pull(sample_col)

if (length(keep_num_cols) == 0) {
  stop("No samples remain after min_reads_per_sample filtering.")
}

raw_filt <- raw %>%
  select(plant_order, all_of(keep_num_cols))

# -----------------------------
# 5) Keep only animals with enough retained samples
# -----------------------------
animal_n <- tibble(sample_col = keep_num_cols) %>%
  mutate(animal = normalise_animal_names(sample_col)) %>%
  count(animal, name = "n_samples")

keep_animals <- animal_n %>%
  filter(n_samples >= min_samples_per_animal) %>%
  pull(animal)

if (length(keep_animals) == 0) {
  stop("No animals remain after min_samples_per_animal filtering.")
}

keep_sample_cols <- tibble(sample_col = keep_num_cols) %>%
  mutate(animal = normalise_animal_names(sample_col)) %>%
  filter(animal %in% keep_animals) %>%
  pull(sample_col)

raw_filt <- raw_filt %>%
  select(plant_order, all_of(keep_sample_cols))

# -----------------------------
# 6) Drop low-abundance plant orders across retained samples
# -----------------------------
order_totals <- raw_filt %>%
  filter(!is.na(plant_order), plant_order != "") %>%
  group_by(plant_order) %>%
  summarise(
    total_reads = sum(across(all_of(keep_sample_cols)), na.rm = TRUE),
    .groups = "drop"
  )

keep_orders <- order_totals %>%
  filter(total_reads >= min_reads_per_order) %>%
  pull(plant_order)

raw_filt <- raw_filt %>%
  filter(plant_order %in% keep_orders)

# -----------------------------
# 7) Build sample x order matrix
# -----------------------------
sample_order <- raw_filt %>%
  filter(!is.na(plant_order), plant_order != "") %>%
  pivot_longer(cols = all_of(keep_sample_cols), names_to = "sample_col", values_to = "reads") %>%
  group_by(sample_col, plant_order) %>%
  summarise(reads = sum(reads, na.rm = TRUE), .groups = "drop") %>%
  pivot_wider(names_from = plant_order, values_from = reads, values_fill = 0)

mat_sample <- sample_order %>%
  column_to_rownames("sample_col") %>%
  as.matrix()

# -----------------------------
# 8) Calculate diversity metrics per sample
# -----------------------------
sample_div <- tibble(
  sample_col = rownames(mat_sample),
  animal = normalise_animal_names(rownames(mat_sample)),
  shannon = vegan::diversity(mat_sample, index = "shannon"),
  simpson = vegan::diversity(mat_sample, index = "simpson"),
  invsimpson = vegan::diversity(mat_sample, index = "invsimpson"),
  richness = rowSums(mat_sample > 0)
)

# optional: order animals by median Shannon
animal_levels <- sample_div %>%
  group_by(animal) %>%
  summarise(med_shannon = median(shannon, na.rm = TRUE), .groups = "drop") %>%
  arrange(desc(med_shannon)) %>%
  pull(animal)

sample_div <- sample_div %>%
  mutate(animal = factor(animal, levels = animal_levels))

# -----------------------------
# 9) Long format for plotting all metrics
# -----------------------------
div_long <- sample_div %>%
  pivot_longer(
    cols = c(shannon, simpson, invsimpson, richness),
    names_to = "metric",
    values_to = "value"
  ) %>%
  mutate(
    metric = factor(
      metric,
      levels = c("shannon", "simpson", "invsimpson", "richness"),
      labels = c("Shannon", "Simpson", "Inverse Simpson", "Richness")
    )
  )

# -----------------------------
# 10) Plot: all diversity metrics by animal
# -----------------------------
p_div <- ggplot(div_long, aes(x = animal, y = value)) +
  geom_boxplot(outlier.shape = NA, width = 0.7) +
  geom_jitter(width = 0.15, alpha = 0.6, size = 1.8) +
  facet_wrap(~ metric, scales = "free_y", ncol = 2) +
  labs(
    x = "Animal",
    y = "Diversity metric value",
    title = "Plant dietary diversity metrics across herbivore species"
  ) +
  theme_bw(base_size = 13) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.grid.minor = element_blank(),
    plot.title = element_text(face = "bold")
  )

p_div

# -----------------------------
# 11) Save plot
# -----------------------------
ggsave(
  filename = "diet_diversity_metrics_by_animal.png",
  plot = p_div,
  width = 12,
  height = 8,
  dpi = 300
)

animal_cols <- c(
  "Bontebok" = "#1b9e77",
  "Buffalo" = "#d95f02",
  "Burchell's zebra" = "#7570b3",
  "Cape zebra" = "#e7298a",
  "Eland" = "#66a61e",
  "Elephant" = "#e6ab02",
  "Giraffe" = "#a6761d",
  "Impala" = "#1f78b4",
  "Plains zebra" = "#6a3d9a",
  "Springbok" = "#b15928",
  "Waterbuck" = "#fb9a99",
  "White rhinoceros" = "#b2df8a",
  "Wildebeest" = "#cab2d6"
)

p_div_nice <- ggplot(div_long, aes(x = animal, y = value)) +
  geom_boxplot(
    aes(fill = animal),
    width = 0.65,
    outlier.shape = NA,
    colour = "black",
    linewidth = 0.5,
    alpha = 0.85
  ) +
  geom_jitter(
    width = 0.2,
    alpha = 0.45,
    size = 1.4,
    colour = "grey45"
  ) +
  facet_wrap(~ metric, scales = "free_y", ncol = 2) +
  scale_fill_manual(values = animal_cols) +
  labs(
    x = NULL,
    y = "Value",
    title = "Dietary diversity across anaimals"
  ) +
  theme_bw(base_size = 13) +
  theme(
    strip.background = element_rect(fill = "grey92", colour = "black"),
    strip.text = element_text(face = "bold"),
    panel.grid.major.x = element_blank(),
    panel.grid.minor = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1),
    plot.title = element_text(face = "bold"),
    legend.position = "none"
  )

p_div_nice

install.packages("patchwork")
library(patchwork)

# assuming:
# p = animal size plot
# p_div_nice = diversity metrics plot

combined_plot <- p_div_nice + p +
  plot_layout(widths = c(1.4, 1)) +
  plot_annotation(tag_levels = "A")

combined_plot
ggsave(
  filename = "combined_diversity_and_body_size_plots.pdf",
  plot = combined_plot,
  width = 16,
  height = 8,
  dpi = 300
)

# -----------------------------
# 12) Optional: export summary stats
# -----------------------------
div_summary <- sample_div %>%
  group_by(animal) %>%
  summarise(
    n_samples = n(),
    shannon_median = median(shannon, na.rm = TRUE),
    shannon_iqr = IQR(shannon, na.rm = TRUE),
    simpson_median = median(simpson, na.rm = TRUE),
    invsimpson_median = median(invsimpson, na.rm = TRUE),
    richness_median = median(richness, na.rm = TRUE),
    .groups = "drop"
  )

write.csv(div_summary, "diet_diversity_metrics_by_animal_summary.csv", row.names = FALSE)