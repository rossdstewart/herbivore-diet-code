# ================================
# Compare Darius eDNA vs RSA diet matrix
# FULL FROM SCRATCH
# ================================

# -----------------------------
# 0) Packages
# -----------------------------
library(readxl)
library(janitor)
library(dplyr)
library(tidyr)
library(tibble)
library(stringr)
library(ggplot2)
library(vegan)
library(writexl)
library(forcats)

# -----------------------------
# 1) User settings
# -----------------------------

setwd("C:/Users/Admin/Desktop/R/Rose et al")

path_darius  <- "Darius_3_formatted_hit_May15_cleanedDS.xlsx"
sheet_darius <- "Darius_3_formatted_hit_May15_KR"

path_rsa  <- "RSA_combined_RRA_order_family.xlsx"
sheet_rsa <- "Order_RRA"

min_reads_per_sample   <- 100   # drop raw sample columns with <100 total reads
min_samples_per_animal <- 3     # drop animals represented by <3 retained samples
min_reads_per_order    <- 100   # drop plant orders with <100 total reads across all animals

# -----------------------------
# 2) Helper functions
# -----------------------------
base_name <- function(x) {
  gsub("[\\._].*$", "", x)
}

clean_animal_name <- function(x) {
  x <- tolower(trimws(x))
  
  dplyr::recode(
    x,
    "bontebok"      = "Bontebok",
    "buffalo"       = "Buffalo",
    "zebra"         = "Burchell's zebra",
    "burchellszebra"= "Burchell's zebra",
    "burchell's zebra" = "Burchell's zebra",
    "capezebra"     = "Cape zebra",
    "cape zebra"    = "Cape zebra",
    "eland"         = "Eland",
    "elephant"      = "Elephant",
    "giraffe"       = "Giraffe",
    "impala"        = "Impala",
    "White"         = "White Rhino",
    "rhino"         = "White Rhino",
    "white_rhino"   = "White Rhino",
    "white_rhinoceros"   = "White Rhino",
    "springbok"     = "Springbok",
    "waterbuck"     = "Waterbuck",
    "wildebeest"    = "Wildebeest",
    .default = str_to_title(x)
  )
}



# final standard that BOTH datasets must use
standardise_animal_name <- function(x) {
  x <- trimws(x)
  
  dplyr::recode(
    x,
    "Burchell's zebra"   = "Plains zebra",
    "White Rhino"        = "White rhinoceros",
    "Rhino"              = "White rhinoceros",
    "White"              = "White rhinoceros",
    "White rhino"        = "White rhinoceros",
    "Plains Zebra"       = "Plains zebra",
    "Cape Zebra"         = "Cape zebra",
    .default = x
  )
}

# applies the full pipeline to vectors of names
normalise_animal_names <- function(x) {
  x %>%
    base_name() %>%
    clean_animal_name() %>%
    standardise_animal_name()
}

# aggregate duplicate columns after renaming
aggregate_duplicate_columns <- function(mat) {
  mat <- as.matrix(mat)
  
  uniq <- unique(colnames(mat))
  
  out <- sapply(uniq, function(nm) {
    cols <- which(colnames(mat) == nm)
    if (length(cols) == 1) {
      mat[, cols]
    } else {
      rowSums(mat[, cols, drop = FALSE], na.rm = TRUE)
    }
  })
  
  out <- as.matrix(out)
  
  if (is.null(colnames(out))) {
    colnames(out) <- uniq
  }
  
  if (is.null(rownames(out))) {
    rownames(out) <- rownames(mat)
  }
  
  out
}

# -----------------------------
# 3) Load Darius raw data
# -----------------------------
raw <- read_excel(path_darius, sheet = sheet_darius) %>%
  clean_names()

if ("order" %in% names(raw)) {
  raw <- raw %>% rename(plant_order = order)
}

num_cols <- raw %>%
  select(where(is.numeric)) %>%
  names()

if (length(num_cols) == 0) stop("No numeric diet columns found in Darius file.")

# -----------------------------
# 4) Drop low-read Darius sample columns
# -----------------------------
col_totals <- raw %>%
  summarise(across(all_of(num_cols), ~ sum(.x, na.rm = TRUE))) %>%
  pivot_longer(everything(), names_to = "sample_col", values_to = "total_reads")

keep_num_cols <- col_totals %>%
  filter(total_reads >= min_reads_per_sample) %>%
  pull(sample_col)

raw_filt <- raw %>%
  select(any_of(c("plant_order", "family", "species")), all_of(keep_num_cols))

num_cols <- keep_num_cols

if (length(num_cols) == 0) {
  stop("No sample columns remain after filtering by min_reads_per_sample.")
}

# -----------------------------
# 5) Count retained Darius samples per animal
# -----------------------------
animal_n <- tibble(col = num_cols) %>%
  mutate(animal = normalise_animal_names(col)) %>%
  count(animal, name = "n_samples")

keep_animals <- animal_n %>%
  filter(n_samples >= min_samples_per_animal) %>%
  pull(animal)

if (length(keep_animals) == 0) {
  stop("No animals remain after filtering by min_samples_per_animal.")
}

# -----------------------------
# 6) Aggregate Darius to ORDER x ANIMAL
# -----------------------------
ord_animal <- raw_filt %>%
  select(plant_order, all_of(num_cols)) %>%
  filter(!is.na(plant_order), plant_order != "") %>%
  pivot_longer(cols = all_of(num_cols), names_to = "col", values_to = "value") %>%
  mutate(animal = normalise_animal_names(col)) %>%
  filter(animal %in% keep_animals) %>%
  group_by(plant_order, animal) %>%
  summarise(value = sum(value, na.rm = TRUE), .groups = "drop") %>%
  pivot_wider(names_from = animal, values_from = value, values_fill = 0) %>%
  arrange(plant_order)

# -----------------------------
# 7) Drop low-abundance plant orders in Darius
# -----------------------------
ord_animal <- ord_animal %>%
  mutate(total_reads = rowSums(across(-plant_order), na.rm = TRUE)) %>%
  filter(total_reads >= min_reads_per_order) %>%
  select(-total_reads)

if (nrow(ord_animal) < 2) {
  stop("After filtering, <2 orders remain in Darius. Lower min_reads_per_order.")
}

# -----------------------------
# 8) Build Darius matrix and convert to relative composition
# -----------------------------
darius_mat <- ord_animal %>%
  column_to_rownames("plant_order") %>%
  as.matrix()

darius_rel <- sweep(darius_mat, 2, colSums(darius_mat), "/")
darius_rel[is.na(darius_rel)] <- 0

# -----------------------------
# 9) Load RSA order-level matrix
# -----------------------------
rsa_raw <- read_excel(path_rsa, sheet = sheet_rsa) %>%
  clean_names()

if ("order" %in% names(rsa_raw)) {
  rsa_raw <- rsa_raw %>% rename(plant_order = order)
}

if (!"plant_order" %in% names(rsa_raw)) {
  stop("RSA sheet must contain a plant_order or order column.")
}

# keep order + animal columns
rsa_df <- rsa_raw %>%
  filter(!is.na(plant_order), plant_order != "")

# identify RSA animal columns
rsa_animal_cols <- setdiff(names(rsa_df), "plant_order")

if (length(rsa_animal_cols) == 0) {
  stop("No animal columns found in RSA file.")
}

# force numeric
rsa_df <- rsa_df %>%
  mutate(across(all_of(rsa_animal_cols), as.numeric))

# -----------------------------
# 10) Clean RSA animal names THE SAME WAY
# -----------------------------
rsa_clean_names <- normalise_animal_names(rsa_animal_cols)

names(rsa_df)[match(rsa_animal_cols, names(rsa_df))] <- rsa_clean_names

# build matrix
rsa_mat <- rsa_df %>%
  column_to_rownames("plant_order") %>%
  as.matrix()

# aggregate duplicate columns created by name cleaning
rsa_mat <- aggregate_duplicate_columns(rsa_mat)

# convert RSA to relative composition just in case
rsa_rel <- sweep(rsa_mat, 2, colSums(rsa_mat), "/")
rsa_rel[is.na(rsa_rel)] <- 0

# -----------------------------
# 11) Restrict to shared animals
# -----------------------------
shared_animals <- intersect(colnames(darius_rel), colnames(rsa_rel))

if (length(shared_animals) < 2) {
  stop(
    paste0(
      "Too few shared animals after cleaning.\nDarius animals: ",
      paste(colnames(darius_rel), collapse = ", "),
      "\nRSA animals: ",
      paste(colnames(rsa_rel), collapse = ", ")
    )
  )
}

darius_rel2 <- darius_rel[, shared_animals, drop = FALSE]
rsa_rel2    <- rsa_rel[, shared_animals, drop = FALSE]

# -----------------------------
# 12) Restrict to shared plant orders
# -----------------------------
shared_orders <- intersect(rownames(darius_rel2), rownames(rsa_rel2))

if (length(shared_orders) < 2) {
  stop("Too few shared plant orders between datasets.")
}

darius_shared <- darius_rel2[shared_orders, , drop = FALSE]
rsa_shared    <- rsa_rel2[shared_orders, , drop = FALSE]

# renormalise after intersecting shared orders
darius_shared <- sweep(darius_shared, 2, colSums(darius_shared), "/")
rsa_shared    <- sweep(rsa_shared, 2, colSums(rsa_shared), "/")

darius_shared[is.na(darius_shared)] <- 0
rsa_shared[is.na(rsa_shared)] <- 0

# -----------------------------
# 13) Bray-Curtis by animal
# -----------------------------
bray_by_animal <- tibble(
  animal = shared_animals,
  bray_curtis = sapply(shared_animals, function(a) {
    as.numeric(
      vegan::vegdist(
        rbind(darius_shared[, a], rsa_shared[, a]),
        method = "bray"
      )
    )
  })
) %>%
  arrange(desc(bray_curtis))

print(bray_by_animal)

# -----------------------------
# 14) Long comparison table
# -----------------------------
comp_long <- bind_rows(
  as.data.frame(darius_shared) %>%
    rownames_to_column("plant_order") %>%
    pivot_longer(-plant_order, names_to = "animal", values_to = "rra") %>%
    mutate(dataset = "Darius"),
  as.data.frame(rsa_shared) %>%
    rownames_to_column("plant_order") %>%
    pivot_longer(-plant_order, names_to = "animal", values_to = "rra") %>%
    mutate(dataset = "RSA")
) %>%
  pivot_wider(names_from = dataset, values_from = rra, values_fill = 0) %>%
  mutate(
    delta = Darius - RSA,
    abs_delta = abs(delta)
  )

# -----------------------------
# 15) Top shifts per animal
# -----------------------------
top_shifts <- comp_long %>%
  group_by(animal) %>%
  slice_max(abs_delta, n = 10, with_ties = FALSE) %>%
  arrange(animal, desc(abs_delta))

# -----------------------------
# 16) Plot: Bray-Curtis by animal
# -----------------------------
p_bray <- ggplot(bray_by_animal,
                 aes(x = reorder(animal, bray_curtis), y = bray_curtis)) +
  geom_col() +
  coord_flip() +
  labs(
    x = NULL,
    y = "Bray-Curtis dissimilarity",
    title = "Diet shift between Darius eDNA and RSA"
  ) +
  theme_bw(base_size = 12)

print(p_bray)

# -----------------------------
# 17) Plot: Shift heatmap
# -----------------------------
heat_dat <- comp_long %>%
  group_by(animal) %>%
  slice_max(abs_delta, n = 8, with_ties = FALSE) %>%
  ungroup()

p_heat <- ggplot(
  heat_dat,
  aes(x = animal,
      y = fct_reorder(plant_order, abs_delta, .fun = max),
      fill = delta)
) +
  geom_tile() +
  labs(
    x = NULL,
    y = "Plant order",
    fill = "Darius - RSA",
    title = "Largest order-level shifts"
  ) +
  theme_bw(base_size = 12) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

print(p_heat)

# -----------------------------
# 18) Optional composition plot
# -----------------------------
top_orders_global <- comp_long %>%
  group_by(plant_order) %>%
  summarise(total_shift = sum(abs_delta, na.rm = TRUE), .groups = "drop") %>%
  slice_max(total_shift, n = 12)

stack_dat <- comp_long %>%
  filter(plant_order %in% top_orders_global$plant_order) %>%
  select(animal, plant_order, Darius, RSA) %>%
  pivot_longer(cols = c(Darius, RSA), names_to = "dataset", values_to = "rra")

p_stack <- ggplot(stack_dat, aes(x = dataset, y = rra, fill = plant_order)) +
  geom_col(position = "fill") +
  facet_wrap(~ animal) +
  labs(
    x = NULL,
    y = "Relative composition",
    fill = "Plant order",
    title = "Diet composition by animal"
  ) +
  theme_bw(base_size = 11)

print(p_stack)

# -----------------------------
# 19) Save outputs
# -----------------------------
ggsave("bray_by_animal.png", p_bray, width = 8, height = 5, dpi = 300)
ggsave("delta_heatmap.png", p_heat, width = 9, height = 7, dpi = 300)
ggsave("diet_composition_facets.png", p_stack, width = 12, height = 8, dpi = 300)

write_xlsx(
  list(
    retained_samples_per_animal = animal_n,
    darius_relative = as.data.frame(darius_shared) %>% rownames_to_column("plant_order"),
    rsa_relative = as.data.frame(rsa_shared) %>% rownames_to_column("plant_order"),
    bray_by_animal = bray_by_animal,
    comparison_long = comp_long,
    top_shifts = top_shifts
  ),
  "edna_vs_rsa_comparison.xlsx"
)

# -----------------------------
# 20) Console checks
# -----------------------------
cat("\nDarius animals after cleaning:\n")
print(colnames(darius_rel))

cat("\nRSA animals after cleaning:\n")
print(colnames(rsa_rel))

cat("\nShared animals used:\n")
print(shared_animals)

cat("\nNumber of shared plant orders:\n")
print(length(shared_orders))


###################################################################################
# =============================
# 15) REAL PAIRED ANALYSIS
# =============================

library(vegan)
library(ggplot2)
library(dplyr)
library(tidyr)
library(tibble)
library(forcats)

# --------------------------------
# Helper: CLR transform
# --------------------------------
clr_transform <- function(mat, pseudocount = 1e-6) {
  mat <- as.matrix(mat)
  mat <- mat + pseudocount
  gm <- exp(rowMeans(log(mat)))
  log(mat / gm)
}

# --------------------------------
# 15A) Build sample x order matrix
# Each row = one animal-community sample
# --------------------------------
sample_df <- bind_rows(
  as.data.frame(t(darius_shared)) %>%
    rownames_to_column("animal") %>%
    mutate(community = "Darius"),
  as.data.frame(t(rsa_shared)) %>%
    rownames_to_column("animal") %>%
    mutate(community = "RSA")
) %>%
  relocate(animal, community)

meta <- sample_df %>%
  select(animal, community)

X_bray <- sample_df %>%
  select(-animal, -community) %>%
  as.matrix()

rownames(X_bray) <- paste(meta$animal, meta$community, sep = "_")

# --------------------------------
# 15B) Bray-Curtis distance
# --------------------------------
dist_bray <- vegdist(X_bray, method = "bray")

# --------------------------------
# 15C) Paired PERMANOVA
# community effect with animal pairing respected
# --------------------------------
set.seed(123)

adon_bray <- adonis2(
  dist_bray ~ community,
  data = meta,
  permutations = 9999,
  strata = meta$animal
)

print(adon_bray)

# --------------------------------
# 15D) Aitchison / CLR check
# Euclidean distance on CLR-transformed compositions
# --------------------------------
X_clr <- clr_transform(X_bray)
dist_aitchison <- dist(X_clr, method = "euclidean")

set.seed(123)

adon_aitchison <- adonis2(
  dist_aitchison ~ community,
  data = meta,
  permutations = 9999,
  strata = meta$animal
)

print(adon_aitchison)

# --------------------------------
# 15E) Ordination plot with paired lines
# --------------------------------
pcoa <- cmdscale(dist_bray, k = 2, eig = TRUE)

ord_df <- as.data.frame(pcoa$points)
colnames(ord_df) <- c("Axis1", "Axis2")
ord_df <- bind_cols(meta, ord_df)

p_ord <- ggplot(ord_df, aes(Axis1, Axis2, colour = community)) +
  geom_line(aes(group = animal), linewidth = 0.5, alpha = 0.7, colour = "grey50") +
  geom_point(size = 3) +
  geom_text(aes(label = animal), hjust = -0.1, vjust = 0.3, size = 3) +
  labs(
    title = "PCoA of diet composition",
    subtitle = "Lines connect the same animal across communities",
    x = "PCoA 1",
    y = "PCoA 2"
  ) +
  theme_bw(base_size = 12)

print(p_ord)

ggsave("ordination_paired_bray.png", p_ord, width = 8, height = 6, dpi = 300)

# --------------------------------
# SUPPORTING TEST FOR THE ORDINATION
# Are animals more similar to themselves across communities
# than expected by chance?
# --------------------------------

# cross-community Bray-Curtis matrix:
# rows = Darius animals, cols = RSA animals
cross_mat <- outer(
  shared_animals,
  shared_animals,
  Vectorize(function(a1, a2) {
    as.numeric(vegdist(rbind(darius_shared[, a1], rsa_shared[, a2]), method = "bray"))
  })
)

rownames(cross_mat) <- shared_animals
colnames(cross_mat) <- shared_animals

# observed mean distance for correct animal matching
obs_mean_same <- mean(diag(cross_mat))

# permutation test
set.seed(123)
nperm <- 9999

perm_mean_same <- replicate(nperm, {
  perm_cols <- sample(shared_animals)
  mean(cross_mat[cbind(shared_animals, perm_cols)])
})

p_value <- (sum(perm_mean_same <= obs_mean_same) + 1) / (nperm + 1)

matching_test <- tibble(
  observed_mean_same_animal_bray = obs_mean_same,
  expected_mean_random_matching = mean(perm_mean_same),
  p_value = p_value
)

print(matching_test)
mean(diag(cross_mat))
mean(cross_mat[row(cross_mat) != col(cross_mat)])

#
within_mean  <- mean(diag(cross_mat))
between_mean <- mean(cross_mat[row(cross_mat) != col(cross_mat)])

summary_dist <- tibble(
  within_animal_mean = within_mean,
  between_animal_mean = between_mean,
  difference = between_mean - within_mean
)

print(summary_dist)
#######################
# -----------------------------
# Build sample x order matrix
# -----------------------------
sample_df <- bind_rows(
  as.data.frame(t(darius_shared)) %>%
    rownames_to_column("animal") %>%
    mutate(community = "Darius"),
  as.data.frame(t(rsa_shared)) %>%
    rownames_to_column("animal") %>%
    mutate(community = "RSA")
)

meta <- sample_df %>%
  select(animal, community)

X <- sample_df %>%
  select(-animal, -community) %>%
  as.matrix()

rownames(X) <- paste(meta$animal, meta$community, sep = "_")

# -----------------------------
# Bray-Curtis distance
# -----------------------------
dist_bray <- vegdist(X, method = "bray")

# -----------------------------
# Paired PERMANOVA
# -----------------------------
set.seed(123)

adon_bray <- adonis2(
  dist_bray ~ community,
  data = meta,
  permutations = 9999,
  strata = meta$animal
)

print(adon_bray)

bd <- betadisper(dist_bray, group = meta$community)
anova(bd)
permutest(bd, permutations = 9999)
##########
#
########
##



library(vegan)
library(ggplot2)
library(dplyr)
library(tibble)
library(ggrepel)

# --------------------------------
# 1) Build sample x order matrix
# --------------------------------
sample_df <- bind_rows(
  as.data.frame(t(darius_shared)) %>%
    rownames_to_column("animal") %>%
    mutate(community = "Fynbos"),
  as.data.frame(t(rsa_shared)) %>%
    rownames_to_column("animal") %>%
    mutate(community = "Savanna")
)

meta <- sample_df %>%
  select(animal, community)

X <- sample_df %>%
  select(-animal, -community) %>%
  as.matrix()

rownames(X) <- paste(meta$animal, meta$community, sep = "_")
  
# --------------------------------
# 2) Bray-Curtis + PCoA
# --------------------------------
dist_bray <- vegdist(X, method = "bray")
pcoa <- cmdscale(dist_bray, k = 2, eig = TRUE)

ord_df <- as.data.frame(pcoa$points)
colnames(ord_df) <- c("PCoA1", "PCoA2")
ord_df <- bind_cols(meta, ord_df)

# axis variance
eig_pos <- pcoa$eig[pcoa$eig > 0]
axis_var <- round(100 * eig_pos / sum(eig_pos), 1)

# --------------------------------
# 3) Feeding guilds
# --------------------------------
guild_df <- tibble(
  animal = c(
    "Buffalo", "Elephant", "Giraffe", "Impala",
    "Plains zebra", "Waterbuck", "White rhinoceros", "Wildebeest"
  ),
  feeding_guild = c(
    "Grazer", "Mixed feeder", "Browser", "Mixed feeder",
    "Grazer", "Grazer", "Grazer", "Grazer"
  )
)

ord_df <- ord_df %>%
  left_join(guild_df, by = "animal")

# --------------------------------
# 4) Label positions
# one label per animal at midpoint of the two communities
# --------------------------------
label_df <- ord_df %>%
  group_by(animal, feeding_guild) %>%
  summarise(
    PCoA1 = mean(PCoA1),
    PCoA2 = mean(PCoA2),
    .groups = "drop"
  )

# --------------------------------
# 5) Manual colours and shapes
# --------------------------------
community_cols <- c(
  "Fynbos" = "#1b9e77",
  "Savanna" = "#d95f02"
)

guild_shapes <- c(
  "Browser" = 17,
  "Grazer" = 16,
  "Mixed feeder" = 15
)

# --------------------------------
# 6) Plot
# --------------------------------
p_ord <- ggplot(ord_df, aes(x = PCoA1, y = PCoA2)) +
  geom_line(
    aes(group = animal),
    colour = "grey40",
    linewidth = 0.8,
    alpha = 0.8
  ) +
  geom_point(
    aes(colour = community, shape = feeding_guild),
    size = 4
  ) +
  ggrepel::geom_text_repel(
    data = label_df,
    aes(label = animal),
    colour = "black",
    size = 3.5,
    show.legend = FALSE,
    max.overlaps = Inf
  ) +
  scale_colour_manual(values = community_cols) +
  scale_shape_manual(values = guild_shapes) +
  labs(
    title = "Diet composition shifts between communities",
    subtitle = "PCoA based on Bray-Curtis dissimilarities of order-level diet composition",
    x = paste0("PCoA 1 (", axis_var[1], "%)"),
    y = paste0("PCoA 2 (", axis_var[2], "%)"),
    colour = "Community",
    shape = "Feeding guild",
    caption = "Matched animal diets differed between communities (paired PERMANOVA: F = 2.13, R² = 0.151, p = 0.031), with no difference in dispersion (p = 0.905)."
  ) +
  theme_bw(base_size = 12) +
  theme(
    panel.grid = element_blank(),
    plot.title = element_text(face = "bold"),
    axis.title = element_text(face = "bold"),
    legend.position = "right"
  )

print(p_ord)

ggsave(
  "pcoa_diet_shift_fynbos_savanna.pdf",
  p_ord,
  width = 10,
  height = 7,
  dpi = 300
)


library(dplyr)
library(tidyr)
library(tibble)

animal_drivers <- bind_rows(lapply(colnames(darius_shared), function(a) {
  tibble(
    animal = a,
    plant_order = rownames(darius_shared),
    fynbos = darius_shared[, a],
    savanna = rsa_shared[, a],
    delta = darius_shared[, a] - rsa_shared[, a]
  )
})) %>%
  group_by(animal) %>%
  slice_max(abs(delta), n = 5, with_ties = FALSE) %>%
  arrange(animal, desc(abs(delta)))

animal_drivers %>%
  arrange(desc(abs(delta)))

animal_delta <- bind_rows(lapply(colnames(darius_shared), function(a) {
  tibble(
    animal = a,
    plant_order = rownames(darius_shared),
    fynbos = darius_shared[, a],
    savanna = rsa_shared[, a],
    delta = darius_shared[, a] - rsa_shared[, a]
  )
})) %>%
  arrange(animal, desc(abs(delta)))

animal_delta

unique(animal_delta$animal)
animal_delta %>%
  filter(animal == "Impala") %>%
  arrange(desc(abs(delta)))
##
axis_drivers <- tibble(
  plant_order = colnames(X),
  PCoA_1 = apply(X, 2, function(v) cor(v, pcoa$points[, 1], method = "spearman")),
  PCoA_2 = apply(X, 2, function(v) cor(v, pcoa$points[, 2], method = "spearman"))
)

# strongest for Axis 1
PCoA_1_top_pos <- axis_drivers %>%
  arrange(desc(PCoA_1)) %>%
  slice(1:5)

PCoA_1_top_neg <- axis_drivers %>%
  arrange(PCoA_1) %>%
  slice(1:5)

# strongest for Axis 2
PCoA_2_top_pos <- axis_drivers %>%
  arrange(desc(PCoA_2)) %>%
  slice(1:5)

PCoA_2_top_neg <- axis_drivers %>%
  arrange(PCoA_2) %>%
  slice(1:5)

PCoA_1_top_pos # positive x axis
PCoA_1_top_neg # negative x axis 
PCoA_2_top_pos # positive y axis
PCoA_2_top_neg # negative y axis
