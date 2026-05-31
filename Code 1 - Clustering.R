# ---- Packages ----
library(readxl)
library(tibble)
library(dplyr)
library(tidyr)
library(janitor)
library(vegan)
library(cluster)   # for silhouette
# optional for a nicer heatmap:
# install.packages("pheatmap"); library(pheatmap)

setwd("C:/Users/Admin/Desktop/R/Rose et al")


# ---- Load and prepare (same as before) ----
path_to_file <- "Darius_3_formatted_hit_May15_cleanedDS.xlsx"
sheet_name   <- "Darius_3_formatted_hit_May15_KR"

raw <- read_excel(path_to_file, sheet = sheet_name) %>% clean_names() 

num_cols <- raw %>% select(where(is.numeric)) %>% colnames()

# total reads per numeric column
col_totals <- raw %>%
  summarise(across(all_of(num_cols), ~ sum(.x, na.rm = TRUE))) %>%
  pivot_longer(everything(), names_to = "sample_col", values_to = "total_reads")

# keep only sample columns with >= 100 reads
keep_num_cols <- col_totals %>%
  filter(total_reads >= 100) %>%
  pull(sample_col)

# drop low-read numeric columns, keep metadata columns
raw <- raw %>%
  select(-all_of(setdiff(num_cols, keep_num_cols)))

# update num_cols after filtering
num_cols <- keep_num_cols

colSums(raw %>% select(where(is.numeric)), na.rm = TRUE)

########

base_name <- function(x) gsub("[\\._].*$", "", x)  # strip after first . or _
diet <- raw %>% select(family, all_of(num_cols))

# how many sample columns belong to each animal
# count how many sample columns belong to each animal
animal_n <- tibble(col = num_cols) %>%
  mutate(animal = base_name(col)) %>%
  count(animal, name = "n_samples")

# only keep animals with at least 3 samples
keep_animals <- animal_n %>%
  filter(n_samples >= 3) %>%
  pull(animal)

fam_animal <- diet %>%
  pivot_longer(cols = all_of(num_cols), names_to = "col", values_to = "value") %>%
  mutate(animal = base_name(col)) %>%
  filter(animal %in% keep_animals) %>%
  group_by(family, animal) %>%
  summarise(value = sum(value, na.rm = TRUE), .groups = "drop") %>%
  pivot_wider(names_from = animal, values_from = value, values_fill = 0) %>%
  arrange(family) %>%
  filter(rowSums(across(-family)) > 0)

# rename columns to include sample counts
new_names <- animal_n %>%
  filter(animal %in% keep_animals) %>%
  mutate(new_name = paste0(animal, " (n=", n_samples, ")")) %>%
  select(animal, new_name) %>%
  tibble::deframe()

fam_animal <- fam_animal %>%
  rename_with(~ ifelse(.x %in% names(new_names), new_names[.x], .x))

mat <- fam_animal %>%
  column_to_rownames("family") %>%
  as.matrix()

# ---- Bray–Curtis on relative composition ----
mat_rel <- sweep(mat, 2, colSums(mat), FUN = "/")
mat_rel[is.na(mat_rel)] <- 0
dist_animals <- vegdist(t(mat_rel), method = "bray")  # animals x animals

# ---- Hierarchical clustering ----
hc <- hclust(dist_animals, method = "average")  # UPGMA

# ---- Choose number of clusters (k) via silhouette ----
# Try k = 2..8 and pick the maximum average silhouette width
ks <- 2:min(8, ncol(mat_rel)-1)
sil_mean <- sapply(ks, function(k) {
  cl <- cutree(hc, k = k)
  mean(silhouette(cl, dist_animals)[, "sil_width"])
})
k_opt <- ks[which.max(sil_mean)]
k_opt

# ---- Final assignments ----
groups <- cutree(hc, k = k_opt)

# Tidy table of clusters
cluster_table <- data.frame(
  animal = names(groups),
  cluster = paste0("C", groups),
  row.names = NULL
) %>% arrange(cluster, animal)
print(cluster_table)

# Optional: save
# write.csv(cluster_table, "animal_clusters.csv", row.names = FALSE)

# ---- Plot dendrogram with cluster rectangles ----
plot(hc, main = sprintf("Animal clustering by plant family (Bray–Curtis)  |  k = %d", k_opt),
     xlab = "", sub = "")
rect.hclust(hc, k = k_opt, border = "red")

# ---- Cluster-ordered heatmap (optional) ----
# Reorder columns (animals) by dendrogram order and cluster
ord <- hc$order
mat_rel_ord <- mat_rel[, ord, drop = FALSE]

# Base R heatmap (no scaling to keep compositions interpretable):
# heatmap requires a distance-based clustering unless we disable. We’ll just image + labels:
par(mar = c(8, 4, 2, 10))
image(t(mat_rel_ord[nrow(mat_rel_ord):1, ]), axes = FALSE, main = "Family-level diet (relative)")
axis(1, at = seq(0, 1, length.out = ncol(mat_rel_ord)),
     labels = colnames(mat_rel_ord), las = 2, cex.axis = 0.7)
axis(2, at = seq(0, 1, length.out = nrow(mat_rel_ord)),
     labels = rev(rownames(mat_rel_ord)), las = 2, cex.axis = 0.5)
box()
# If you prefer a nice heatmap: uncomment pheatmap lines at the top and use:
# pheatmap::pheatmap(mat_rel_ord, cluster_rows = TRUE, cluster_cols = FALSE,
#                    main = sprintf("Family-level diet (columns ordered by hc, k=%d)", k_opt))

# ---- Quick sanity checks you might find useful ----
# How similar are animals within each cluster?
by(groups, groups, function(g) NULL)  # prints cluster IDs
aggregate(as.matrix(dist_animals), list(cluster = groups[col(dist_animals)[,1]]), mean)

#######################################################################################################
######## ORDER

# ---- Packages ----
library(readxl)
library(dplyr)
library(tidyr)
library(janitor)
library(vegan)
library(cluster)

# ---- Bray–Curtis on relative composition ----
mat_rel <- sweep(mat, 2, colSums(mat), FUN = "/")
mat_rel[is.na(mat_rel)] <- 0
dist_animals <- vegdist(t(mat_rel), method = "bray")  # animals x animals

# ---- Hierarchical clustering (UPGMA) ----
hc <- hclust(dist_animals, method = "average")

# ---- Pick k via silhouette (2..8 or up to n-1) ----
ks <- 2:min(8, ncol(mat_rel) - 1)
sil_mean <- sapply(ks, function(k) {
  cl <- cutree(hc, k = k)
  mean(silhouette(cl, dist_animals)[, "sil_width"])
})
k_opt <- ks[which.max(sil_mean)]

# ---- Final groups & outputs ----
groups <- cutree(hc, k = k_opt)

cluster_table <- data.frame(
  animal  = names(groups),
  cluster = paste0("C", groups),
  row.names = NULL
) %>% arrange(cluster, animal)

print(cluster_table)
# write.csv(cluster_table, "animal_clusters_order_level.csv", row.names = FALSE)

# ---- Plot dendrogram with cluster boxes ----
plot(hc, main = sprintf("Animals clustered by ORDER-level diet (Bray–Curtis) | k = %d", k_opt),
     xlab = "", sub = "")
rect.hclust(hc, k = k_opt, border = "red")

# ---- (Optional) Heatmap, columns ordered by hc ----
# install.packages("pheatmap")
library(pheatmap)
ord <- hc$order
pheatmap::pheatmap(mat_rel[, ord, drop = FALSE],
                   cluster_rows = TRUE, cluster_cols = FALSE,
                   main = sprintf("Order-level diet composition (k = %d)", k_opt))


# ---- (Optional) k-means on PCoA of Bray–Curtis for hard partitions ----
# pcoa <- cmdscale(dist_animals, k = (ncol(mat_rel) - 1), eig = TRUE)
# coords <- as.data.frame(pcoa$points); rownames(coords) <- colnames(mat_rel)
# set.seed(123); km <- kmeans(coords, centers = k_opt, nstart = 50)
# data.frame(animal = rownames(coords), kmeans_cluster = paste0("K", km$cluster)) %>% arrange(kmeans_cluster, animal)


###################################################################################################################3
## method 2
# ---- Packages ----
library(readxl)
library(dplyr)
library(tidyr)
library(janitor)
library(vegan)     # vegdist, metaMDS, simper, decostand
library(cluster)   # pam, silhouette
library(ggplot2)

# =========================
# Parameters you may tweak:
# =========================
force_k <- 3            # set to 2 for browser–grazer; NA = choose by silhouette (2..6)
min_reads_per_sample <- 100   # drop raw sample columns with <100 total reads
min_samples_per_animal <- 3   # drop animals represented by <3 retained samples
min_reads_per_order <- 100    # drop plant orders with <100 total reads across all animals

# -----------------------------
# 1) Load data
# -----------------------------
raw <- read_excel(path_to_file, sheet = sheet_name) %>%
  clean_names()

# avoid clash with base::order
if ("order" %in% names(raw)) {
  raw <- raw %>% rename(plant_order = order)
}

# identify numeric sample columns
num_cols <- raw %>% select(where(is.numeric)) %>% names()

if (length(num_cols) == 0) stop("No numeric diet columns found.")

# helper to collapse replicate columns to animal name
base_name <- function(x) gsub("[\\._].*$", "", x)

colnames(ord_animal)

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
# -----------------------------
# 2) Drop low-read sample columns
# -----------------------------
col_totals <- raw %>%
  summarise(across(all_of(num_cols), ~ sum(.x, na.rm = TRUE))) %>%
  pivot_longer(everything(), names_to = "sample_col", values_to = "total_reads")

keep_num_cols <- col_totals %>%
  filter(total_reads >= min_reads_per_sample) %>%
  pull(sample_col)

raw <- raw %>%
  select(-all_of(setdiff(num_cols, keep_num_cols)))

num_cols <- keep_num_cols

if (length(num_cols) == 0) stop("No sample columns remain after filtering by min_reads_per_sample.")

# -----------------------------
# 3) Count samples per animal and keep only animals with enough samples
# -----------------------------
animal_n <- tibble(col = num_cols) %>%
  mutate(
    animal = base_name(col),
    animal = clean_animal_name(animal)
  ) %>%
  count(animal, name = "n_samples")

keep_animals <- animal_n %>%
  filter(n_samples >= min_samples_per_animal) %>%
  pull(animal)

if (length(keep_animals) == 0) stop("No animals remain after filtering by min_samples_per_animal.")

# -----------------------------
# 4) Aggregate to ORDER × ANIMAL
# -----------------------------
ord_animal <- raw %>%
  select(plant_order, all_of(num_cols)) %>%
  filter(!is.na(plant_order), plant_order != "") %>%
  pivot_longer(cols = all_of(num_cols), names_to = "col", values_to = "value") %>%
  mutate(
    animal = base_name(col),
    animal = clean_animal_name(animal)
  ) %>%
  filter(animal %in% keep_animals) %>%
  group_by(plant_order, animal) %>%
  summarise(value = sum(value, na.rm = TRUE), .groups = "drop") %>%
  pivot_wider(names_from = animal, values_from = value, values_fill = 0) %>%
  arrange(plant_order)

# -----------------------------
# 5) Drop low-abundance plant orders
# -----------------------------
ord_animal <- ord_animal %>%
  mutate(total_reads = rowSums(across(-plant_order), na.rm = TRUE)) %>%
  filter(total_reads >= min_reads_per_order) %>%
  select(-total_reads)

if (nrow(ord_animal) < 2) {
  stop("After filtering, <2 orders remain. Lower 'min_reads_per_order'.")
}

# -----------------------------
# 6) Rename animal columns to include number of retained samples
# -----------------------------
new_names <- animal_n %>%
  filter(animal %in% keep_animals) %>%
  mutate(new_name = paste0(animal, " (n=", n_samples, ")")) %>%
  select(animal, new_name) %>%
  tibble::deframe()

ord_animal <- ord_animal %>%
  rename_with(~ ifelse(.x %in% names(new_names), new_names[.x], .x))

# view result
ord_animal

# --------------------------------
# 3) Build matrix (orders x animals) and clean
# --------------------------------
mat <- ord_animal %>%
  column_to_rownames("plant_order") %>%
  as.matrix()

# Drop animals (columns) with zero total
col_tot <- colSums(mat, na.rm = TRUE)
if (any(col_tot == 0)) {
  message("Dropping animals with 0 total reads: ",
          paste(names(col_tot)[col_tot == 0], collapse = ", "))
  mat <- mat[, col_tot > 0, drop = FALSE]
}
if (ncol(mat) < 2) stop("Need at least 2 animals (columns) with non-zero totals.")

# --------------------------------
# 4) Relative composition + Bray–Curtis
# --------------------------------
mat_rel <- sweep(mat, 2, colSums(mat), "/")
mat_rel[is.na(mat_rel)] <- 0

bc <- vegdist(t(mat_rel), method = "bray")   # animals x animals Bray–Curtis


# --------------------------------
# 4) Composition + jaccard
# --------------------------------
mat_pa <- ifelse(mat > 0, 1, 0)

jac <- vegdist(t(mat_pa), method = "jaccard", binary = TRUE)

# --------------------------------
# 5) PAM (k-medoids) clustering on Bray–Curtis
# --------------------------------
choose_k <- function(dist_obj, k_range = 2:6) {
  sils <- sapply(k_range, function(k) {
    pam(dist_obj, k = k, diss = TRUE)$silinfo$avg.width
  })
  data.frame(k = k_range, avg_sil = sils) %>% arrange(desc(avg_sil))
}

if (is.na(force_k)) {
  sil_tbl <- choose_k(bc, k_range = 2:min(6, length(attr(bc, "Labels")) - 1))
  print(sil_tbl)
  k_opt <- sil_tbl$k[1]
} else {
  k_opt <- force_k
}

pam_fit <- pam(bc, k = k_opt, diss = TRUE)

clusters <- data.frame(
  animal  = names(pam_fit$clustering),
  cluster = paste0("C", pam_fit$clustering),
  row.names = NULL
) %>% arrange(cluster, animal)

cat("\n== PAM cluster assignments (k =", k_opt, ") ==\n")
print(clusters, row.names = FALSE)

# --------------------------------
# 6) What separates clusters? (unsupervised)
#    SIMPER between groups on the original relative matrix
# --------------------------------
grp <- factor(pam_fit$clustering)
comm_animals <- t(mat_rel)   # rows = animals, cols = orders

if (nlevels(grp) >= 2) {
  sim <- simper(comm_animals, grp, permutations = 0)
  cat("\n== SIMPER: top orders contributing to between-cluster differences ==\n")
  # Print top contributors for each pair; limit to first few lines for readability
  print(lapply(sim, function(tb) head(tb[order(-tb$average), ], 10)))
}

# --------------------------------
# 7) Quick ordination plot (visual only)
# --------------------------------
set.seed(1)
nmds <- metaMDS(comm_animals, distance = "bray", k = 2, trymax = 50, autotransform = FALSE, trace = FALSE)

scores_df <- data.frame(scores(nmds, display = "sites"))
scores_df$animal  <- rownames(scores_df)
scores_df$cluster <- factor(pam_fit$clustering)

ggplot(scores_df, aes(NMDS1, NMDS2, color = cluster, label = animal)) +
  geom_point(size = 3) +
  ggrepel::geom_text_repel(show.legend = FALSE) +
  theme_minimal() +
  labs(title = paste0("NMDS (Bray–Curtis) with PAM clusters (k=", k_opt, ")"))

##
# Reorder animals by cluster assignment
ord_cols <- clusters$animal
mat_rel_ord <- mat_rel[, ord_cols]

pheatmap(mat_rel_ord,
         cluster_rows = TRUE,    # still cluster orders
         cluster_cols = FALSE,   # keep your cluster order
         annotation_col = data.frame(Cluster = clusters$cluster,
                                     row.names = clusters$animal),
         main = "Order-level diet by PAM clusters",
         fontsize_row = 6,
         fontsize_col = 10)


##
topN <- 10  # change as needed

# Keep only animals present in mat_rel and clusters, ordered by PAM cluster
common_animals <- intersect(colnames(mat_rel), clusters$animal)
clust_df <- clusters %>%
  filter(animal %in% common_animals) %>%
  arrange(cluster, animal)
ord_cols <- clust_df$animal

# Pick top-N orders by *counts* and subset proportions for the heatmap
top_orders <- names(sort(rowSums(mat[, ord_cols, drop = FALSE]), decreasing = TRUE))
top_orders <- head(top_orders, min(topN, length(top_orders)))

mat_top <- mat_rel[top_orders, ord_cols, drop = FALSE]

# Add "Other" row so columns sum to 1
other_row <- pmax(0, 1 - colSums(mat_top))            # numeric remainder per animal
mat_top2  <- rbind(mat_top, Other = other_row)

# Create % labels; force each column to sum exactly 100 by assigning rounding diff to "Other"
lab <- round(100 * mat_top2)                           # initial rounded %
for (j in seq_len(ncol(lab))) {
  diff <- 100 - sum(lab[, j], na.rm = TRUE)           # could be -2..+2 typically
  lab["Other", j] <- lab["Other", j] + diff           # adjust "Other" to make exact 100
}
num_labels <- matrix(
  sprintf("%d%%", lab),
  nrow = nrow(lab), ncol = ncol(lab),
  dimnames = dimnames(lab)
)

# Annotation & gaps between clusters
ann_col <- data.frame(Cluster = clust_df$cluster, row.names = clust_df$animal)
cluster_sizes <- as.integer(table(clust_df$cluster))
gaps_col <- cumsum(cluster_sizes); gaps_col <- gaps_col[-length(gaps_col)]

# Heatmap (columns fixed by PAM; rows clustered for readability)
pheatmap(mat_top2,
         cluster_rows = F,
         cluster_cols = FALSE,
         annotation_col = ann_col,
         gaps_col = gaps_col,
         display_numbers = num_labels,
         number_color = "black",
         main = sprintf("Top %d orders (+ Other) by diet composition (PAM order)", topN),
         fontsize_row = 7,
         fontsize_col = 10,
         fontsize_number = 6)

library(ggrepel)
library(gridExtra)
library(grid)
graphics.off()
# --- build heatmap as a grob ---
cluster_cols <- c(
  "1" = "#3A7D44",
  "2" = "#9C6644",
  "3" = "#7570b3"
)
ann_colors <- list(
  Cluster = c(
    "C1" = "#3A7D44",
    "C2" = "#9C6644",
    "C3" = "#7570b3"
  )
)

ph <- pheatmap(
  mat_top2,
  cluster_rows = FALSE,
  cluster_cols = FALSE,
  annotation_col = ann_col,
  annotation_colors = ann_colors,
  gaps_col = gaps_col,
  display_numbers = num_labels,
  number_color = "black",
  fontsize_row = 15,
  fontsize_col = 18,
  fontsize_number = 12,
  silent = TRUE
)


p_nmds <- ggplot(scores_df, aes(NMDS1, NMDS2, color = cluster, label = animal)) +
  geom_point(size = 5) +
  ggrepel::geom_text_repel(size = 5, show.legend = FALSE) +
  scale_color_manual(values = cluster_cols) +
  theme_minimal()


gA <- arrangeGrob(textGrob("A", x = 0.02, y = 0.98, just = c("left", "top"),
                           gp = gpar(fontsize = 16, fontface = "bold")),
                  ph$gtable,
                  ncol = 1,
                  heights = c(0.05, 0.95))

gB <- arrangeGrob(textGrob("B", x = 0.02, y = 0.98, just = c("left", "top"),
                           gp = gpar(fontsize = 16, fontface = "bold")),
                  ggplotGrob(p_nmds),
                  ncol = 1,
                  heights = c(0.05, 0.95))

grid.arrange(gA, gB, ncol = 2, widths = c(1.2, 1))

ggsave("diet_heatmap_nmds_landscape - BC.pdf", arrangeGrob(gA, gB, ncol = 2, widths = c(1.2, 1)), width = 16, height = 9, units = "in")
##################################################################################################################

adonis2(t(mat_rel) ~ guild_df$guild, method = "bray")
adonis2(t(mat_pa) ~ guild_df$guild, method = "jaccard", binary = TRUE)