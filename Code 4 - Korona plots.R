library(readxl)
library(dplyr)
library(tidyr)
library(janitor)
library(stringr)
library(purrr)
library(tibble)
library(webshot)
library(webshot2)
library(pagedown)

webshot::install_phantomjs()

setwd("C:/Users/Admin/Desktop/R/Rose et al")

Sys.which("ktImportXML")

Sys.setenv(PATH = paste("C:/Tools/krona_tools/bin", Sys.getenv("PATH"), sep = ";"))
Sys.which("ktImportXML")

kronatools_dir <- "C:/Tools/krona_tools"

# -----------------------------
# 1) Paths and settings
# -----------------------------
path_darius  <- "Darius_3_formatted_hit_May15_cleanedDS.xlsx"
sheet_darius <- "Darius_3_formatted_hit_May15_KR"

min_reads_per_sample   <- 100   # drop sample columns with <100 total reads
min_samples_per_animal <- 3     # keep animals with >=3 retained samples
min_reads_per_taxon    <- 1     # optional: drop tiny taxa within an animal Krona
animal_summary_fn      <- sum   # use sum to combine samples per animal
# animal_summary_fn    <- mean  # switch to mean if you want average composition

out_dir <- "krona_by_animal"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# -----------------------------
# 2) Load embed_krona function
# -----------------------------
source("https://raw.githubusercontent.com/markschl/embed_krona/main/embed_krona.R")

# If ktImportXML is not on your PATH, set this instead:
kronatools_dir <- "C:/Tools/krona_tools"

# -----------------------------
# 3) Helper functions
# -----------------------------
base_name <- function(x) {
  gsub("[\\._].*$", "", x)
}

clean_animal_name <- function(x) {
  x <- tolower(trimws(x))
  
  dplyr::recode(
    x,
    "bontebok"            = "Bontebok",
    "buffalo"             = "Buffalo",
    "zebra"               = "Burchell's zebra",
    "burchellszebra"      = "Burchell's zebra",
    "burchell's zebra"    = "Burchell's zebra",
    "capezebra"           = "Cape zebra",
    "cape zebra"          = "Cape zebra",
    "eland"               = "Eland",
    "elephant"            = "Elephant",
    "giraffe"             = "Giraffe",
    "impala"              = "Impala",
    "white"               = "White Rhino",
    "rhino"               = "White Rhino",
    "white_rhino"         = "White Rhino",
    "white_rhinoceros"    = "White Rhino",
    "springbok"           = "Springbok",
    "waterbuck"           = "Waterbuck",
    "wildebeest"          = "Wildebeest",
    .default = str_to_title(x)
  )
}

standardise_animal_name <- function(x) {
  x <- trimws(x)
  
  dplyr::recode(
    x,
    "Burchell's zebra" = "Burchell's zebra",
    "White Rhino"      = "White rhinoceros",
    "Rhino"            = "White rhinoceros",
    "White"            = "White rhinoceros",
    "White rhino"      = "White rhinoceros",
    "Plains Zebra"     = "Burchell's zebra",
    "Cape Zebra"       = "Cape zebra",
    .default = x
  )
}

normalise_animal_names <- function(x) {
  x %>%
    base_name() %>%
    clean_animal_name() %>%
    standardise_animal_name()
}

safe_taxon <- function(x, unknown = "Unknown") {
  x <- trimws(as.character(x))
  x[x == "" | is.na(x)] <- unknown
  x
}

safe_file_name <- function(x) {
  x %>%
    str_replace_all("[^A-Za-z0-9]+", "_") %>%
    str_replace_all("^_+|_+$", "")
}

# -----------------------------
# 4) Load raw data
# -----------------------------
raw <- read_excel(path_darius, sheet = sheet_darius) %>%
  clean_names()

# rename taxonomy columns if needed
if ("order" %in% names(raw))   raw <- raw %>% rename(plant_order = order)
if ("family" %in% names(raw))  raw <- raw %>% rename(plant_family = family)
if ("species" %in% names(raw)) raw <- raw %>% rename(plant_species = species)

# identify numeric sample columns
num_cols <- raw %>%
  select(where(is.numeric)) %>%
  names()

if (length(num_cols) == 0) {
  stop("No numeric diet columns found.")
}

# -----------------------------
# 5) Drop low-read sample columns
# -----------------------------
col_totals <- raw %>%
  summarise(across(all_of(num_cols), ~ sum(.x, na.rm = TRUE))) %>%
  pivot_longer(everything(), names_to = "sample_col", values_to = "total_reads")

keep_num_cols <- col_totals %>%
  filter(total_reads >= min_reads_per_sample) %>%
  pull(sample_col)

if (length(keep_num_cols) == 0) {
  stop("No sample columns remain after min_reads_per_sample filtering.")
}

# -----------------------------
# 6) Keep animals with enough retained samples
# -----------------------------
animal_n <- tibble(sample_col = keep_num_cols) %>%
  mutate(animal = normalise_animal_names(sample_col)) %>%
  count(animal, name = "n_samples") %>%
  arrange(desc(n_samples))

print(animal_n)

keep_animals <- animal_n %>%
  filter(n_samples >= min_samples_per_animal) %>%
  pull(animal)

if (length(keep_animals) == 0) {
  stop("No animals remain after min_samples_per_animal filtering.")
}

# -----------------------------
# 7) Build long table at full taxonomy depth
# -----------------------------
tax_long <- raw %>%
  select(
    plant_order,
    plant_family,
    plant_species,
    all_of(keep_num_cols)
  ) %>%
  pivot_longer(
    cols = all_of(keep_num_cols),
    names_to = "sample_col",
    values_to = "reads"
  ) %>%
  mutate(
    animal       = normalise_animal_names(sample_col),
    plant_order  = safe_taxon(plant_order,  "Unknown_order"),
    plant_family = safe_taxon(plant_family, "Unknown_family"),
    plant_species= safe_taxon(plant_species,"Unknown_species")
  ) %>%
  filter(
    animal %in% keep_animals,
    !is.na(reads),
    reads > 0
  )

# -----------------------------
# 8) Summarise to one profile per animal
# -----------------------------
animal_taxa <- tax_long %>%
  group_by(animal, plant_order, plant_family, plant_species) %>%
  summarise(
    abundance = animal_summary_fn(reads, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  filter(abundance >= min_reads_per_taxon)

# -----------------------------
# 9) Make one Krona chart per animal
# -----------------------------
animals <- sort(unique(animal_taxa$animal))

for (sp in animals) {
  
  dat_sp <- animal_taxa %>%
    filter(animal == sp) %>%
    arrange(desc(abundance), plant_order, plant_family, plant_species)
  
  if (nrow(dat_sp) == 0) next
  
  classification <- as.matrix(
    dat_sp[, c("plant_order", "plant_family", "plant_species")]
  )
  
  magnitude <- matrix(dat_sp$abundance, ncol = 1)
  colnames(magnitude) <- "abundance"
  
  out_file <- file.path(
    out_dir,
    paste0(safe_file_name(sp), "_krona.html")
  )
  
  plot_krona(
    classification = classification,
    magnitude      = magnitude,
    output         = out_file,
    display        = FALSE,
    root_label     = sp,
    total_label    = "Total reads",
    snapshot       = FALSE
  )
  
  message("Saved: ", out_file)
}

# -----------------------------
# 10) Optional: save the animal x taxon table used for plotting
# -----------------------------
write.csv(
  animal_taxa,
  file.path(out_dir, "animal_taxa_used_for_krona.csv"),
  row.names = FALSE
)

###############################################################################
# Testing 1 animal

kronatools_dir <- "C:/Tools/Krona/KronaTools"

sp <- unique(animal_taxa$animal)[1]

dat_sp <- animal_taxa %>%
  filter(animal == sp) %>%
  arrange(desc(abundance), plant_order, plant_family, plant_species)

classification <- dat_sp %>%
  select(plant_order, plant_family, plant_species) %>%
  as.matrix()

magnitude <- matrix(dat_sp$abundance, ncol = 1)
colnames(magnitude) <- "abundance"

out_file <- file.path(out_dir, paste0(safe_file_name(sp), "_krona.html"))

plot_krona(
  classification = classification,
  magnitude      = magnitude,
  output         = out_file,
  display        = FALSE,
  root_label     = sp,
  total_label    = "Total reads",
  snapshot       = FALSE
)

#
Sys.which("ktImportXML")
system2("ktImportXML", stdout = TRUE, stderr = TRUE)
plot_krona(
  classification = classification,
  magnitude      = magnitude,
  output         = out_file,
  display        = FALSE,
  root_label     = sp,
  total_label    = "Total reads",
  snapshot       = FALSE
)


##################################################################
# Loop with export as PNG
png_dir <- file.path(out_dir, "png")
dir.create(png_dir, showWarnings = FALSE, recursive = TRUE)

animals <- sort(unique(animal_taxa$animal))

for (sp in animals) {
  
  dat_sp <- animal_taxa %>%
    filter(animal == sp) %>%
    arrange(desc(abundance), plant_order, plant_family, plant_species)
  
  if (nrow(dat_sp) == 0) next
  
  classification <- as.matrix(
    dat_sp[, c("plant_order", "plant_family", "plant_species")]
  )
  
  magnitude <- matrix(dat_sp$abundance, ncol = 1)
  colnames(magnitude) <- "abundance"
  
  out_file_png <- file.path(
    png_dir,
    paste0(safe_file_name(sp), "_krona.png")
  )
  
  plot_krona(
    classification      = classification,
    magnitude           = magnitude,
    output              = out_file_png,
    display             = FALSE,
    root_label          = sp,
    total_label         = "Total reads",
    snapshot            = TRUE,
    snapshot_format     = "png",
    snapshot_fontsize   = 13,
    snapshot_chart_size = 0.7,     # slightly smaller inside canvas
    snapshot_dim        = c(13, 12),
    snapshot_res        = 300
  )
  
  message("Saved: ", out_file_png)
}

#
sp <- animals[1]

dat_sp <- animal_taxa %>%
  filter(animal == sp) %>%
  arrange(desc(abundance), plant_order, plant_family, plant_species)

classification <- as.matrix(
  dat_sp[, c("plant_order", "plant_family", "plant_species")]
)

magnitude <- matrix(dat_sp$abundance, ncol = 1)
colnames(magnitude) <- "abundance"

out_file_png <- file.path(
  png_dir,
  paste0(safe_file_name(sp), "_krona.png")
)

plot_krona(
  classification      = classification,
  magnitude           = magnitude,
  output              = out_file_png,
  display             = FALSE,
  root_label          = sp,
  total_label         = "Total reads",
  snapshot            = TRUE,
  snapshot_format     = "png",
  snapshot_fontsize   = 11,
  snapshot_chart_size = 1.2,
  snapshot_dim        = c(10, 8),
  snapshot_res        = 300
)

#############################################
# PDF export
install.packages("pagedown")


html_files <- list.files(out_dir, pattern = "\\.html$", full.names = TRUE)

for (f in html_files) {
  pdf_file <- sub("\\.html$", ".pdf", f)
  
  chrome_print(
    input  = normalizePath(f),
    output = pdf_file
  )
  
  message("Saved: ", pdf_file)
}


############################
# PNG export

html_files <- list.files(out_dir, pattern = "\\.html$", full.names = TRUE)

for (f in html_files) {
  png_file <- sub("\\.html$", ".png", f)
  
  webshot(
    url = paste0("file:///", normalizePath(f, winslash = "/")),
    file = png_file,
    vwidth = 1600,
    vheight = 1200,
    zoom = 2
  )
  
  message("Saved: ", png_file)
}
