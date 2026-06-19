library(readr)
library(dplyr)
library(ComplexHeatmap)
library(circlize)
library(grid)

# -------------------------
# 1) Load data
# -------------------------
df <- read_csv("For Heatmap plot2.csv", show_col_types = FALSE)

# -------------------------
# 2) Numeric feature matrix
# -------------------------
num_df <- df %>% select(where(is.numeric))
if (ncol(num_df) == 0) stop("No numeric columns found in the CSV.")
mat <- as.matrix(num_df)

# -------------------------
# 3) Scale + clamp (publication-friendly)
# -------------------------
mat_z <- scale(mat)

# Clamp to reduce dominance of outliers (good for printed interpretation)
# Use [-2.5, 2.5] for smoother print contrast; change to [-2, 3] if you insist.
mat_z <- pmax(pmin(mat_z, 2.5), -2.5)

# -------------------------
# 4) Publication heatmap palette (diverging, neutral midpoint)
#    - Blue -> light neutral -> red
#    - Midpoint is near-white (but not washed out)
# -------------------------
col_fun <- colorRamp2(
  c(-2.5, 0, 2.5),
  c("#2B6CB0", "#F7F7F7", "#C53030")  # deep blue, neutral light gray, deep red
)

# -------------------------
# 5) Left annotations with per-row grid
# -------------------------
n <- nrow(df)

# Tier always Strong (soft purple for print)
tier_vec  <- rep("Strong", n)
tier_cols <- c("Strong" = "#C9A7FF")  # softer/lighter purple than neon

# Group (if exists)
has_group <- "Group" %in% names(df)
group_vec <- if (has_group) as.character(df$Group) else rep("P_A3K", n)

# Print-friendly group colors (still green/salmon/teal but slightly muted)
group_cols <- c(
  "P_A3K"  = "#4E9A06",  # green (muted)
  "P_A30K" = "#E76F51",  # salmon (muted)
  "P_P30K" = "#00A6A6"   # teal (muted)
)

# Grid line style (subtle for print)
grid_gp <- gpar(col = "#B0B0B0", lwd = 0.7)

anno_tier <- anno_simple(
  tier_vec,
  which = "row",
  col = tier_cols,
  border = TRUE,
  gp = grid_gp
)

anno_group <- anno_simple(
  group_vec,
  which = "row",
  col = group_cols,
  border = TRUE,
  gp = grid_gp
)

row_anno <- rowAnnotation(
  Tier  = anno_tier,
  Group = anno_group,
  annotation_name_gp = gpar(fontsize = 12, fontface = "bold"),
  show_annotation_name = TRUE,
  simple_anno_size = unit(6, "mm")
)

# -------------------------
# 6) Heatmap (publication layout)
# -------------------------
ht <- Heatmap(
  mat_z,
  name = "Z",
  col = col_fun,
  left_annotation = row_anno,

  cluster_rows = TRUE,
  cluster_columns = TRUE,

  show_row_names = FALSE,
  show_column_names = TRUE,
  column_names_rot = 90,
  column_names_gp = gpar(fontsize = 11, fontface = "bold"),

  # Clean title
  column_title = "Heatmap of Strong AMPs (z-scored features)",
  column_title_gp = gpar(fontsize = 18, fontface = "bold"),

  # Subtle grid on heatmap cells (print-friendly)
  rect_gp = gpar(col = "#B0B0B0", lwd = 0.55),

  # Legend tuned for print
  heatmap_legend_param = list(
    at = c(-2.5, -1.5, -0.5, 0, 0.5, 1.5, 2.5),
    title = "Z-score",
    labels_gp = gpar(fontsize = 10),
    title_gp = gpar(fontsize = 11, fontface = "bold")
  ),

  # Keep vectors crisp (important for print/PDF)
  use_raster = FALSE
)

# -------------------------
# 7) Save outputs
# -------------------------
# High-res PNG (good for Word/PowerPoint)
png("Heatmap_publication.png", width = 3600, height = 2400, res = 300)
draw(ht, heatmap_legend_side = "right", annotation_legend_side = "right")
dev.off()

# Optional: PDF (best for journals; fully vector)
pdf("Heatmap_publication.pdf", width = 12, height = 8)
draw(ht, heatmap_legend_side = "right", annotation_legend_side = "right")
dev.off()

message("Saved: Heatmap_publication.png and Heatmap_publication.pdf")
