library(ggpubr)
library(rbiom)
library(glue)
library(vegan)
library(phyloseq)
library(tidyverse)

######################################################################################
########## Distinct characteristics of early-stage PD-RBD(+) gut microbiome ##########
######################################################################################

# figures 1A, B
# Examining the beta-diversity between early-stage PD patients and HC
load("16S.RData") # Load 16S phyloseq
load("WGS.RData") # Load WGS phyloseq

ps <- subset_samples(ps, D_duration_y < 2 | is.na(D_duration_y))
set.seed(111) # Set seed
ord.cols <- c(HC = "#F2D57E", `PD-RBD(-)` = "#F21D56", `PD-RBD(+)` = "#1F4C73")

# Rarefy table to minimum depth
# Skip if data type is shotgun metagenom
min.depth <- min(sample_sums(ps))
ps <-  rarefy_even_depth(ps, min.depth, replace = FALSE, rngseed = 111)

# Calcaulte distance and dimensionality reduction
d <- unifrac(ps@otu_table, weighted = FALSE, ps@phy_tree)
dm <- stats::cmdscale(d, k = nrow(as.matrix(d))-1, eig = TRUE)

# Draw PCoA
pts <- data.frame(as.matrix(dm$points))
names(pts) <- paste0("Axis", seq(1, ncol(pts)))
pts <- pts[,c("Axis1", "Axis2")]
df <- data.frame(pts, sample_data(ps))

# Calculate centroids
cent <- aggregate(reformulate("RBD", response = paste0("cbind(Axis1, Axis2)")), data = df, FUN = mean) %>%
  setNames(c("RBD", paste0("end_Axis1"), paste0("end_Axis2")))
df <- merge(df, cent, by = "RBD", sort = FALSE)
    
# Calculate relative eigen values
eig <- dm$eig
rel.eig <- eig / sum(eig)
rel.eig <- rel.eig[c(1, 2)]

# Axis parameters
axis.labels <- round(rel.eig, 4)*100
axis.labels <- paste0("PCo ", c(1, 2), " (", axis.labels, " %)")

# Group labeling
labs <- table(df$RBD)
labs <- str_c(names(labs), " (n = ", labs, ")")
names(labs) <- c("HC", "PD-RBD(-)", "PD-RBD(+)")

# Draw
p <- ggplot(df, aes(x = Axis1, y = Axis2)) + 
  geom_hline(yintercept = 0, linetype = "dashed", color = "black", linewidth = .5) + 
  geom_vline(xintercept = 0, linetype = "dashed", color = "black", linewidth = .5) +
  stat_ellipse(aes(group = RBD, color = RBD), lwd = 0.5, linetype = 2, show.legend = FALSE, segments = 101) +
  geom_point(aes(color = RBD), size = 2.5) +
  geom_segment(aes(xend = end_Axis1, yend = end_Axis2, color = RBD), lwd = 0.5, alpha = 0.75) +
  geom_point(data = cent, aes(x = end_Axis1, y = end_Axis2, fill = RBD), size = 3.5, shape = 21, 
             color = "black", stroke = 0.75, alpha = 0.75) + 
  theme_bw() + 
  labs(x = axis.labels[1], y = axis.labels[2]) +
  theme(
    panel.border = element_rect(linewidth = 1),
    text = element_text(size = 18, face = "bold", color = "black"),
    axis.text = element_text(size = 18, color = "black"),
    panel.grid = element_blank(),
    legend.text = element_text(face = "bold", size = 20),
    legend.title = element_text(face = "bold", size = 20),
    legend.key.size = unit(2, "line")
  ) + scale_y_continuous(expand = expansion(mult = 0.15)) + 
  scale_x_continuous(expand = expansion(mult = 0.15)) + 
  scale_fill_manual(guide = NULL, values = ord.cols) +
  scale_color_manual(name = "Group", values = ord.cols, breaks = names(ord.cols), labels = labs) +
  guides(color = guide_legend(override.aes = list(shape = 16, alpha = 1, linetype = "blank", size =6, fill = ord.cols)))

########################################################################
########## Gut microbiome of patients with PD according to    ##########
########## the presence of premotor RBD and stage of the disease  ######
########################################################################

# Figures 1C, D
# Examining the beta-diversity between PD patients and HC, considering the duration of PD.
load("16S.RData") # Load 16S phyloseq
load("WGS.RData") # Load WGS phyloseq
load("utilities.RData") # Load Function

# Rarefy table to minimum depth
# Skip if data type is shotgun metagenom
min.depth <- min(sample_sums(ps))
ps <-  rarefy_even_depth(ps, min.depth, replace = FALSE, rngseed = 111)

# Calcaulte distance and dimensionality reduction
d <- unifrac(ps@otu_table, weighted = FALSE, ps@phy_tree)
dm <- stats::cmdscale(d, k = nrow(as.matrix(d))-1, eig = TRUE)

# Draw PCoA
pts <- data.frame(as.matrix(dm$points))
names(pts) <- paste0("Axis", seq(1, ncol(pts)))
pts <- pts[,c("Axis1", "Axis2")]
df <- data.frame(pts, sample_data(ps))

# Calculate centroids
cent <- aggregate(reformulate("RBD", response = paste0("cbind(Axis1, Axis2)")), data = df, FUN = mean) %>%
  setNames(c("RBD", paste0("end_Axis1"), paste0("end_Axis2")))
df <- merge(df, cent, by = "RBD", sort = FALSE)

max.d.duration <- max(df$D_duration_y, na.rm = TRUE)+1

min.year.list <- c(0, 2, "all")
max.year.list <- c(2, max.d.duration, "all")
year.comb <- rbind(min.year.list, max.year.list)

lres.df <- list()
for (idx in 1:ncol(year.comb)){
  min_year <- year.comb[1, idx]
  max_year <- year.comb[2, idx]
  
  if (min_year == "all") min_year <- NULL
  if (max_year == "all") max_year <- NULL
  
  lres.df[[idx]] <- calc.cent(df, min_year, max_year)
  
}

mean.cl.df <- bind_rows(lres.df)
mean.cl.df <- mean.cl.df %>% 
  mutate(
    year.cut = str_c(min_year, "≤", RBD , "<", max_year,  sep = ""),
    year.cut = ifelse(RBD == "HC", "HC", year.cut),
    year.cut = ifelse(str_detect(year.cut, "All"), "All", year.cut),
    year.cut = ifelse(str_detect(year.cut, as.character(max.d.duration)), 
                      str_replace(year.cut, paste0("<", max.d.duration), ""), year.cut),
    year.cut = ifelse(str_detect(year.cut, "0≤"), 
                      str_replace(year.cut, "0≤", ""), year.cut)
  )

mean.cl.df <- mean.cl.df %>% filter(year.cut != "All")

RBD_minus <- c("#F21D56", "#74162F")
RBD_plus <- c("#1F4C73", "#283B56")
RBD_combined <- c(rbind(RBD_plus, RBD_minus))
if ("HC" %in% df$RBD){
  HC.color <- "#F2D57E"
} else{
  HC.color <- NULL
}

gradual.colors <- c(RBD_combined, HC.color)
mean.cl.df$colors <- gradual.colors

mean.cl.df <- mean.cl.df %>% 
  mutate(cur.group = ifelse(year.cut == "All", paste(year.cut, RBD, sep = " "), year.cut)) %>% 
  mutate(cur.group = ifelse(RBD == "HC", "HC", cur.group)) %>% 
  mutate(cur.group = glue("{cur.group} (n={n})"))
group.levels <- c(str_subset(mean.cl.df$cur.group, "PD\\-RBD\\(\\+\\)"), 
                  str_subset(mean.cl.df$cur.group, "PD\\-RBD\\(\\-\\)"),
                  str_subset(mean.cl.df$cur.group, "HC"))
mean.cl.df$cur.group <- factor(mean.cl.df$cur.group, levels = group.levels)
mean.cl.df <- mean.cl.df %>% arrange(cur.group)

# Calculate relative eigen values
eig <- dm$eig
rel.eig <- eig / sum(eig)
rel.eig <- rel.eig[c(1, 2)]

# Axis parameters
axis.labels <- round(rel.eig, 4)*100
axis.labels <- paste0("PCo ", c(1, 2), " (", axis.labels, " %)")

# Draw
p <- ggplot(mean.cl.df, aes(x = Axis1_median, y = Axis2_median, color = cur.group, 
                            shape = cur.group)) + 
  geom_hline(yintercept = 0, linetype = "dashed", color = "black", linewidth = .5) + 
  geom_vline(xintercept = 0, linetype = "dashed", color = "black", linewidth = .5) +
  geom_point(size = 6) +
  geom_errorbar(aes(ymin = Axis2_min, ymax = Axis2_max), width = 0) + 
  geom_errorbar(aes(xmin = Axis1_min, xmax = Axis1_max), width = 0) + 
  scale_color_manual(name = "Groups", labels = mean.cl.df$cur.group, 
                     values = mean.cl.df$colors) +
  scale_shape_manual(name = "Groups", labels = mean.cl.df$cur.group, 
                     values = c(16, 16, 17, 17, 15)) + 
  theme_bw() + 
  labs(x = axis.labels[1], y = axis.labels[2]) +
  theme(
    panel.border = element_rect(linewidth = 1),
    text = element_text(size = 20, color = "black"),
    axis.text = element_text(size = 20, color = "black"),
    axis.title = element_text(face = "bold"),
    panel.grid = element_blank(),
    legend.text = element_text(size = 20, color = "black"),
    legend.title = element_text(size = 20, face = "bold", color = "black"),
    legend.key.size = unit(2, "line")
  ) + scale_y_continuous(expand = expansion(mult = 0.15)) + 
  scale_x_continuous(expand = expansion(mult = 0.15))


############################################################################
########## Difference gut microbiome dynamics in patients with PD ##########
##########       according to the presence of premotor RBD        ##########       
############################################################################

# Figures 1G-J
load("16S.RData") # Load 16S phyloseq
load("WGS.RData") # Load WGS phyloseq
load("utilities.RData") # Load Function

# Rarefy table to minimum depth
# Skip if data type is shotgun metagenom
min.depth <- min(sample_sums(ps))
ps <-  rarefy_even_depth(ps, min.depth, replace = FALSE, rngseed = 111)

# Calcaulte distance and dimensionality reduction
d <- unifrac(ps@otu_table, weighted = FALSE, ps@phy_tree)
md <- data.frame(ps@sam_data)
md$sample_id <- rownames(md)
cols <- c(HC = "#F2D57E", `PD-RBD(-)` = "#F21D56", `PD-RBD(+)` = "#1F4C73")

# Current distance to melted format
m <- as.matrix(d)
m <- data.frame(m)
m <- m %>% rownames_to_column("id.1")
df <- m %>% pivot_longer(!id.1, names_to = "id.2", values_to = "distance")

# Remove duplicate distances originated from symmetric distance matrix
df <- df %>% filter(id.1 != id.2)
df <- df %>% rowwise %>% mutate(id = str_c(sort(c(id.1, id.2)), collapse = "vs"))
df <- df %>% distinct(id, distance) %>% ungroup
df <- df %>% mutate(id.1 = str_extract(id, "(^.+?)vs", group = 1),
                              id.2 = str_extract(id, "vs(.+$)", group = 1))
df <- df %>% select(-id)
df <- df %>% left_join(md, by = join_by(id.1 == sample_id)) %>% 
  left_join(md, by = join_by(id.2 == sample_id), suffix = c(".1", ".2"))

# Make column representing the distance type
var.order <- c("HC", "PD-RBD(-)", "PD-RBD(+)")
df <- df %>% rowwise() %>% mutate(type = check.type(RBD.1, RBD.2)) %>% ungroup()

# Figures 1G, H
gdf <- df %>% filter(!str_detect(type, "^Within"))
gdf <- gdf %>% filter(str_detect(type, "^HC"))
gdf <- gdf %>% mutate(D_duration_y = coalesce(D_duration_y.1, D_duration_y.2),
                      RBD = case_when(
                        RBD.1 == "HC" ~ RBD.2,
                        RBD.2 == "HC" ~ RBD.1,
                        .default = NA))

p <- ggplot(gdf, aes(x = ifelse(RBD == "PD-RBD(+)", D_duration_y+0.1, D_duration_y-0.1),
                     y = distance, color = RBD)) + 
    geom_jitter(width = 0.02, alpha = 0.7, shape = 21) +
    geom_smooth(method = "lm") +
    stat_cor(method = "pearson", cor.coef.name = "r", 
             geom = "label", alpha = 0.5, show.legend = FALSE,
             label.x.npc = 0.4, label.y.npc = 0.2) +
    theme_bw() +
    scale_color_manual(name = "Group", values = cols) +
    labs(x = "Disease duration(year)", y = "Weighted UniFrac distance to HC") +
    theme(
      axis.title = element_text(color = "black", face = "bold"),
      legend.title = element_text(color = "black", face = "bold"),
      axis.text = element_text(color = "black")
    )

  
# Figures 1I, J
gdf <- df %>% filter(type %in% c("Within PD-RBD(+)", "Within PD-RBD(-)"))
gdf <- gdf %>% mutate(stage.1 = ifelse(D_duration_y.1 < 2, "Early", "Late"),
                            stage.2 = ifelse(D_duration_y.2 < 2, "Early", "Late"))
gdf <- gdf %>% filter(stage.1 != stage.2)
gdf <- gdf %>% rowwise %>%
  mutate(D_duration_y = max(D_duration_y.1, D_duration_y.2)) %>% ungroup

p <- ggplot(gdf, aes(x = ifelse(RBD.2 == "PD-RBD(+)", as.numeric(D_duration_y)+0.1, as.numeric(D_duration_y)-0.1), y = distance, color = RBD.1)) + 
  geom_jitter(width = 0.02, alpha = 0.7, shape = 21) +
  geom_smooth(method = "lm") +
  stat_cor(method = "pearson", cor.coef.name = "r", 
           geom = "label", alpha = 0.5, show.legend = FALSE,
           label.x.npc = 0.4, label.y.npc = 0.2) +
  theme_bw() +
  scale_color_manual(name = "Group", values = cols) +
  scale_x_continuous(breaks = seq(2, max(gdf$D_duration_y.2), by = 2)) +
  labs(x = "Disease duration(year)", y = glue("Unweighted UniFrac distance\nfrom <2 years")) +
  theme(
    axis.title = element_text(color = "black", face = "bold"),
    legend.title = element_text(color = "black", face = "bold"),
    axis.text = element_text(color = "black")
  )
