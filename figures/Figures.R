#Figures.R

# ---- packages ----------------------------------------------------------------
suppressPackageStartupMessages({
library(tidyverse)
library(patchwork)
library(rnaturalearth)
library(rnaturalearthdata)
library(sf)
library(RColorBrewer)
library(ggnewscale)
library(scales)
})

# ---- set working directory ---------------------------------------------------
setwd("") 

###########
## Fig 1 ##
###########

##panel 1

# ---- map data ----------------------------------------------------------------
world_ll  <- ne_countries(scale = "medium", returnclass = "sf")
europe_ll <- world_ll %>% filter(region_un == "Europe")

highlight_countries <- c("Norway", "Finland", "United Kingdom", "Netherlands", "Germany")
europe_ll <- europe_ll %>%
  mutate(fill_country = ifelse(name %in% highlight_countries, "Highlighted", "Other"))

#project to LAEA Europe
europe <- st_transform(europe_ll, 3035)

# ---- adjusted label positions ------------------------------------------------
cohorts <- tibble::tribble(
  ~cohort,    ~lon,    ~lat,
  "ALSPAC",   -2.5,    51.6,  
  "GenR",      4.4,    52.1,
  "MCS",      -1.0,    54.5,
  "MoBa",      9.0,    61.4,  
  "NTR",       5.4,    53.5,  
  "POSEIDON", 11.7,    51.8,
  "PREDO",    26.0,    62.0,  
  "TwinLife",  8.3,    49.0
)

#alphabetical legend/order
cohorts$cohort <- factor(cohorts$cohort, levels = sort(unique(cohorts$cohort)))
pts_label_ll <- st_as_sf(cohorts, coords = c("lon","lat"), crs = 4326, remove = FALSE)
pts_label <- st_transform(pts_label_ll, crs_eu)

# ---- palette -----------------------------------------------------------------
pal_cohorts <- c(
  "#2c7fb8",  # 1 ALSPAC
  "#d95f02",  # 2 GenR
  "#7570b3",  # 3 MCS
  "#1b9e9e",  # 4 MoBa
  "#66a61e",  # 5 NTR
  "#e6ab02",  # 6 POSEIDON
  "#a6761d",  # 7 PREDO
  "#8c2d04"   # 8 TwinLife
)

# --- crop ---------------------------------------------------------------------
crop_ll  <- st_as_sfc(st_bbox(c(
  xmin = -22,
  xmax =  25,
  ymin =  45.8,
  ymax =  72
), crs = st_crs(4326)))
crop_prj <- st_transform(crop_ll, crs_eu) |> st_bbox()

# ---- plot map ----------------------------------------------------------------
ggplot() +
  geom_sf(data = europe, aes(fill = fill_country), color = "grey30", linewidth = 0.25, show.legend = FALSE) +
  scale_fill_manual(values = c(Highlighted = "grey70", Other = "grey90")) +
  ggnewscale::new_scale_fill() +
  geom_sf_label(
    data = pts_label,
    aes(label = cohort, fill = cohort),
    color = "white", label.size = 0.35,
    size = 7.25, 
    label.padding = unit(3.8, "pt"),
    label.r = unit(2, "pt"),
    fontface = "bold"
  ) +
  scale_fill_manual(values = pal_cohorts) +
  coord_sf(
    xlim = c(crop_prj["xmin"], crop_prj["xmax"]),
    ylim = c(crop_prj["ymin"], crop_prj["ymax"]),
    expand = FALSE, crs = st_crs(crs_eu), clip = "on"
  ) +
  labs(x = NULL, y = NULL) +
  theme_minimal() +
  theme(
    axis.text  = element_blank(),
    axis.ticks = element_blank(),
    axis.title = element_blank(),
    panel.grid = element_blank(),
    legend.position = "none",
    plot.title   = element_blank()
  )

##panel 2

# ---- data --------------------------------------------------------------------
data <- tibble::tribble(
  ~Cohort,        ~Group,                   ~N,
  "MoBa",         "Preschool age (1–5)",    23505,
  "GenR",         "Preschool age (1–5)",     1599,
  "MCS",          "Preschool age (1–5)",     2243,
  "NTR",          "Preschool age (1–5)",     2094,
  "POSEIDON",     "Preschool age (1–5)",      221,
  "PREDO",        "Preschool age (1–5)",      249,
  "TwinLife",     "Preschool age (1–5)",      160,
  "MoBa",         "School age (6–12)",        896,
  "GenR",         "School age (6–12)",       1699,
  "MCS",          "School age (6–12)",       2365,
  "NTR",          "School age (6–12)",       1761,
  "POSEIDON",     "School age (6–12)",        137,
  "PREDO",        "School age (6–12)",        230,
  "ALSPAC",       "School age (6–12)",        996
)

# ---- alphabetical cohort order -----------------------------------------------
cohort_levels <- c("ALSPAC","GenR","MCS","MoBa","NTR","POSEIDON","PREDO","TwinLife")
data <- data %>%
  mutate(
    Cohort = factor(Cohort, levels = cohort_levels),
    Group  = factor(Group, levels = c("Preschool age (1–5)", "School age (6–12)"))
  )

# ---- palette consistent with panel 1 -----------------------------------------
dark2 <- brewer.pal(8, "Dark2")
dark2[4] <- "#1b9e9e"   # MoBa teal
dark2[1] <- "#2c7fb8"   # ALSPAC blue
dark2[8] <- "#8c2d04"   # TwinLife burgundy
names(dark2) <- cohort_levels

# ---- tighten y-limit ---------------------------------------------------------
y_max <- data %>%
  group_by(Group) %>%
  summarise(total = sum(N), .groups = "drop") %>%
  summarise(max_total = max(total)) %>%
  pull(max_total)
y_lim <- ceiling(y_max / 1000) * 1000 * 1.04

# ---- convert x to numeric to control gap -------------------------------------
group_pos <- c("Preschool age (1–5)" = 1.00,
               "School age (6–12)"   = 1.70)

data2 <- data %>%
  mutate(
    pos   = unname(group_pos[Group]), 
    Group = factor(Group, levels = names(group_pos))
  )

# ---- final barplot -----------------------------------------------------------
ggplot(data2, aes(x = pos, y = N, fill = Cohort)) +
  geom_col(width = 0.6) +  
  scale_x_continuous(
    breaks = unname(group_pos),
    labels = names(group_pos),
    expand = expansion(mult = c(0.01, 0.01))
  ) +
  scale_y_continuous(
    limits = c(0, y_lim),
    expand = expansion(mult = c(0, 0.01)),
    labels = function(x) paste0(formatC(x/1000, format = "f", digits = 0), "K")
  ) +
  scale_fill_manual(values = dark2, breaks = cohort_levels) +
  labs(
    x = NULL,
    y = "Complete trios"
  ) +
  theme_minimal(base_size = 22) +
  theme(
    text = element_text(face = "bold", colour = "black"),
    axis.title.y = element_text(size = 24, margin = margin(r = 10), face = "plain"),
    axis.text.x  = element_text(size = 20),
    axis.text.y  = element_text(size = 20),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    legend.position = "right",
    legend.title = element_blank(),
    legend.text  = element_text(size = 22),
    legend.key.height = unit(20, "pt"),
    legend.key.width  = unit(28, "pt"),
    plot.margin = margin(t = 10, r = 10, b = 10, l = 10)
  )

###########
## Fig 2 ##
###########

#read in results files age 1-5
data <- readr::read_csv("MVMA_RE_1step_lrm_trio_all_1-5y.CSV")

# ---- prepare data ------------------------------------------------------------
data <- data %>%
  filter(str_detect(PRS_predictor, "^prs_[cmf]_qc_")) %>%
  mutate(
    Role = case_when(
      str_detect(PRS_predictor, "^prs_m_") ~ "Maternal effect",
      str_detect(PRS_predictor, "^prs_f_") ~ "Paternal effect",
      str_detect(PRS_predictor, "^prs_c_") ~ "Direct effect"
    ),
    PRS = PRS_predictor %>%
      str_remove("^prs_[cmf]_qc_") %>%
      str_remove("_trio$") %>%
      str_replace("^psych_", "") %>%
      str_replace_all("pc1", "PC1") %>%
      str_replace_all("mean", "MEAN") %>%
      str_replace_all("EDU", "EA") %>%
      toupper(),
    Sig = if_else(FDR_p < 0.05, "FDR significant", "Not FDR significant"),
    
    PRS_group = case_when(
      PRS %in% c("CDG", "MEAN", "PC1")                            ~ "General psych",
      PRS %in% c("ANX", "MDD", "NEUROT", "PPD", "PTSD")           ~ "Affective",
      PRS %in% c("AN","OCD")                                      ~ "Compulsive",
      PRS %in% c("ADHD", "ASD")                                   ~ "Neuro",
      PRS %in% c("BD", "SCZ")                                     ~ "Psychotic",
      PRS %in% c("ALC", "CIG")                                    ~ "Substance",
      PRS %in% c("EA", "INSOM")                                   ~ "Other"
    )
  )

# ---- define PRS levels -------------------------------------------------------
prs_levels <- c(
  "CDG", "MEAN", "PC1",                            # General psych
  "ANX", "MDD", "NEUROT", "PPD", "PTSD",           # Affective
  "AN", "OCD",                                     # Compulsive
  "ADHD", "ASD",                                   # Neurodev
  "BD", "SCZ",                                     # Psychotic
  "ALC", "CIG",                                    # Substance
  "EA", "INSOM"                                    # Other
)

#reversing for top-to-bottom plotting
data$PRS <- factor(data$PRS, levels = rev(prs_levels))

#ordered PRS groups (for facet order)
data$PRS_group <- factor(data$PRS_group,
                             levels = c("General psych", "Affective", "Compulsive", "Neuro",  "Psychotic", "Substance", "Other")
)

#reorder role levels 
data$Role <- factor(data$Role, levels = c("Direct effect", "Paternal effect", "Maternal effect"))

# ---- label outcomes ----------------------------------------------------------
data$CBCL_outcome <- factor(data$CBCL_outcome, levels = c(
  "h_c_total_1_5y", "h_c_extern_1_5y", "h_c_intern_1_5y"
), labels = c(
  "Total problems", "Externalising", "Internalising"
))

#create outline var
data$Outline <- ifelse(data$Sig == "FDR significant", "black", as.character(data$Role))

# ---- colour palette ----------------------------------------------------------
role_colours <- c(
  "Maternal effect" = "#AD002AFF",     
  "Paternal effect" = "#00468BFF",    
  "Direct effect"  = "#0099B4FF"      
)

# ---- plot total difficulties age 1-5 -----------------------------------------
tot <- ggplot(filter(data, CBCL_outcome == "Total problems"), 
                aes(x = Beta, y = PRS, shape = Role)) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey50") +
  geom_errorbarh(
    aes(xmin = CI_Lower, xmax = CI_Upper, alpha = Sig, colour = Role, group = Role),
    height = 0,
    linewidth = 1.06,
    position = position_dodge(width = 0.8)
  ) +
  geom_point(
    aes(alpha = Sig, fill = Role, colour = Outline, group=Role),
    size = 2.45,
    stroke = 0.8,
    position = position_dodge(width = 0.8)
  ) +
  scale_alpha_manual(values = c("FDR significant" = 1, "Not FDR significant" = 0.3)) +
  scale_shape_manual(values = c("Maternal effect" = 24, "Paternal effect" = 22, "Direct effect" = 21)) +
  scale_fill_manual(values = role_colours) +
  scale_colour_manual(values = c(role_colours, "black" = "black"), guide = "none") +
  coord_cartesian(xlim = c(-0.12, 0.12), clip = "off") +
  guides(
    shape = guide_legend(order = 1, reverse = TRUE),
    fill = guide_legend(order = 1, reverse = TRUE),
    alpha = guide_legend(order = 2)
  ) +
  xlab(expression(beta ~ "total difficulties")) +
  facet_grid(
    rows = vars(PRS_group),
    scales = "free_y",
    space = "free"
  ) +
  theme_minimal(base_size = 18) +
  theme(
    panel.grid = element_blank(),
    panel.spacing.y = unit(0.2, "lines"), 
    panel.border = element_rect(color = "grey60", fill = NA, linewidth = 0.6),
    axis.title.y = element_blank(),
    legend.title = element_blank(),
    legend.position = "top",
    legend.justification = "center",
    legend.text = element_text(size = 14),
    legend.key.size = unit(1.1, "lines"),
    strip.text = element_blank(),
    strip.placement = "outside",
    plot.margin = margin(0, 0, 0, 0)
  )

# ---- plot internalising age 1-5 ----------------------------------------------
int <- ggplot(filter(data, CBCL_outcome == "Internalising"), 
              aes(x = Beta, y = PRS, shape = Role)) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey50") +
  geom_errorbarh(
    aes(xmin = CI_Lower, xmax = CI_Upper, alpha = Sig, colour = Role, group = Role),
    height = 0,
    linewidth = 1.06,
    position = position_dodge(width = 0.8)
  ) +
  geom_point(
    aes(alpha = Sig, fill = Role, colour = Outline, group=Role),
    size = 2.45,
    stroke = 0.8,
    position = position_dodge(width = 0.8)
  ) +
  scale_alpha_manual(values = c("FDR significant" = 1, "Not FDR significant" = 0.3)) +
  scale_shape_manual(values = c("Maternal effect" = 24, "Paternal effect" = 22, "Direct effect" = 21)) +
  scale_fill_manual(values = role_colours) +
  scale_colour_manual(values = c(role_colours, "black" = "black"), guide = "none") +
  coord_cartesian(xlim = c(-0.12, 0.12), clip = "off") +
  guides(
    shape = guide_legend(order = 1, reverse = TRUE),
    fill = guide_legend(order = 1, reverse = TRUE),
    alpha = guide_legend(order = 2)
  ) +
  xlab(expression(beta ~ "internalising")) +
  facet_grid(
    rows = vars(PRS_group),
    scales = "free_y",
    space = "free"
  ) +
  theme_minimal(base_size = 18) +
  theme(
    panel.grid = element_blank(),
    panel.spacing.y = unit(0.2, "lines"), 
    panel.border = element_rect(color = "grey60", fill = NA, linewidth = 0.6),
    axis.title.y = element_blank(),
    axis.text.y = element_blank(),
    legend.title = element_blank(),
    legend.position = "top",
    legend.justification = "center",
    legend.text = element_text(size = 14),
    legend.key.size = unit(1.1, "lines"),
    strip.text = element_text(face = "bold", size = 10),
    strip.placement = "outside",
    plot.margin = margin(0, 0, 0, 0)
  )

# ---- plot for externalising age 1-5 ------------------------------------------
ext <- ggplot(filter(data, CBCL_outcome == "Externalising"), 
              aes(x = Beta, y = PRS, shape = Role)) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey50") +
  geom_errorbarh(
    aes(xmin = CI_Lower, xmax = CI_Upper, alpha = Sig, colour = Role, group = Role),
    height = 0,
    linewidth = 1.06,
    position = position_dodge(width = 0.8)
  ) +
  geom_point(
    aes(alpha = Sig, fill = Role, colour = Outline, group=Role),
    size = 2.45,
    stroke = 0.8,
    position = position_dodge(width = 0.8)
  ) +
  scale_alpha_manual(values = c("FDR significant" = 1, "Not FDR significant" = 0.3)) +
  scale_shape_manual(values = c("Maternal effect" = 24, "Paternal effect" = 22, "Direct effect" = 21)) +
  scale_fill_manual(values = role_colours) +
  scale_colour_manual(values = c(role_colours, "black" = "black"), guide = "none") +
  coord_cartesian(xlim = c(-0.12, 0.12), clip = "off") +
  guides(
    shape = guide_legend(order = 1, reverse = TRUE),
    fill = guide_legend(order = 1, reverse = TRUE),
    alpha = guide_legend(order = 2)
  ) +
  xlab(expression(beta ~ "externalising")) +
  facet_grid(
    rows = vars(PRS_group),
    scales = "free_y",
    space = "free",
  ) +
  theme_minimal(base_size = 18) +
  theme(
    panel.grid = element_blank(),
    panel.spacing.y = unit(0.2, "lines"), 
    panel.border = element_rect(color = "grey60", fill = NA, linewidth = 0.6),
    axis.title.y = element_blank(),
    axis.text.y = element_blank(),
    legend.title = element_blank(),
    legend.position = "top",
    legend.justification = "center",
    legend.text = element_text(size = 14),
    legend.key.size = unit(1.1, "lines"),
    strip.text = element_blank(),
    strip.placement = "outside",
    plot.margin = margin(0, 0, 0, 0)
  )

# ---- stick plots together/save out -------------------------------------------
combined <- (tot + ext + int) +
  plot_layout(guides = "collect") &
  theme(legend.position = "top")

#print/save out
print(combined)
ggsave("combined_Fig2.png", combined, width = 10.5, height = 8.4, dpi = 300)

###########
## Fig 3 ##
###########

#read in results files age 6-18
data2 <- readr::read_csv("MVMA_RE_1step_lrm_trio_all_6-18y.CSV")

# ---- prepare data ------------------------------------------------------------
data2 <- data2 %>%
  filter(str_detect(PRS_predictor, "^prs_[cmf]_qc_")) %>%
  mutate(
    Role = case_when(
      str_detect(PRS_predictor, "^prs_m_") ~ "Maternal effect",
      str_detect(PRS_predictor, "^prs_f_") ~ "Paternal effect",
      str_detect(PRS_predictor, "^prs_c_") ~ "Direct effect"
    ),
    PRS = PRS_predictor %>%
      str_remove("^prs_[cmf]_qc_") %>%
      str_remove("_trio$") %>%
      str_replace("^psych_", "") %>%
      str_replace_all("pc1", "PC1") %>%
      str_replace_all("mean", "MEAN") %>%
      str_replace_all("EDU", "EA") %>%
      toupper(),
    Sig = if_else(FDR_p < 0.05, "FDR significant", "Not FDR significant"),
    
    PRS_group = case_when(
      PRS %in% c("CDG", "MEAN", "PC1")                            ~ "General psych",
      PRS %in% c("ANX", "MDD", "NEUROT", "PPD", "PTSD")           ~ "Affective",
      PRS %in% c("AN","OCD")                                      ~ "Compulsive",
      PRS %in% c("ADHD", "ASD")                                   ~ "Neuro",
      PRS %in% c("BD", "SCZ")                                     ~ "Psychotic",
      PRS %in% c("ALC", "CIG")                                    ~ "Substance",
      PRS %in% c("EA", "INSOM")                                   ~ "Other"
    )
  )

# ---- define PRS levels -------------------------------------------------------
prs_levels <- c(
  "CDG", "MEAN", "PC1",                            # General psych
  "ANX", "MDD", "NEUROT", "PPD", "PTSD",           # Affective
  "AN", "OCD",                                     # Compulsive
  "ADHD", "ASD",                                   # Neurodev
  "BD", "SCZ",                                     # Psychotic
  "ALC", "CIG",                                    # Substance
  "EA", "INSOM"                                    # Other
)

#reversing for top-to-bottom plotting
data2$PRS <- factor(data2$PRS, levels = rev(prs_levels))

#ordered PRS groups (for facet order)
data2$PRS_group <- factor(data2$PRS_group,
                          levels = c("General psych", "Affective", "Compulsive", "Neuro", "Psychotic", "Substance", "Other")
)

#reorder role
data2$Role <- factor(data2$Role, levels = c("Direct effect", "Paternal effect", "Maternal effect"))

# ---- label outcomes ----------------------------------------------------------
data2$CBCL_outcome <- factor(data2$CBCL_outcome, levels = c(
  "h_c_total_6_18y", "h_c_extern_6_18y", "h_c_intern_6_18y"
), labels = c(
  "Total problems", "Externalising", "Internalising"
))

#create outline var
data2$Outline <- ifelse(data2$Sig == "FDR significant", "black", as.character(data2$Role))

# ---- colour palette ----------------------------------------------------------
role_colours <- c(
  "Maternal effect" = "#AD002AFF",     
  "Paternal effect" = "#00468BFF",    
  "Direct effect"  = "#0099B4FF"      
)

# ---- plot total difficulties age 6-18 ----------------------------------------
tot2 <- ggplot(filter(data2, CBCL_outcome == "Total problems"), 
              aes(x = Beta, y = PRS, shape = Role)) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey50") +
  geom_errorbarh(
    aes(xmin = CI_Lower, xmax = CI_Upper, alpha = Sig, colour = Role, group = Role),
    height = 0,
    linewidth = 1.06,
    position = position_dodge(width = 0.8)
  ) +
  geom_point(
    aes(alpha = Sig, fill = Role, colour = Outline, group=Role),
    size = 2.45,
    stroke = 0.8,
    position = position_dodge(width = 0.8)
  ) +
  scale_alpha_manual(values = c("FDR significant" = 1, "Not FDR significant" = 0.3)) +
  scale_shape_manual(values = c("Maternal effect" = 24, "Paternal effect" = 22, "Direct effect" = 21)) +
  scale_fill_manual(values = role_colours) +
  scale_colour_manual(values = c(role_colours, "black" = "black"), guide = "none") +
  coord_cartesian(xlim = c(-0.12, 0.12), clip = "off") +
  guides(
    shape = guide_legend(order = 1, reverse = TRUE),
    fill = guide_legend(order = 1, reverse = TRUE),
    alpha = guide_legend(order = 2)
  ) +
  xlab(expression(beta ~ "total difficulties")) +
  facet_grid(
    rows = vars(PRS_group),
    scales = "free_y",
    space = "free"
  ) +
  theme_minimal(base_size = 18) +
  theme(
    panel.grid = element_blank(),
    panel.spacing.y = unit(0.2, "lines"), 
    panel.border = element_rect(color = "grey60", fill = NA, linewidth = 0.6),
    axis.title.y = element_blank(),
    legend.title = element_blank(),
    legend.position = "top",
    legend.justification = "center",
    legend.text = element_text(size = 14),
    legend.key.size = unit(1.1, "lines"),
    strip.text = element_blank(),
    strip.placement = "outside",
    plot.margin = margin(0, 0, 0, 0)
  )

# ---- plot internalising age 6-18 ---------------------------------------------
int2 <- ggplot(filter(data2, CBCL_outcome == "Internalising"), 
              aes(x = Beta, y = PRS, shape = Role)) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey50") +
  geom_errorbarh(
    aes(xmin = CI_Lower, xmax = CI_Upper, alpha = Sig, colour = Role, group = Role),
    height = 0,
    linewidth = 1.06,
    position = position_dodge(width = 0.8)
  ) +
  geom_point(
    aes(alpha = Sig, fill = Role, colour = Outline, group=Role),
    size = 2.45,
    stroke = 0.8,
    position = position_dodge(width = 0.8)
  ) +
  scale_alpha_manual(values = c("FDR significant" = 1, "Not FDR significant" = 0.3)) +
  scale_shape_manual(values = c("Maternal effect" = 24, "Paternal effect" = 22, "Direct effect" = 21)) +
  scale_fill_manual(values = role_colours) +
  scale_colour_manual(values = c(role_colours, "black" = "black"), guide = "none") +
  coord_cartesian(xlim = c(-0.12, 0.12), clip = "off") +
  guides(
    shape = guide_legend(order = 1, reverse = TRUE),
    fill = guide_legend(order = 1, reverse = TRUE),
    alpha = guide_legend(order = 2)
  ) +
  xlab(expression(beta ~ "internalising")) +
  facet_grid(
    rows = vars(PRS_group),
    scales = "free_y",
    space = "free"
  ) +
  theme_minimal(base_size = 18) +
  theme(
    panel.grid = element_blank(),
    panel.spacing.y = unit(0.2, "lines"), 
    panel.border = element_rect(color = "grey60", fill = NA, linewidth = 0.6),
    axis.title.y = element_blank(),
    axis.text.y = element_blank(),
    legend.title = element_blank(),
    legend.position = "top",
    legend.justification = "center",
    legend.text = element_text(size = 14),
    legend.key.size = unit(1.1, "lines"),
    strip.text = element_text(face = "bold", size = 10),
    strip.placement = "outside",
    plot.margin = margin(0, 0, 0, 0)
  )

# ---- plot externalising age 6-18 ---------------------------------------------
ext2 <- ggplot(filter(data2, CBCL_outcome == "Externalising"), 
              aes(x = Beta, y = PRS, shape = Role)) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey50") +
  geom_errorbarh(
    aes(xmin = CI_Lower, xmax = pmin(CI_Upper,  0.132), alpha = Sig, colour = Role, group = Role),
    height = 0,
    linewidth = 1.06,
    position = position_dodge(width = 0.8)
  ) +
  geom_point(
    aes(alpha = Sig, fill = Role, colour = Outline, group=Role),
    size = 2.45,
    stroke = 0.8,
    position = position_dodge(width = 0.8)
  ) +
  scale_alpha_manual(values = c("FDR significant" = 1, "Not FDR significant" = 0.3)) +
  scale_shape_manual(values = c("Maternal effect" = 24, "Paternal effect" = 22, "Direct effect" = 21)) +
  scale_fill_manual(values = role_colours) +
  scale_colour_manual(values = c(role_colours, "black" = "black"), guide = "none") +
  coord_cartesian(xlim = c(-0.12, 0.12), clip = "off") +
  guides(
    shape = guide_legend(order = 1, reverse = TRUE),
    fill = guide_legend(order = 1, reverse = TRUE),
    alpha = guide_legend(order = 2)
  ) +
  xlab(expression(beta ~ "externalising")) +
  facet_grid(
    rows = vars(PRS_group),
    scales = "free_y",
    space = "free"
  ) +
  theme_minimal(base_size = 18) +
  theme(
    panel.grid = element_blank(),
    panel.spacing.y = unit(0.2, "lines"), 
    panel.border = element_rect(color = "grey60", fill = NA, linewidth = 0.6),
    axis.title.y = element_blank(),
    axis.text.y = element_blank(),
    legend.title = element_blank(),
    legend.position = "top",
    legend.justification = "center",
    legend.text = element_text(size = 14),
    legend.key.size = unit(1.1, "lines"),
    strip.text = element_blank(),
    strip.placement = "outside",
    plot.margin = margin(0, 0, 0, 0)
  )

# ---- stick plots together/save out -------------------------------------------
combined2 <- (tot2 + ext2 + int2) +
  plot_layout(guides = "collect") &
  theme(legend.position = "top")

#print/save out
print(combined2)
ggsave("combined_Fig3.png", combined2, width = 10.5, height = 8.4, dpi = 300)

############
## Fig S1 ##
############

# ---- paths -------------------------------------------------------------------
dir_in  <- "./Cohort-level_PRScmf_correlations/" 
out_png <- file.path(dir_in, "heatmap_mother_father_PGS_18x18.png")
out_csv <- file.path(dir_in, "matrix_mother_father_PGS_18x18.csv")

# ---- helpers -----------------------------------------------------------------

#extract age / reporter / cohort from filename:
parse_fname <- function(x) {
  m <- stringr::str_match(basename(x), "^PRScmf_corr_([0-9]+y)_(mo|fa|mixed)_([^\\.]+)\\.csv$")
  if (any(is.na(m))) warning("Filename didn't match expected pattern: ", basename(x))
  tibble(file = x,
         age = m[,2],
         reporter = m[,3],
         cohort_raw = m[,4])
}

#identify mother / father PRS columns
is_f <- function(x) stringr::str_detect(x, "^prs_f_qc_.*_trio$")
is_m <- function(x) stringr::str_detect(x, "^prs_m_qc_.*_trio$")

#extract trait label from PRS string
trait_from_str <- function(x)  stringr::str_match(x, "prs_[fm]_qc_([A-Za-z0-9_]+)_trio")[,2]

#programmatic cohort grouping: keep only letters, uppercase
cohort_stem <- function(x) x %>% toupper() %>% gsub("[^A-Z]+", "", .)

#display label harmonization
label_map <- c(
  "EDU" = "EA",
  "psych_pc1" = "PC1",
  "psych_mean" = "MEAN"
)

prs_levels <- rev(c("CDG", "MEAN", "PC1",                   # General psych 
                    "ANX", "MDD", "NEUROT", "PPD", "PTSD",  # Affective 
                    "AN", "OCD",                            # Compulsive 
                    "ADHD", "ASD",                          # Neurodev 
                    "BD", "SCZ",                            # Psychotic 
                    "ALC", "CIG",                           # Substance 
                    "EA", "INSOM"                           # Other
))

# ---- load files --------------------------------------------------------------
files <- list.files(dir_in, pattern = "^PRScmf_corr_.*\\.csv$", full.names = TRUE)

meta <- purrr::map_dfr(files, parse_fname) %>%
  mutate(age_num = suppressWarnings(as.numeric(gsub("[^0-9.]", "", age))))

dat <- purrr::map2_dfr(meta$file, seq_len(nrow(meta)), function(f, i) {
  readr::read_csv(f, show_col_types = FALSE) %>% mutate(.file_row_id = i)
}) %>%
  left_join(meta %>% mutate(.file_row_id = row_number()), by = ".file_row_id") %>%
  select(-.file_row_id)

# ---- mother–father pairs -----------------------------------------------------
mf_pairs_all <- dat %>%
  filter( (is_f(Predictor1) & is_m(Predictor2)) | (is_m(Predictor1) & is_f(Predictor2)) ) %>%
  mutate(
    f_trait_raw = if_else(is_f(Predictor1), trait_from_str(Predictor1), trait_from_str(Predictor2)),
    m_trait_raw = if_else(is_m(Predictor1), trait_from_str(Predictor1), trait_from_str(Predictor2)),
    cohort_group = cohort_stem(cohort_raw),
    w_row = if_else(is.na(N) | N <= 0, 1, N)  # row weight for reporter averaging
  ) %>%
  select(cohort_raw, cohort_group, age, age_num, reporter,
         f_trait_raw, m_trait_raw, r = Correlation, w_row)

# ---- keep only same-trait pairs ----------------------------------------------
mf_same_trait <- mf_pairs_all %>%
  filter(!is.na(f_trait_raw), !is.na(m_trait_raw), f_trait_raw == m_trait_raw) %>%
  rename(trait_raw = f_trait_raw) %>%
  select(cohort_raw, cohort_group, age, age_num, reporter, trait_raw, r, w_row)

# ---- youngest age per cohort_raw ---------------------------------------------
youngest_age_lookup <- mf_same_trait %>%
  distinct(cohort_raw, age_num) %>%
  group_by(cohort_raw) %>%
  summarise(age_min = min(age_num, na.rm = TRUE), .groups = "drop")

mf_youngest <- mf_same_trait %>%
  inner_join(youngest_age_lookup, by = "cohort_raw") %>%
  filter(age_num == age_min) %>%
  select(-age_min)

# ---- average across reporters at the youngest age ----------------------------
within_cohort <- mf_youngest %>%
  group_by(cohort_group, trait_raw) %>%
  summarise(
    r = weighted.mean(r, w_row, na.rm = TRUE),
    N_total = sum(w_row[!is.na(r)]),
    .groups = "drop"
  )

agg <- within_cohort %>%
  group_by(trait_raw) %>%
  summarise(
    r_mean = weighted.mean(r, N_total, na.rm = TRUE),
    r_sd   = sd(r, na.rm = TRUE),
    k      = n(),
    r_se   = r_sd / sqrt(k),
    r_lo   = r_mean - 1.96 * r_se,
    r_hi   = r_mean + 1.96 * r_se,
    .groups = "drop"
  ) %>%
  mutate(
    trait = recode(trait_raw, !!!label_map, .default = trait_raw)
  )

# ---- enforce trait order -----------------------------------------------------
levs <- c(prs_levels, setdiff(unique(agg$trait), prs_levels))

agg_ord <- agg %>%
  mutate(trait = factor(trait, levels = levs)) %>%
  arrange(trait)

# ---- save tidy estimates -----------------------------------------------------
agg_ord %>%
  select(trait, k, r_mean, r_lo, r_hi, r_sd) %>%
  write_csv(out_csv)

# ---- plot dot + CI plot ------------------------------------------------------
p <- ggplot(agg_ord, aes(x = r_mean, y = trait)) +
    geom_vline(xintercept = 0, colour = "grey30", linewidth = 0.5, linetype = "dashed") +
    geom_point(size = 3.9, colour = "#2A9D8F") +
    geom_errorbarh(aes(xmin = r_lo, xmax = r_hi), alpha = 0.65, height = 0, linewidth = 1, colour = "#2A9D8F") +
    scale_x_continuous(limits = c(-0.05, max(0.15, max(agg_ord$r_hi, na.rm = TRUE))),
                       breaks = scales::pretty_breaks()) +
    labs(
      x = "Mean mother–father PGS correlation (r)",
      y = NULL,
    ) +
    theme_minimal(base_size = 18) +
    theme(panel.grid.major.y = element_blank(),
          plot.title = element_text(face = "bold"))

#print/save out
print(p)
ggsave("mother_father_PGS_correlations.png", p, width = 6.6, height = 6.6, dpi = 300)

