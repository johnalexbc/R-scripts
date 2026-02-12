library(readxl)
library(dplyr)
library(janitor)
library(meta)
library(stringr)

#### 0. Path and data import --------------------------------------------

file_path <- "C:/Users/johna/Desktop/Research/Meta-análisis Rosario/Hoja de extracción meta análisis rosario.xlsx"

dat_all <- read_xlsx(file_path) |>
  clean_names()

# Output folder = same as Excel file
out_dir <- dirname(file_path)


#### 0b. Diagnostic report: why rows are excluded -----------------------

dat_diagnostic <- dat_all |>
  mutate(
    excl_single_arm = design == "single_arm",
    excl_no_ctrl_n  = is.na(n_ctrl) | n_ctrl <= 0,
    excl_no_int_n   = is.na(n_int)  | n_int  <= 0,
    excl_not_post   = follow_up_momento != "post" | is.na(follow_up_momento),
    exclusion_reason = case_when(
      excl_single_arm ~ "single_arm_design",
      excl_no_ctrl_n  ~ "missing_or_zero_n_ctrl",
      excl_no_int_n   ~ "missing_or_zero_n_int",
      excl_not_post   ~ "not_post_follow_up",
      TRUE            ~ "included"
    )
  )

diag_file <- file.path(out_dir, "diagnostic_exclusion_report.csv")
write.csv(dat_diagnostic, diag_file, row.names = FALSE)
message("Saved diagnostic exclusion report: ", diag_file)


#### 1. Filter: only controlled studies, POST measurements --------------

dat_control_post <- dat_all |>
  filter(
    design != "single_arm",
    !is.na(n_ctrl), n_ctrl > 0,
    !is.na(n_int),  n_int  > 0,
    follow_up_momento == "post"
  )


#### 2. Harmonise direction: higher = WORSE -----------------------------

dat_control_post <- dat_control_post |>
  mutate(
    mean_ctrl_harm = if_else(direction_outcome == "higher=better",
                             -mean_ctrl, mean_ctrl),
    mean_int_harm  = if_else(direction_outcome == "higher=better",
                             -mean_int,  mean_int)
  )


#### 3. Map fine outcomes to 11 outcome groups --------------------------

dat_control_post <- dat_control_post |>
  mutate(
    outcome_group = case_when(
      # 1) Caregiver burden / distress
      outcome_domain %in% c("burden","malestar_cuidador") ~
        "caregiver_burden_distress",
      
      # 2) Caregiver emotional symptoms / mental health
      outcome_domain %in% c("depresion","ansiedad",
                            "ansiedad_depresion","malestar_psicologico",
                            "afecto_positivo",
                            "sintomas_neuropsiquiatricos_cuidador",
                            "salud_mental") ~
        "caregiver_emotional_symptoms_mental_health",
      
      # 3) Caregiver quality of life (including sleep)
      outcome_domain %in% c("calidad_vida","calidad_sueno") ~
        "caregiver_quality_of_life",
      
      # 4) Patient neuropsychiatric symptoms / behaviour
      outcome_domain %in% c("sintomas_neuropsiquiatricos",
                            "problemas_memoria_conducta") ~
        "patient_neuropsychiatric_symptoms_behavior",
      
      # 5) Caregiver stress
      outcome_domain %in% c("estres_percibido","estres_reaccion") ~
        "caregiver_stress",
      
      # 6) Caregiver physical / general health
      outcome_domain %in% c("salud_fisica","salud_general") ~
        "caregiver_physical_general_health",
      
      # 7) Caregiver self-efficacy / competence
      outcome_domain %in% c("autoeficacia","autoeficacia_apoyo",
                            "autoeficacia_bpsd",
                            "autoeficacia_cuidado_rutina",
                            "autoeficacia_malestar",
                            "competencia_cuidador","dominio") ~
        "caregiver_self_efficacy_competence",
      
      # 8) Caregiver coping and resources
      outcome_domain %in% c("afrontamiento_apoyo",
                            "afrontamiento_problema",
                            "afrontamiento_resignado",
                            "apoyo_social",
                            "resiliencia") ~
        "caregiver_coping_resources",
      
      # 9) Dementia knowledge and attitudes
      outcome_domain %in% c("conocimiento",
                            "conocimiento_AD",
                            "actitudes_demencia") ~
        "dementia_knowledge_attitudes",
      
      # 10) Patient daily functioning (ADL / IADL)
      outcome_domain %in% c("ADL",
                            "funcionamiento_diario",
                            "funcionamiento_instrumental") ~
        "patient_daily_functioning",
      
      # 11) Patient cognitive function
      outcome_domain %in% c("tamizaje_cognitivo") ~
        "patient_cognitive_function",
      
      TRUE ~ NA_character_
    )
  ) |>
  filter(!is.na(outcome_group))


#### 4. Study label: "Xxx et al., 2012" ---------------------------------

make_studlab <- function(author_core, year) {
  auth <- ifelse(nchar(author_core) == 0,
                 "Study",
                 str_to_title(author_core))
  paste0(auth, " et al., ", year)
}


#### 5. Collapsed data frame (1 row per study × outcome_group) ----------

dat_for_meta <- dat_control_post |>
  mutate(
    # first block of letters from 'd' = first author
    author_core = gsub("[^a-z]+.*$", "", tolower(d)),
    studlab     = make_studlab(author_core, year)
  ) |>
  group_by(outcome_group, studlab) |>
  summarise(
    n_ctrl = first(n_ctrl),
    n_int  = first(n_int),
    mean_ctrl_harm = mean(mean_ctrl_harm, na.rm = TRUE),
    mean_int_harm  = mean(mean_int_harm,  na.rm = TRUE),
    sd_ctrl = mean(sd_ctrl, na.rm = TRUE),
    sd_int  = mean(sd_int,  na.rm = TRUE),
    .groups = "drop"
  )


#### 6. Meta-analysis per outcome group ---------------------------------

run_meta_group <- function(data, group, sm = "SMD") {
  
  dat <- data |>
    filter(outcome_group == group) |>
    filter(
      !is.na(mean_int_harm), !is.na(mean_ctrl_harm),
      !is.na(sd_int), !is.na(sd_ctrl),
      n_int  > 0, n_ctrl > 0
    )
  
  if (nrow(dat) == 0L) {
    stop(paste0("No rows for outcome_group = '", group,
                "' in post measurements with control group."))
  }
  
  m <- metacont(
    n.e    = dat$n_ctrl,
    mean.e = dat$mean_ctrl_harm,
    sd.e   = dat$sd_ctrl,
    n.c    = dat$n_int,
    mean.c = dat$mean_int_harm,
    sd.c   = dat$sd_int,
    studlab = dat$studlab,
    sm      = sm,
    common  = FALSE,
    random  = TRUE,
    method.tau       = "REML",
    method.random.ci = "classic",
    prediction       = FALSE
  )
  
  return(m)
}


#### 7. List of meta objects by outcome_group ---------------------------

outcome_groups <- sort(unique(dat_for_meta$outcome_group))

meta_list <- lapply(outcome_groups, function(g) {
  run_meta_group(dat_for_meta, group = g)
})
names(meta_list) <- outcome_groups


#### 8. Nice English labels for outcomes -------------------------------

label_group_en <- function(group) {
  dplyr::case_when(
    group == "caregiver_burden_distress" ~
      "Caregiver burden / distress",
    group == "caregiver_emotional_symptoms_mental_health" ~
      "Caregiver emotional symptoms / mental health",
    group == "caregiver_quality_of_life" ~
      "Caregiver quality of life (including sleep)",
    group == "patient_neuropsychiatric_symptoms_behavior" ~
      "Patient neuropsychiatric symptoms / behaviour",
    group == "caregiver_stress" ~
      "Caregiver stress",
    group == "caregiver_physical_general_health" ~
      "Caregiver physical / general health",
    group == "caregiver_self_efficacy_competence" ~
      "Caregiver self-efficacy / competence",
    group == "caregiver_coping_resources" ~
      "Caregiver coping and resources",
    group == "dementia_knowledge_attitudes" ~
      "Dementia knowledge and attitudes",
    group == "patient_daily_functioning" ~
      "Patient daily functioning (ADL / IADL)",
    group == "patient_cognitive_function" ~
      "Patient cognitive function",
    TRUE ~ group
  )
}

sanitize_name <- function(x) {
  x %>%
    tolower() %>%
    gsub("[^a-z0-9]+", "_", .) %>%
    gsub("_+", "_", .) %>%
    gsub("^_|_$", "", .)
}


#### 9. Forest + funnel -------------------------------------------------

save_plots_for_group <- function(m, group, out_dir) {
  
  group_clean <- sanitize_name(group)
  group_label <- label_group_en(group)   # plot title
  
  k <- length(m$TE)
  # algo más alto por estudio para que haya aire
  height_in <- max(5.5, 3.0 + 0.30 * k)
  
  ## FOREST -------------------------------------------------------------
  
  forest_file <- file.path(out_dir, paste0("forest_", group_clean, ".tiff"))
  
  message("Saving forest for ", group_label,
          " (", group, "), k = ", k,
          " → ", forest_file)
  
  tiff(
    filename    = forest_file,
    width       = 8.5,
    height      = height_in,
    units       = "in",
    res         = 600,
    compression = "lzw"
  )
  
  par(mar = c(5.5, 4.5, 4, 2))
  
  forest(
    m,
    main        = group_label,
    xlab        = "Standardized mean difference (Hedges g)",
    smlab       = "",
    prediction  = FALSE,  # sin intervalo extra: solo diamante (IC95%)
    common      = FALSE,
    random      = TRUE,
    leftcols    = c("studlab"),
    leftlabs    = c("Study"),
    
    ## --- SOLO I2, tau2 y p de heterogeneidad ---
    print.I2       = TRUE,
    print.I2.ci    = FALSE,
    print.tau2     = TRUE,
    print.tau2.ci  = FALSE,
    print.Q        = FALSE,
    print.pval.Q   = TRUE,
    
    ## tamaño general
    fontsize       = 8,
    ## fuente de la línea de heterogeneidad
    fs.hetstat     = 4,
    
    ## más espacio entre columnas y el bosque
    colgap.forest.left  = "0.10cm",
    colgap.forest.right = "0.60cm"
    # sin "favours" text
  )
  
  dev.off()
  
  ## FUNNEL -------------------------------------------------------------
  
  funnel_file <- file.path(out_dir, paste0("funnel_", group_clean, ".tiff"))
  
  message("Saving funnel for ", group_label,
          " (", group, "), k = ", k,
          " → ", funnel_file)
  
  tiff(
    filename    = funnel_file,
    width       = 7,
    height      = 7,
    units       = "in",
    res         = 600,
    compression = "lzw"
  )
  
  par(mar = c(5, 4.5, 4, 2))
  
  funnel(
    m,
    main    = group_label,
    xlab    = "Standardized mean difference (Hedges g)",
    ylab    = "Standard error",
    studlab = FALSE
  )
  
  dev.off()
}


#### 10. Loop: save all plots -------------------------------------------

for (g in outcome_groups) {
  m <- meta_list[[g]]
  save_plots_for_group(m, g, out_dir)
}

#### 11. Global forest: prioritized outcome (1 per study) -------------

# Column that contains the instrument / scale name
instrument_col <- "outcome_scale"

# Safety check
if (!instrument_col %in% names(dat_control_post)) {
  stop("The column specified in instrument_col ('", instrument_col,
       "') is not in dat_control_post. Run names(dat_control_post) to see available columns.")
}

# Priority requested for the global plot: 1 outcome per study
outcome_priority <- c(
  "caregiver_burden_distress",
  "caregiver_emotional_symptoms_mental_health",
  "caregiver_quality_of_life",
  "patient_neuropsychiatric_symptoms_behavior",
  "caregiver_self_efficacy_competence"
)

# Region mapping for study-label shading
region_to_group <- function(region_value) {
  reg <- tolower(ifelse(is.na(region_value), "", region_value))
  dplyr::case_when(
    str_detect(reg, "asia") ~ "asia",
    str_detect(reg, "latin") | str_detect(reg, "latam") ~ "latin_america",
    str_detect(reg, "africa") ~ "africa",
    TRUE ~ "other"
  )
}

region_color_map <- c(
  asia = "#1B9E77",
  latin_america = "#7570B3",
  africa = "#D95F02",
  other = "#333333"
)

# 1) Keep exactly one prioritized outcome row per study-year
dat_all_effects <- dat_control_post |>
  filter(
    !is.na(mean_int_harm), !is.na(mean_ctrl_harm),
    !is.na(sd_int), !is.na(sd_ctrl),
    n_int  > 0, n_ctrl > 0,
    !is.na(outcome_group),
    outcome_group %in% outcome_priority
  ) |>
  mutate(
    study_id          = paste0(tolower(d), "_", year),
    outcome_rank      = match(outcome_group, outcome_priority),
    author_core      = gsub("[^a-z]+.*$", "", tolower(d)),
    author_label     = str_to_title(author_core),
    instrument_label = .data[[instrument_col]],
    instrument_label = if_else(is.na(instrument_label) | instrument_label == "",
                               "Unknown scale", instrument_label),
    region_group     = region_to_group(region),
    studlab_instr    = paste0(author_label, " et al., ", year,
                              " (", instrument_label, ")")
  ) |>
  arrange(study_id, outcome_rank, instrument_label) |>
  group_by(study_id) |>
  slice(1) |>
  ungroup()

# 2) Meta-analysis using one prioritized row per study
m_all <- metacont(
  n.e    = dat_all_effects$n_ctrl,
  mean.e = dat_all_effects$mean_ctrl_harm,
  sd.e   = dat_all_effects$sd_ctrl,
  n.c    = dat_all_effects$n_int,
  mean.c = dat_all_effects$mean_int_harm,
  sd.c   = dat_all_effects$sd_int,
  studlab = dat_all_effects$studlab_instr,
  sm      = "SMD",
  common  = FALSE,
  random  = TRUE,
  method.tau       = "REML",
  method.random.ci = "classic",
  prediction       = FALSE
)

study_label_colors <- unname(region_color_map[dat_all_effects$region_group])

# 3) Plot height based on k
k_all      <- length(m_all$TE)
height_all <- max(6, 3.0 + 0.30 * k_all)

forest_all_file <- file.path(out_dir, "forest_all_outcomes_instruments.tiff")

message("Saving global forest with prioritized outcomes (1 per study) → ", forest_all_file)

tiff(
  filename    = forest_all_file,
  width       = 8.5,
  height      = height_all,
  units       = "in",
  res         = 600,
  compression = "lzw"
)

par(mar = c(5.5, 4.5, 4, 2))

forest(
  m_all,
  main        = "Prioritized outcomes by study (1 outcome per study)",
  xlab        = "Standardized mean difference (Hedges g)",
  smlab       = "",
  leftcols    = c("studlab", "w.random"),
  leftlabs    = c("Study (instrument)", "Weight"),
  overall     = TRUE,
  prediction  = FALSE,
  common      = FALSE,
  random      = TRUE,
  text.random = "Random effects model",
  col.study   = study_label_colors,
  
  print.I2       = TRUE,
  print.I2.ci    = FALSE,
  print.tau2     = TRUE,
  print.tau2.ci  = FALSE,
  print.Q        = FALSE,
  print.pval.Q   = TRUE,
  
  fontsize       = 8,
  fs.hetstat     = 4,
  colgap.forest.left  = "0.10cm",
  colgap.forest.right = "0.60cm"
)

dev.off()
library(readxl)
library(dplyr)
library(janitor)
library(meta)
library(stringr)

#### 0. Path and data import --------------------------------------------

file_path <- "C:/Users/johna/Desktop/Research/Meta-análisis Rosario/Hoja de extracción meta análisis rosario.xlsx"

dat_all <- read_xlsx(file_path) |>
  clean_names()

# Output folder = same as Excel file
out_dir <- dirname(file_path)


#### 0b. Diagnostic report: why rows are excluded -----------------------

dat_diagnostic <- dat_all |>
  mutate(
    excl_single_arm = design == "single_arm",
    excl_no_ctrl_n  = is.na(n_ctrl) | n_ctrl <= 0,
    excl_no_int_n   = is.na(n_int)  | n_int  <= 0,
    excl_not_post   = follow_up_momento != "post" | is.na(follow_up_momento),
    exclusion_reason = case_when(
      excl_single_arm ~ "single_arm_design",
      excl_no_ctrl_n  ~ "missing_or_zero_n_ctrl",
      excl_no_int_n   ~ "missing_or_zero_n_int",
      excl_not_post   ~ "not_post_follow_up",
      TRUE            ~ "included"
    )
  )

diag_file <- file.path(out_dir, "diagnostic_exclusion_report.csv")
write.csv(dat_diagnostic, diag_file, row.names = FALSE)
message("Saved diagnostic exclusion report: ", diag_file)


#### 1. Filter: only controlled studies, POST measurements --------------

dat_control_post <- dat_all |>
  filter(
    design != "single_arm",
    !is.na(n_ctrl), n_ctrl > 0,
    !is.na(n_int),  n_int  > 0,
    follow_up_momento == "post"
  )


#### 2. Harmonise direction: higher = WORSE -----------------------------

dat_control_post <- dat_control_post |>
  mutate(
    mean_ctrl_harm = if_else(direction_outcome == "higher=better",
                             -mean_ctrl, mean_ctrl),
    mean_int_harm  = if_else(direction_outcome == "higher=better",
                             -mean_int,  mean_int)
  )


#### 3. Map fine outcomes to 11 outcome groups --------------------------

dat_control_post <- dat_control_post |>
  mutate(
    outcome_group = case_when(
      # 1) Caregiver burden / distress
      outcome_domain %in% c("burden","malestar_cuidador") ~
        "caregiver_burden_distress",
      
      # 2) Caregiver emotional symptoms / mental health
      outcome_domain %in% c("depresion","ansiedad",
                            "ansiedad_depresion","malestar_psicologico",
                            "afecto_positivo",
                            "sintomas_neuropsiquiatricos_cuidador",
                            "salud_mental") ~
        "caregiver_emotional_symptoms_mental_health",
      
      # 3) Caregiver quality of life (including sleep)
      outcome_domain %in% c("calidad_vida","calidad_sueno") ~
        "caregiver_quality_of_life",
      
      # 4) Patient neuropsychiatric symptoms / behaviour
      outcome_domain %in% c("sintomas_neuropsiquiatricos",
                            "problemas_memoria_conducta") ~
        "patient_neuropsychiatric_symptoms_behavior",
      
      # 5) Caregiver stress
      outcome_domain %in% c("estres_percibido","estres_reaccion") ~
        "caregiver_stress",
      
      # 6) Caregiver physical / general health
      outcome_domain %in% c("salud_fisica","salud_general") ~
        "caregiver_physical_general_health",
      
      # 7) Caregiver self-efficacy / competence
      outcome_domain %in% c("autoeficacia","autoeficacia_apoyo",
                            "autoeficacia_bpsd",
                            "autoeficacia_cuidado_rutina",
                            "autoeficacia_malestar",
                            "competencia_cuidador","dominio") ~
        "caregiver_self_efficacy_competence",
      
      # 8) Caregiver coping and resources
      outcome_domain %in% c("afrontamiento_apoyo",
                            "afrontamiento_problema",
                            "afrontamiento_resignado",
                            "apoyo_social",
                            "resiliencia") ~
        "caregiver_coping_resources",
      
      # 9) Dementia knowledge and attitudes
      outcome_domain %in% c("conocimiento",
                            "conocimiento_AD",
                            "actitudes_demencia") ~
        "dementia_knowledge_attitudes",
      
      # 10) Patient daily functioning (ADL / IADL)
      outcome_domain %in% c("ADL",
                            "funcionamiento_diario",
                            "funcionamiento_instrumental") ~
        "patient_daily_functioning",
      
      # 11) Patient cognitive function
      outcome_domain %in% c("tamizaje_cognitivo") ~
        "patient_cognitive_function",
      
      TRUE ~ NA_character_
    )
  ) |>
  filter(!is.na(outcome_group))


#### 4. Study label: "Xxx et al., 2012" ---------------------------------

make_studlab <- function(author_core, year) {
  auth <- ifelse(nchar(author_core) == 0,
                 "Study",
                 str_to_title(author_core))
  paste0(auth, " et al., ", year)
}


#### 5. Collapsed data frame (1 row per study × outcome_group) ----------

dat_for_meta <- dat_control_post |>
  mutate(
    # first block of letters from 'd' = first author
    author_core = gsub("[^a-z]+.*$", "", tolower(d)),
    studlab     = make_studlab(author_core, year)
  ) |>
  group_by(outcome_group, studlab) |>
  summarise(
    n_ctrl = first(n_ctrl),
    n_int  = first(n_int),
    mean_ctrl_harm = mean(mean_ctrl_harm, na.rm = TRUE),
    mean_int_harm  = mean(mean_int_harm,  na.rm = TRUE),
    sd_ctrl = mean(sd_ctrl, na.rm = TRUE),
    sd_int  = mean(sd_int,  na.rm = TRUE),
    .groups = "drop"
  )


#### 6. Meta-analysis per outcome group ---------------------------------

run_meta_group <- function(data, group, sm = "SMD") {
  
  dat <- data |>
    filter(outcome_group == group) |>
    filter(
      !is.na(mean_int_harm), !is.na(mean_ctrl_harm),
      !is.na(sd_int), !is.na(sd_ctrl),
      n_int  > 0, n_ctrl > 0
    )
  
  if (nrow(dat) == 0L) {
    stop(paste0("No rows for outcome_group = '", group,
                "' in post measurements with control group."))
  }
  
  m <- metacont(
    n.e    = dat$n_ctrl,
    mean.e = dat$mean_ctrl_harm,
    sd.e   = dat$sd_ctrl,
    n.c    = dat$n_int,
    mean.c = dat$mean_int_harm,
    sd.c   = dat$sd_int,
    studlab = dat$studlab,
    sm      = sm,
    common  = FALSE,
    random  = TRUE,
    method.tau       = "REML",
    method.random.ci = "classic",
    prediction       = FALSE
  )
  
  return(m)
}


#### 7. List of meta objects by outcome_group ---------------------------

outcome_groups <- sort(unique(dat_for_meta$outcome_group))

meta_list <- lapply(outcome_groups, function(g) {
  run_meta_group(dat_for_meta, group = g)
})
names(meta_list) <- outcome_groups


#### 8. Nice English labels for outcomes -------------------------------

label_group_en <- function(group) {
  dplyr::case_when(
    group == "caregiver_burden_distress" ~
      "Caregiver burden / distress",
    group == "caregiver_emotional_symptoms_mental_health" ~
      "Caregiver emotional symptoms / mental health",
    group == "caregiver_quality_of_life" ~
      "Caregiver quality of life (including sleep)",
    group == "patient_neuropsychiatric_symptoms_behavior" ~
      "Patient neuropsychiatric symptoms / behaviour",
    group == "caregiver_stress" ~
      "Caregiver stress",
    group == "caregiver_physical_general_health" ~
      "Caregiver physical / general health",
    group == "caregiver_self_efficacy_competence" ~
      "Caregiver self-efficacy / competence",
    group == "caregiver_coping_resources" ~
      "Caregiver coping and resources",
    group == "dementia_knowledge_attitudes" ~
      "Dementia knowledge and attitudes",
    group == "patient_daily_functioning" ~
      "Patient daily functioning (ADL / IADL)",
    group == "patient_cognitive_function" ~
      "Patient cognitive function",
    TRUE ~ group
  )
}

sanitize_name <- function(x) {
  x %>%
    tolower() %>%
    gsub("[^a-z0-9]+", "_", .) %>%
    gsub("_+", "_", .) %>%
    gsub("^_|_$", "", .)
}


#### 9. Forest + funnel -------------------------------------------------

save_plots_for_group <- function(m, group, out_dir) {
  
  group_clean <- sanitize_name(group)
  group_label <- label_group_en(group)   # plot title
  
  k <- length(m$TE)
  # algo más alto por estudio para que haya aire
  height_in <- max(5.5, 3.0 + 0.30 * k)
  
  ## FOREST -------------------------------------------------------------
  
  forest_file <- file.path(out_dir, paste0("forest_", group_clean, ".tiff"))
  
  message("Saving forest for ", group_label,
          " (", group, "), k = ", k,
          " → ", forest_file)
  
  tiff(
    filename    = forest_file,
    width       = 8.5,
    height      = height_in,
    units       = "in",
    res         = 600,
    compression = "lzw"
  )
  
  par(mar = c(5.5, 4.5, 4, 2))
  
  forest(
    m,
    main        = group_label,
    xlab        = "Standardized mean difference (Hedges g)",
    smlab       = "",
    prediction  = FALSE,  # sin intervalo extra: solo diamante (IC95%)
    common      = FALSE,
    random      = TRUE,
    leftcols    = c("studlab"),
    leftlabs    = c("Study"),
    
    ## --- SOLO I2, tau2 y p de heterogeneidad ---
    print.I2       = TRUE,
    print.I2.ci    = FALSE,
    print.tau2     = TRUE,
    print.tau2.ci  = FALSE,
    print.Q        = FALSE,
    print.pval.Q   = TRUE,
    
    ## tamaño general
    fontsize       = 8,
    ## fuente de la línea de heterogeneidad
    fs.hetstat     = 4,
    
    ## más espacio entre columnas y el bosque
    colgap.forest.left  = "0.10cm",
    colgap.forest.right = "0.60cm"
    # sin "favours" text
  )
  
  dev.off()
  
  ## FUNNEL -------------------------------------------------------------
  
  funnel_file <- file.path(out_dir, paste0("funnel_", group_clean, ".tiff"))
  
  message("Saving funnel for ", group_label,
          " (", group, "), k = ", k,
          " → ", funnel_file)
  
  tiff(
    filename    = funnel_file,
    width       = 7,
    height      = 7,
    units       = "in",
    res         = 600,
    compression = "lzw"
  )
  
  par(mar = c(5, 4.5, 4, 2))
  
  funnel(
    m,
    main    = group_label,
    xlab    = "Standardized mean difference (Hedges g)",
    ylab    = "Standard error",
    studlab = FALSE
  )
  
  dev.off()
}


#### 10. Loop: save all plots -------------------------------------------

for (g in outcome_groups) {
  m <- meta_list[[g]]
  save_plots_for_group(m, g, out_dir)
}

#### 11. Global forest: all outcomes, by study and instrument ----------

# Column that contains the instrument / scale name
instrument_col <- "outcome_scale"

# Safety check
if (!instrument_col %in% names(dat_control_post)) {
  stop("The column specified in instrument_col ('", instrument_col,
       "') is not in dat_control_post. Run names(dat_control_post) to see available columns.")
}

# 1) Build a long data frame: one row = study × instrument × outcome
dat_all_effects <- dat_control_post |>
  filter(
    !is.na(mean_int_harm), !is.na(mean_ctrl_harm),
    !is.na(sd_int), !is.na(sd_ctrl),
    n_int  > 0, n_ctrl > 0
  ) |>
  mutate(
    author_core      = gsub("[^a-z]+.*$", "", tolower(d)),
    author_label     = str_to_title(author_core),
    instrument_label = .data[[instrument_col]],
    studlab_instr    = paste0(author_label, " et al., ", year,
                              " (", instrument_label, ")")
  )

# 2) Meta-analysis using ALL study × instrument rows
m_all <- metacont(
  n.e    = dat_all_effects$n_ctrl,
  mean.e = dat_all_effects$mean_ctrl_harm,
  sd.e   = dat_all_effects$sd_ctrl,
  n.c    = dat_all_effects$n_int,
  mean.c = dat_all_effects$mean_int_harm,
  sd.c   = dat_all_effects$sd_int,
  studlab = dat_all_effects$studlab_instr,
  sm      = "SMD",
  common  = FALSE,
  random  = TRUE,
  method.tau       = "REML",
  method.random.ci = "classic",
  prediction       = FALSE
)

# 3) Plot height based on k
k_all      <- length(m_all$TE)
height_all <- max(6, 3.0 + 0.30 * k_all)

forest_all_file <- file.path(out_dir, "forest_all_outcomes_instruments.tiff")

message("Saving global forest with all outcomes → ", forest_all_file)

tiff(
  filename    = forest_all_file,
  width       = 8.5,
  height      = height_all,
  units       = "in",
  res         = 600,
  compression = "lzw"
)

par(mar = c(5.5, 4.5, 4, 2))

forest(
  m_all,
  sortvar     = -m_all$TE,   # largest positive effects at the top
  main        = "All outcomes by study and instrument",
  xlab        = "Standardized mean difference (Hedges g)",
  smlab       = "",
  leftcols    = c("studlab"),
  leftlabs    = c("Study (instrument)"),
  overall     = FALSE,       # sin diamante global
  prediction  = FALSE,
  common      = FALSE,
  random      = TRUE,
  
  print.I2       = TRUE,
  print.I2.ci    = FALSE,
  print.tau2     = TRUE,
  print.tau2.ci  = FALSE,
  print.Q        = FALSE,
  print.pval.Q   = TRUE,
  
  fontsize       = 8,
  fs.hetstat     = 4,
  colgap.forest.left  = "0.10cm",
  colgap.forest.right = "0.60cm"
)

dev.off()

