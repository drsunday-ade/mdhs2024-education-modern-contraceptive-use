#!/usr/bin/env Rscript
# ======================================================================
# Malawi DHS 2024 (Standard DHS) — H3 Manuscript Pack (FINAL FIXED)
# Study: Maternal education gradients in modern contraceptive use among women in union
# Target journals: Lancet Global Health / BMJ Sexual & Reproductive Health
#
# INPUT (per your README; folders beside this script):
#   MWIR81DT/MWIR81FL.dta   (Individual Recode, IR)  <-- primary analytic file
#
# OUTPUTS (OUT_ROOT):
#   figures_png/ (650 dpi, auto-upsize to >=1.5MB)
#     Figure1_Flow.png
#     Figure2_UnadjustedPrev_byEducation.png
#     Figure3_AdjustedPrev_byEducation.png
#     Figure4_aPR_Forest_Education.png
#     Figure5_AdjustedPrev_Education_byResidence.png
#   tables_csv/
#     Table1_SampleCharacteristics_byEducation.csv
#     Table2_ModernPrev_byEducation_andWealth.csv
#     Table3_Regression_aPR_Models.csv
#     Table4_AdjustedPrev_Standardized.csv
#     Table5_EffectModification_Residence.csv
#   tables_html/ (same tables as HTML)
#   logs/runlog.txt + logs/sessionInfo.txt
#
# Key DHS definitions (IR):
# - women in union: v502 in {1 married, 2 living with partner}
# - modern method use: v313==3 (modern method type)
# - survey design: weight v005/1e6, PSU v021, strata v022 (fallback v023)
#
# CORE FIXES:
# - Do NOT use survey::update() (not exported in some installations).
# - Add variables to a design by editing design$variables (safe across versions).
# - Ensure all CI arithmetic is numeric BEFORE percent formatting.
# ======================================================================

# ---------------------------
# 0) stability + args
# ---------------------------
Sys.setenv(OMP_NUM_THREADS="1", OPENBLAS_NUM_THREADS="1", MKL_NUM_THREADS="1", VECLIB_MAXIMUM_THREADS="1")
options(stringsAsFactors = FALSE, scipen = 999, width = 140, warn = 1)
set.seed(20260301)

args <- commandArgs(trailingOnly = TRUE)
get_arg <- function(flag, default=NULL) {
  hit <- which(args == flag)
  if (!length(hit) || hit == length(args)) return(default)
  args[[hit + 1]]
}
get_script_dir <- function() {
  ca <- commandArgs(trailingOnly = FALSE)
  file_arg <- ca[grepl("^--file=", ca)]
  if (length(file_arg) == 1) {
    p <- sub("^--file=", "", file_arg)
    return(normalizePath(dirname(p), winslash="/", mustWork=FALSE))
  }
  normalizePath(getwd(), winslash="/", mustWork=FALSE)
}

PROJECT_ROOT <- normalizePath(get_arg("--project_root", get_script_dir()), winslash="/", mustWork=FALSE)
OUT_ROOT <- normalizePath(get_arg("--out_root", file.path(PROJECT_ROOT, "outputs_H3_manuscript_pack_FINAL")),
                          winslash="/", mustWork=FALSE)
AUTO_INSTALL <- tolower(get_arg("--auto_install","true")) %in% c("true","1","yes","y")

setwd(PROJECT_ROOT)

# ---------------------------
# 1) packages
# ---------------------------
install_if_missing <- function(pkgs) {
  miss <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)]
  if (length(miss) && AUTO_INSTALL) {
    message("Installing: ", paste(miss, collapse=", "))
    install.packages(miss, repos="https://cloud.r-project.org")
  }
}

install_if_missing(c(
  "fs","readr","stringr","purrr","dplyr","tidyr","tibble",
  "haven","janitor","labelled",
  "survey","broom",
  "ggplot2","scales","ggrepel",
  "splines",
  "gt"
))

suppressPackageStartupMessages({
  library(fs); library(readr); library(stringr); library(purrr); library(dplyr); library(tidyr); library(tibble)
  library(haven); library(janitor); library(labelled)
  library(survey); library(broom)
  library(ggplot2); library(scales); library(ggrepel)
  library(splines)
  library(gt)
})

HAS_RAGG <- requireNamespace("ragg", quietly = TRUE)
if (!HAS_RAGG) install_if_missing("ragg")
HAS_RAGG <- requireNamespace("ragg", quietly = TRUE)

# ---------------------------
# 2) outputs + logging
# ---------------------------
dir_create(OUT_ROOT)
dir_create(path(OUT_ROOT,"logs"))
dir_create(path(OUT_ROOT,"figures_png"))
dir_create(path(OUT_ROOT,"tables_csv"))
dir_create(path(OUT_ROOT,"tables_html"))

LOG <- path(OUT_ROOT,"logs","runlog.txt")
log_line <- function(...) {
  msg <- paste0(..., collapse="")
  cat(format(Sys.time(), "%Y-%m-%d %H:%M:%S"), " | ", msg, "\n", file=LOG, append=TRUE)
  invisible(msg)
}

with_step <- function(step_name, expr, fatal=TRUE) {
  log_line("STEP START: ", step_name)
  exprq <- substitute(expr)
  out <- tryCatch(
    eval(exprq, envir=parent.frame()),
    error=function(e){
      log_line("STEP ERROR: ", step_name, " | ", conditionMessage(e))
      writeLines(c(
        paste0("Step: ", step_name),
        paste0("Message: ", conditionMessage(e)),
        "",
        paste(capture.output(traceback()), collapse="\n")
      ), path(OUT_ROOT,"logs", paste0("error_", gsub("[^A-Za-z0-9]+","_",step_name), ".txt")))
      if (fatal) stop(e)
      structure(list(ok=FALSE, error=conditionMessage(e)), class="step_failed")
    }
  )
  log_line("STEP END: ", step_name)
  out
}

# ---------------------------
# 3) helpers (robust)
# ---------------------------
req_vars <- function(df, vars, where="dataset") {
  miss <- vars[!(vars %in% names(df))]
  if (length(miss)) stop(where, " missing required variables: ", paste(miss, collapse=", "))
  TRUE
}
to_num <- function(x) suppressWarnings(as.numeric(as.character(x)))

# Robustly add a variable to a survey design WITHOUT update()/survey::update()
design_add_var <- function(des, name, vec) {
  if (is.null(des$variables) || !is.data.frame(des$variables)) {
    stop("Design object does not contain a data.frame in $variables.")
  }
  if (length(vec) != nrow(des$variables)) {
    stop("Length mismatch adding '", name, "': vec length=", length(vec),
         " but design n=", nrow(des$variables))
  }
  des2 <- des
  des2$variables[[name]] <- vec
  des2
}

# DHS survey design for IR
make_design_ir <- function(df) {
  strata_var <- if ("v022" %in% names(df)) "v022" else if ("v023" %in% names(df)) "v023" else NA_character_
  if (is.na(strata_var)) stop("Cannot find strata variable (need v022 or v023).")
  
  df <- df %>%
    mutate(
      wt = .data[["v005"]] / 1e6,
      psu = .data[["v021"]],
      strata = .data[[strata_var]]
    ) %>%
    filter(is.finite(wt), !is.na(wt), wt > 0, !is.na(psu), !is.na(strata))
  
  options(survey.lonely.psu = "adjust")
  survey::svydesign(ids=~psu, strata=~strata, weights=~wt, nest=TRUE, data=df)
}

# Survey-weighted PR preferred; OR fallback
safe_svyglm <- function(des, fml) {
  warn <- character(0)
  handler <- function(w) { warn <<- c(warn, conditionMessage(w)); invokeRestart("muffleWarning") }
  
  m1 <- withCallingHandlers(
    tryCatch(svyglm(fml, design=des, family=quasipoisson(link="log"), control=glm.control(maxit=100)),
             error=function(e) NULL),
    warning=handler
  )
  if (!is.null(m1)) return(list(model=m1, scale="PR", link="log", warnings=unique(warn)))
  
  warn <- character(0)
  m2 <- withCallingHandlers(
    tryCatch(svyglm(fml, design=des, family=quasibinomial(link="logit"), control=glm.control(maxit=100)),
             error=function(e) NULL),
    warning=handler
  )
  if (!is.null(m2)) return(list(model=m2, scale="OR", link="logit", warnings=unique(warn)))
  
  NULL
}

tidy_exp <- function(fit) {
  if (is.null(fit) || is.null(fit$model)) return(tibble())
  broom::tidy(fit$model, conf.int=TRUE) %>%
    mutate(effect = exp(estimate), lcl = exp(conf.low), ucl = exp(conf.high))
}

wald_p <- function(fit, term) {
  if (is.null(fit) || is.null(fit$model)) return(NA_real_)
  tryCatch({
    rt <- survey::regTermTest(fit$model, as.formula(paste0("~", term)))
    as.numeric(rt$p)
  }, error=function(e) NA_real_)
}

# Marginal standardization (numeric CI first; format later)
std_prev_levels <- function(des, fit, varname, levels_vec) {
  stopifnot(!is.null(fit$model))
  base_dat <- des$variables
  out <- vector("list", length(levels_vec))
  
  for (i in seq_along(levels_vec)) {
    lv <- levels_vec[[i]]
    nd <- base_dat
    nd[[varname]] <- factor(lv, levels=levels_vec)
    
    pred <- as.numeric(predict(fit$model, newdata=nd, type="response"))
    pred <- pmin(1, pmax(0, pred))  # keep in [0,1] for plotting/CI
    
    des2 <- design_add_var(des, "pred_tmp", pred)
    m <- tryCatch(svymean(~pred_tmp, des2, na.rm=TRUE), error=function(e) NULL)
    
    if (is.null(m)) {
      out[[i]] <- tibble(level=lv, est=NA_real_, se=NA_real_, lcl=NA_real_, ucl=NA_real_)
    } else {
      est <- as.numeric(m)[1]
      se  <- sqrt(as.numeric(vcov(m))[1])
      z <- qnorm(0.975)
      out[[i]] <- tibble(level=lv, est=est, se=se, lcl=pmax(0, est - z*se), ucl=pmin(1, est + z*se))
    }
  }
  bind_rows(out)
}

# High-DPI PNG save w/ minimum size assurance
save_png650 <- function(filename, plot, width_in=8.2, height_in=5.6, dpi=650, min_bytes=1.5e6) {
  fp <- path(OUT_ROOT, "figures_png", filename)
  
  save_once <- function(w, h) {
    if (HAS_RAGG) {
      ragg::agg_png(filename=fp, width=w, height=h, units="in", res=dpi, background="white")
      print(plot)
      dev.off()
    } else {
      ggsave(fp, plot=plot, width=w, height=h, dpi=dpi, units="in",
             bg="white", limitsize=FALSE, device="png", type="cairo-png")
    }
    file.size(fp)
  }
  
  sz <- save_once(width_in, height_in)
  if (!is.na(sz) && sz < min_bytes) sz <- save_once(width_in + 1.4, height_in + 1.0)
  if (!is.na(sz) && sz < min_bytes) sz <- save_once(width_in + 2.8, height_in + 2.0)
  
  log_line("Saved figure: ", fp, " | bytes=", ifelse(is.na(sz), "NA", format(sz, scientific=FALSE)))
  invisible(fp)
}

# Journal theme (no clipping/overlap)
theme_journal <- function() {
  theme_minimal(base_size = 12, base_family = "sans") +
    theme(
      plot.title.position = "plot",
      plot.caption.position = "plot",
      plot.title = element_text(face="bold", size=14, margin=margin(b=6)),
      plot.subtitle = element_text(size=11, margin=margin(b=8)),
      plot.caption = element_text(size=9, color="grey30", margin=margin(t=10)),
      axis.title = element_text(size=11),
      axis.text = element_text(size=10),
      legend.title = element_text(size=10, face="bold"),
      legend.text = element_text(size=9),
      panel.grid.minor = element_blank(),
      panel.grid.major.x = element_blank(),
      plot.margin = margin(16, 20, 16, 20)
    )
}

# Palettes (print-safe, journal-like)
pal_edu <- c("No education"="#243B53","Primary"="#3C6E71","Secondary"="#C8553D","Higher"="#8E44AD")
pal_res <- c("Urban"="#1F77B4","Rural"="#2CA02C")

# ---------------------------
# 4) Read IR file
# ---------------------------
IR_PATH <- file.path(PROJECT_ROOT, "MWIR81DT", "MWIR81FL.dta")
if (!file.exists(IR_PATH)) {
  hits <- tryCatch(fs::dir_ls(file.path(PROJECT_ROOT, "MWIR81DT"), recurse=TRUE, type="file"),
                   error=function(e) character(0))
  hit <- hits[tolower(basename(hits)) == "mwir81fl.dta"]
  if (length(hit)) IR_PATH <- hit[[1]]
}
if (!file.exists(IR_PATH)) stop("Cannot find MWIR81FL.dta under MWIR81DT/")

IR <- with_step("Read IR", janitor::clean_names(haven::read_dta(IR_PATH)), fatal=TRUE)
req_vars(IR, c("v005","v021","v012","v024","v025","v106","v190","v313","v502"), where="IR")

# ---------------------------
# 5) Build analytic cohort (women in union) + recodes
# ---------------------------
dat0 <- with_step("Build analytic dataset", {
  d <- IR
  
  # Eligibility: women in union (v502 1=married, 2=living with partner)
  v502 <- to_num(d$v502)
  d$in_union <- ifelse(is.na(v502), NA_integer_, ifelse(v502 %in% c(1,2), 1L, 0L))
  
  # Outcome: modern contraception (v313==3)
  v313 <- to_num(d$v313)
  d$modern_use <- ifelse(is.na(v313), NA_integer_,
                         ifelse(v313 == 3, 1L,
                                ifelse(v313 %in% c(0,1,2), 0L, NA_integer_)))
  
  # Education (v106: 0 none, 1 primary, 2 secondary, 3 higher)
  v106 <- to_num(d$v106)
  d$educ <- dplyr::case_when(
    v106 == 0 ~ "No education",
    v106 == 1 ~ "Primary",
    v106 == 2 ~ "Secondary",
    v106 == 3 ~ "Higher",
    TRUE ~ NA_character_
  )
  d$educ <- factor(d$educ, levels=c("No education","Primary","Secondary","Higher"))
  
  # Wealth quintile (v190: 1 poorest ... 5 richest)
  v190 <- to_num(d$v190)
  d$wealth_q <- dplyr::case_when(
    v190 == 1 ~ "Q1 Poorest",
    v190 == 2 ~ "Q2 Poorer",
    v190 == 3 ~ "Q3 Middle",
    v190 == 4 ~ "Q4 Richer",
    v190 == 5 ~ "Q5 Richest",
    TRUE ~ NA_character_
  )
  d$wealth_q <- factor(d$wealth_q, levels=c("Q1 Poorest","Q2 Poorer","Q3 Middle","Q4 Richer","Q5 Richest"))
  
  # Residence (v025: 1 urban, 2 rural)
  v025 <- to_num(d$v025)
  d$residence <- dplyr::case_when(
    v025 == 1 ~ "Urban",
    v025 == 2 ~ "Rural",
    TRUE ~ NA_character_
  )
  d$residence <- factor(d$residence, levels=c("Urban","Rural"))
  
  # Region label (v024), keep DHS labels if present
  d$region <- labelled::to_factor(d$v024, levels="labels", sort_levels="none")
  d$region <- as.factor(d$region)
  
  # Age
  d$age <- to_num(d$v012)
  
  # Keep target group only
  d <- d %>% filter(in_union == 1)
  
  d
}, fatal=TRUE)

# Flow counts (unweighted)
flow_counts <- with_step("Compute flow counts", {
  n_total <- nrow(IR)
  n_union <- nrow(IR %>%
                    mutate(in_union = ifelse(to_num(v502) %in% c(1,2), 1L, 0L)) %>%
                    filter(in_union==1))
  n_union_nonmiss_outcome <- nrow(dat0 %>% filter(!is.na(modern_use)))
  n_cc <- nrow(dat0 %>% filter(!is.na(modern_use), !is.na(educ), !is.na(wealth_q),
                               !is.na(residence), !is.na(region), !is.na(age)))
  tibble(
    step=c("All women (IR, 15–49)",
           "Women in union (v502∈{1,2})",
           "Non-missing outcome (modern method status from v313)",
           "Complete-case analytic sample"),
    n=c(n_total, n_union, n_union_nonmiss_outcome, n_cc)
  )
}, fatal=TRUE)

write_csv(flow_counts, path(OUT_ROOT,"tables_csv","FlowCounts_unweighted.csv"))

# Complete-case analytic dataset
dat_cc <- dat0 %>%
  filter(!is.na(modern_use), !is.na(educ), !is.na(wealth_q),
         !is.na(residence), !is.na(region), !is.na(age))

# Survey design
des_cc <- with_step("Create survey design (complete-case)", make_design_ir(dat_cc), fatal=TRUE)

# ---------------------------
# 6) Models (match screening + publishable sensitivity)
# ---------------------------
fits <- with_step("Fit models", {
  # Model S0: EXACT screening-like spec (age linear) -> should match prior significance
  fml_S0 <- modern_use ~ educ + wealth_q + residence + region + age
  fit_S0 <- safe_svyglm(des_cc, fml_S0); if (is.null(fit_S0)) stop("Model S0 failed.")
  
  # Model S1: flexible age (spline), publishable robustness
  fml_S1 <- modern_use ~ educ + wealth_q + residence + region + splines::ns(age, df=4)
  fit_S1 <- safe_svyglm(des_cc, fml_S1); if (is.null(fit_S1)) stop("Model S1 failed.")
  
  # Interaction: education × residence (age linear for stability)
  fml_INT <- modern_use ~ educ * residence + wealth_q + region + age
  fit_INT <- safe_svyglm(des_cc, fml_INT); if (is.null(fit_INT)) stop("Model INT failed.")
  
  list(fit_S0=fit_S0, fit_S1=fit_S1, fit_INT=fit_INT)
}, fatal=TRUE)

fit_S0 <- fits$fit_S0
fit_S1 <- fits$fit_S1
fit_INT <- fits$fit_INT

p_educ_S0 <- wald_p(fit_S0, "educ")
p_wealth_S0 <- wald_p(fit_S0, "wealth_q")

p_educ_S1 <- wald_p(fit_S1, "educ")
p_wealth_S1 <- wald_p(fit_S1, "wealth_q")

p_int <- wald_p(fit_INT, "educ:residence")

log_line("Global Wald p-values | S0 educ=", signif(p_educ_S0,3), " wealth=", signif(p_wealth_S0,3),
         " | S1 educ=", signif(p_educ_S1,3), " wealth=", signif(p_wealth_S1,3),
         " | INT educ×residence=", signif(p_int,3))

# ---------------------------
# 7) TABLES (5) — CSV + HTML
# ---------------------------

# ---- Table 1: Sample characteristics by education (weighted)
Table1 <- with_step("Table 1", {
  tab_educ <- svytable(~educ, des_cc)
  share <- prop.table(tab_educ)
  
  # Weighted prevalence of modern_use by education
  prev <- svyby(~modern_use, ~educ, des_cc, svymean, na.rm=TRUE, vartype=c("se","ci"))
  prev_df <- as_tibble(prev) %>%
    transmute(educ=as.character(educ),
              modern_prev=modern_use, modern_lcl=ci_l, modern_ucl=ci_u)
  
  # Weighted mean age by education
  age_m <- svyby(~age, ~educ, des_cc, svymean, na.rm=TRUE, vartype=c("se","ci"))
  age_df <- as_tibble(age_m) %>%
    transmute(educ=as.character(educ), age_mean=age, age_lcl=ci_l, age_ucl=ci_u)
  
  # Urban share by education (residence==Urban)
  urban_vec <- ifelse(as.character(des_cc$variables$residence) == "Urban", 1, 0)
  des_urb <- design_add_var(des_cc, "urban", urban_vec)
  urb <- svyby(~urban, ~educ, des_urb, svymean, na.rm=TRUE, vartype=c("se","ci"))
  urb_df <- as_tibble(urb) %>%
    transmute(educ=as.character(educ), urban_share=urban, urban_lcl=ci_l, urban_ucl=ci_u)
  
  share_df <- tibble(
    educ = names(tab_educ),
    n_unw = as.integer(table(dat_cc$educ)[names(tab_educ)]),
    pop_share_wt = as.numeric(share)
  )
  
  out <- share_df %>%
    left_join(age_df, by="educ") %>%
    left_join(urb_df, by="educ") %>%
    left_join(prev_df, by="educ") %>%
    mutate(
      pop_share_wt = percent(pop_share_wt, accuracy=0.1),
      age_mean = round(age_mean, 1),
      age_lcl = round(age_lcl, 1),
      age_ucl = round(age_ucl, 1),
      urban_share = percent(urban_share, accuracy=0.1),
      urban_lcl = percent(urban_lcl, accuracy=0.1),
      urban_ucl = percent(urban_ucl, accuracy=0.1),
      modern_prev = percent(modern_prev, accuracy=0.1),
      modern_lcl = percent(modern_lcl, accuracy=0.1),
      modern_ucl = percent(modern_ucl, accuracy=0.1)
    )
  
  out
}, fatal=TRUE)

write_csv(Table1, path(OUT_ROOT,"tables_csv","Table1_SampleCharacteristics_byEducation.csv"))
gtsave(
  Table1 %>%
    gt() %>%
    tab_header(
      title="Table 1. Sample characteristics by maternal education among women in union, Malawi DHS 2024",
      subtitle="Survey-weighted population shares and indicators; 95% CIs shown for age, urban share, and modern contraception prevalence."
    ) %>%
    cols_label(
      educ="Education",
      n_unw="Unweighted N",
      pop_share_wt="Population share (weighted)",
      age_mean="Mean age (years)",
      age_lcl="Age 95% L",
      age_ucl="Age 95% U",
      urban_share="Urban (%)",
      urban_lcl="Urban 95% L",
      urban_ucl="Urban 95% U",
      modern_prev="Modern contraception (%)",
      modern_lcl="Modern 95% L",
      modern_ucl="Modern 95% U"
    ) %>%
    opt_table_lines("none") %>%
    tab_options(table.font.size=12),
  path(OUT_ROOT,"tables_html","Table1_SampleCharacteristics_byEducation.html")
)

# ---- Table 2: Modern contraception prevalence by education × wealth (weighted)
Table2 <- with_step("Table 2", {
  prev2 <- svyby(~modern_use, ~educ + wealth_q, des_cc, svymean, na.rm=TRUE, vartype=c("se","ci"))
  as_tibble(prev2) %>%
    transmute(
      educ=as.character(educ),
      wealth_q=as.character(wealth_q),
      modern_prev=percent(modern_use, accuracy=0.1),
      lcl=percent(ci_l, accuracy=0.1),
      ucl=percent(ci_u, accuracy=0.1)
    )
}, fatal=TRUE)

write_csv(Table2, path(OUT_ROOT,"tables_csv","Table2_ModernPrev_byEducation_andWealth.csv"))
gtsave(
  Table2 %>%
    gt(groupname_col="educ") %>%
    tab_header(
      title="Table 2. Modern contraception prevalence by education and wealth quintile (women in union), Malawi DHS 2024",
      subtitle="Survey-weighted prevalence with 95% confidence intervals."
    ) %>%
    cols_label(
      wealth_q="Wealth quintile",
      modern_prev="Modern contraception (%)",
      lcl="95% L",
      ucl="95% U"
    ) %>%
    opt_table_lines("none") %>%
    tab_options(table.font.size=12),
  path(OUT_ROOT,"tables_html","Table2_ModernPrev_byEducation_andWealth.html")
)

# ---- Table 3: Regression results (aPR) — screening model + spline sensitivity
Table3 <- with_step("Table 3", {
  td0 <- tidy_exp(fit_S0) %>% mutate(model=paste0("Model S0 (screening spec; ", fit_S0$scale, ")"))
  td1 <- tidy_exp(fit_S1) %>% mutate(model=paste0("Model S1 (age spline; ", fit_S1$scale, ")"))
  
  td <- bind_rows(td0, td1) %>%
    filter(term != "(Intercept)") %>%
    mutate(
      effect=round(effect, 3),
      lcl=round(lcl, 3),
      ucl=round(ucl, 3),
      p_value=signif(p.value, 3),
      term_clean = term
    ) %>%
    mutate(
      term_clean = str_replace(term_clean, "^educ", "Education: "),
      term_clean = str_replace(term_clean, "^wealth_q", "Wealth: "),
      term_clean = str_replace(term_clean, "^residence", "Residence: "),
      term_clean = str_replace(term_clean, "^region", "Region: "),
      term_clean = str_replace(term_clean, "splines::ns\\(age, df = 4\\)", "Age (spline basis)"),
      term_clean = str_replace(term_clean, "^age$", "Age (years)")
    ) %>%
    select(model, term=term_clean, effect, lcl, ucl, p_value)
  
  td
}, fatal=TRUE)

write_csv(Table3, path(OUT_ROOT,"tables_csv","Table3_Regression_aPR_Models.csv"))
gtsave(
  Table3 %>%
    gt(groupname_col="model") %>%
    tab_header(
      title="Table 3. Association of maternal education and wealth with modern contraceptive use (women in union), Malawi DHS 2024",
      subtitle=paste0(
        "Survey-weighted models (PR preferred; OR fallback). Global Wald p (S0): educ=",
        signif(p_educ_S0,3), ", wealth=", signif(p_wealth_S0,3),
        ". Global Wald p (S1): educ=", signif(p_educ_S1,3), ", wealth=", signif(p_wealth_S1,3), "."
      )
    ) %>%
    cols_label(term="Predictor level", effect="Effect", lcl="95% L", ucl="95% U", p_value="p-value") %>%
    opt_table_lines("none") %>%
    tab_options(table.font.size=12),
  path(OUT_ROOT,"tables_html","Table3_Regression_aPR_Models.html")
)

# ---- Table 4: Standardized adjusted prevalence by education (S0 vs S1)
Table4_num <- with_step("Table 4 (numeric)", {
  edu_levels <- levels(dat_cc$educ)
  std0 <- std_prev_levels(des_cc, fit_S0, "educ", edu_levels) %>% mutate(model="Model S0 (screening spec)")
  std1 <- std_prev_levels(des_cc, fit_S1, "educ", edu_levels) %>% mutate(model="Model S1 (age spline)")
  bind_rows(std0, std1) %>% rename(educ=level)
}, fatal=TRUE)

Table4 <- Table4_num %>%
  mutate(
    est = percent(est, accuracy=0.1),
    lcl = percent(lcl, accuracy=0.1),
    ucl = percent(ucl, accuracy=0.1)
  ) %>%
  select(model, educ, est, lcl, ucl)

write_csv(Table4, path(OUT_ROOT,"tables_csv","Table4_AdjustedPrev_Standardized.csv"))
gtsave(
  Table4 %>%
    gt(groupname_col="model") %>%
    tab_header(
      title="Table 4. Standardized adjusted prevalence of modern contraception by education (women in union), Malawi DHS 2024",
      subtitle="Predicted prevalence standardized to the analytic distribution; 95% CIs via survey linearization."
    ) %>%
    cols_label(educ="Education", est="Adjusted prevalence (%)", lcl="95% L", ucl="95% U") %>%
    opt_table_lines("none") %>%
    tab_options(table.font.size=12),
  path(OUT_ROOT,"tables_html","Table4_AdjustedPrev_Standardized.html")
)

# ---- Table 5: Effect modification by residence (interaction model; standardized)
Table5_num <- with_step("Table 5 (numeric)", {
  edu_levels <- levels(dat_cc$educ)
  res_levels <- levels(dat_cc$residence)
  
  base_dat <- des_cc$variables
  out <- list()
  k <- 1
  
  for (r0 in res_levels) {
    nd_base <- base_dat
    nd_base$residence <- factor(r0, levels=res_levels)
    
    for (e0 in edu_levels) {
      nd <- nd_base
      nd$educ <- factor(e0, levels=edu_levels)
      
      pred <- as.numeric(predict(fit_INT$model, newdata=nd, type="response"))
      pred <- pmin(1, pmax(0, pred))
      
      des2 <- design_add_var(des_cc, "pred_tmp", pred)
      m <- tryCatch(svymean(~pred_tmp, des2, na.rm=TRUE), error=function(e) NULL)
      
      if (is.null(m)) {
        out[[k]] <- tibble(residence=r0, educ=e0, est=NA_real_, lcl=NA_real_, ucl=NA_real_)
      } else {
        est <- as.numeric(m)[1]
        se  <- sqrt(as.numeric(vcov(m))[1])
        z <- qnorm(0.975)
        out[[k]] <- tibble(
          residence=r0, educ=e0,
          est=est,
          lcl=pmax(0, est - z*se),
          ucl=pmin(1, est + z*se)
        )
      }
      k <- k + 1
    }
  }
  
  bind_rows(out) %>% mutate(interaction_p = p_int)
}, fatal=TRUE)

Table5 <- Table5_num %>%
  mutate(
    est = percent(est, accuracy=0.1),
    lcl = percent(lcl, accuracy=0.1),
    ucl = percent(ucl, accuracy=0.1),
    interaction_p = signif(interaction_p, 3)
  ) %>%
  select(residence, educ, est, lcl, ucl, interaction_p)

write_csv(Table5, path(OUT_ROOT,"tables_csv","Table5_EffectModification_Residence.csv"))
gtsave(
  Table5 %>%
    gt(groupname_col="residence") %>%
    tab_header(
      title="Table 5. Education gradients in modern contraceptive use by residence (interaction model), Malawi DHS 2024",
      subtitle=paste0("Standardized adjusted prevalence within residence strata. Global Wald p for educ×residence interaction: ",
                      signif(p_int,3), ".")
    ) %>%
    cols_label(
      educ="Education",
      est="Adjusted prevalence (%)",
      lcl="95% L",
      ucl="95% U",
      interaction_p="Interaction p"
    ) %>%
    opt_table_lines("none") %>%
    tab_options(table.font.size=12),
  path(OUT_ROOT,"tables_html","Table5_EffectModification_Residence.html")
)

# ---------------------------
# 8) FIGURES (5) — 650 dpi, >=1.5MB, no clipping
# ---------------------------

# Figure 1 — flow diagram
Fig1 <- with_step("Figure 1", {
  fc <- flow_counts
  box <- tibble(
    x = 0,
    y = c(4,3,2,1),
    label = paste0(fc$step, "\nN = ", format(fc$n, big.mark=","))
  )
  
  ggplot() +
    geom_label(
      data=box,
      aes(x=x, y=y, label=label),
      label.size=0.35,
      label.r=unit(0.25, "lines"),
      size=3.6,
      fill="white"
    ) +
    geom_segment(aes(x=0, xend=0, y=3.65, yend=3.35), linewidth=0.7, lineend="round") +
    geom_segment(aes(x=0, xend=0, y=2.65, yend=2.35), linewidth=0.7, lineend="round") +
    geom_segment(aes(x=0, xend=0, y=1.65, yend=1.35), linewidth=0.7, lineend="round") +
    annotate(
      "text", x=0, y=0.35,
      label=paste0(
        "Target: women in union (v502∈{1,2}). Outcome: modern method use (v313==3).\n",
        "Design: weights=v005/1e6; PSU=v021; strata=v022 (fallback v023)."
      ),
      size=3.1, color="grey20", vjust=0
    ) +
    scale_y_continuous(limits=c(0.1, 4.7), expand=expansion(mult=c(0,0))) +
    scale_x_continuous(limits=c(-1.4, 1.4), expand=expansion(mult=c(0,0))) +
    labs(
      title="Figure 1. Study population flow (women in union), Malawi DHS 2024",
      subtitle="Unweighted counts at each selection step.",
      x=NULL, y=NULL,
      caption="DHS Individual Recode (MWIR81FL)."
    ) +
    theme_void(base_family="sans") +
    theme(
      plot.title = element_text(face="bold", size=14, margin=margin(b=6)),
      plot.subtitle = element_text(size=11, margin=margin(b=10)),
      plot.caption = element_text(size=9, color="grey30", margin=margin(t=10)),
      plot.margin = margin(16, 20, 16, 20)
    )
}, fatal=TRUE)
save_png650("Figure1_Flow.png", Fig1, width_in=8.6, height_in=6.6, dpi=650)

# Figure 2 — unadjusted weighted prevalence by education
Fig2 <- with_step("Figure 2", {
  prev <- svyby(~modern_use, ~educ, des_cc, svymean, na.rm=TRUE, vartype=c("se","ci")) %>%
    as_tibble() %>%
    transmute(
      educ=factor(as.character(educ), levels=levels(dat_cc$educ)),
      est=modern_use, lcl=ci_l, ucl=ci_u
    )
  
  ggplot(prev, aes(x=educ, y=est, color=educ)) +
    geom_point(size=3.2) +
    geom_errorbar(aes(ymin=lcl, ymax=ucl), width=0.10, linewidth=0.7) +
    scale_color_manual(values=pal_edu, guide="none") +
    scale_y_continuous(labels=percent_format(accuracy=1), limits=c(0,1),
                       expand=expansion(mult=c(0.02,0.06))) +
    labs(
      title="Figure 2. Unadjusted modern contraceptive prevalence by maternal education (women in union), Malawi DHS 2024",
      subtitle="Survey-weighted prevalence with 95% CI.",
      x="Maternal education",
      y="Modern contraception prevalence",
      caption="Outcome: v313==3 (modern method type)."
    ) +
    theme_journal()
}, fatal=TRUE)
save_png650("Figure2_UnadjustedPrev_byEducation.png", Fig2, width_in=8.2, height_in=5.6, dpi=650)

# Figure 3 — standardized adjusted prevalence by education (S0 vs S1)
Fig3 <- with_step("Figure 3", {
  tbl <- Table4_num %>%
    mutate(
      model = factor(model, levels=c("Model S0 (screening spec)","Model S1 (age spline)")),
      educ  = factor(educ, levels=levels(dat_cc$educ))
    )
  
  ggplot(tbl, aes(x=educ, y=est, group=model, color=model)) +
    geom_line(linewidth=1.0) +
    geom_point(size=2.8) +
    geom_errorbar(aes(ymin=lcl, ymax=ucl), width=0.10, linewidth=0.6) +
    scale_y_continuous(labels=percent_format(accuracy=1), limits=c(0,1),
                       expand=expansion(mult=c(0.02,0.06))) +
    labs(
      title="Figure 3. Standardized adjusted prevalence of modern contraceptive use by education (women in union), Malawi DHS 2024",
      subtitle=paste0("Global Wald p (S0): educ=", signif(p_educ_S0,3), ", wealth=", signif(p_wealth_S0,3),
                      ". Global Wald p (S1): educ=", signif(p_educ_S1,3), ", wealth=", signif(p_wealth_S1,3), "."),
      x="Maternal education",
      y="Adjusted prevalence (standardized)",
      color="Model",
      caption="Standardization averages model predictions over the analytic covariate distribution."
    ) +
    theme_journal()
}, fatal=TRUE)
save_png650("Figure3_AdjustedPrev_byEducation.png", Fig3, width_in=8.6, height_in=5.8, dpi=650)

# Figure 4 — forest plot: education effects (Model S0)
Fig4 <- with_step("Figure 4", {
  td <- tidy_exp(fit_S0) %>%
    filter(str_detect(term, "^educ")) %>%
    mutate(level = str_replace(term, "^educ", "")) %>%
    mutate(level = ifelse(level=="", "No education (ref)", level))
  
  # Add explicit ref
  ref <- tibble(level="No education (ref)", effect=1, lcl=1, ucl=1, p.value=NA_real_)
  td2 <- bind_rows(ref, td %>% select(level, effect, lcl, ucl, p.value))
  
  # Order using the education factor order
  # (Model terms will be "Primary","Secondary","Higher" relative to "No education")
  order_levels <- rev(c("No education (ref)","Primary","Secondary","Higher"))
  td2$level <- factor(td2$level, levels=order_levels)
  
  ggplot(td2, aes(y=level, x=effect)) +
    geom_vline(xintercept=1, linetype="dashed", linewidth=0.6, color="grey50") +
    geom_point(size=2.9, color="#111111") +
    geom_errorbarh(aes(xmin=lcl, xmax=ucl), height=0.18, linewidth=0.7, color="#111111") +
    scale_x_continuous(trans="log10",
                       breaks=c(0.5,0.75,1,1.25,1.5,2,3),
                       labels=function(x) format(x, trim=TRUE),
                       expand=expansion(mult=c(0.02,0.10))) +
    labs(
      title="Figure 4. Adjusted prevalence ratios (aPR) for modern contraceptive use by education (women in union), Malawi DHS 2024",
      subtitle="Model S0 (screening spec): educ + wealth + residence + region + age (linear).",
      x=paste0("Adjusted ratio (", fit_S0$scale, "; log scale)"),
      y=NULL,
      caption="Estimates from survey-weighted GLM (quasi-Poisson log link for PR when available)."
    ) +
    theme_journal()
}, fatal=TRUE)
save_png650("Figure4_aPR_Forest_Education.png", Fig4, width_in=8.8, height_in=5.8, dpi=650)

# Figure 5 — standardized adjusted prevalence by education and residence (interaction)
Fig5 <- with_step("Figure 5", {
  tbl <- Table5_num %>%
    mutate(
      residence=factor(residence, levels=levels(dat_cc$residence)),
      educ=factor(educ, levels=levels(dat_cc$educ))
    )
  
  ggplot(tbl, aes(x=educ, y=est, color=residence, group=residence)) +
    geom_line(linewidth=1.0) +
    geom_point(size=2.8) +
    geom_errorbar(aes(ymin=lcl, ymax=ucl), width=0.10, linewidth=0.6) +
    scale_color_manual(values=pal_res) +
    scale_y_continuous(labels=percent_format(accuracy=1), limits=c(0,1),
                       expand=expansion(mult=c(0.02,0.06))) +
    labs(
      title="Figure 5. Education gradients in modern contraceptive use by residence (women in union), Malawi DHS 2024",
      subtitle=paste0("Standardized adjusted prevalence from educ×residence model. Global interaction p = ", signif(p_int,3), "."),
      x="Maternal education",
      y="Adjusted prevalence (standardized)",
      color="Residence",
      caption="Interaction model: educ×residence + wealth + region + age (linear)."
    ) +
    theme_journal()
}, fatal=TRUE)
save_png650("Figure5_AdjustedPrev_Education_byResidence.png", Fig5, width_in=8.8, height_in=5.8, dpi=650)

# ---------------------------
# 9) Session info
# ---------------------------
with_step("Write session info", {
  writeLines(capture.output(sessionInfo()), path(OUT_ROOT,"logs","sessionInfo.txt"))
}, fatal=FALSE)

log_line("DONE.")
message("H3 manuscript pack complete. Outputs in: ", OUT_ROOT)