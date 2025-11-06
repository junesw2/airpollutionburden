############################################################
#  Libraries
############################################################
library(dplyr)
library(tidyr)
library(stringr)
library(purrr)
library(mixmeta)
library(MASS)
library(foreach)
library(Matrix)
library(readr)
##### helper: country → macro-region -----------------------------------------
region_obtain <- function(country) {
  region <- country
  region[country %in% c("LI","BE","AT","CH","CZ","DE","FR","GB","UK","IE","LU","NL")] <- "Western Europe"
  region[country %in% c("BA","BG","HR","RO","SI","SK","HU","TR","AL","RS","ME","MK","PL")] <- "Eastern Europe"
  region[country %in% c("CY","EL","ES","GI","IT","MT","PT","GR")]                     <- "Southern Europe"
  region[country %in% c("FI","DK","IS","NO","SE","EE","LT","LV")]                     <- "Northern Europe"
  region
}

##### constants ---------------------------------------------------------------
var_df      <- c("pm25v0","pm2p5","pm10v0","pm10","no2pred","no2","o3_8hpred","o3_8h",
                 "bscpm25","bscpm10","bscno2","pmcoarse","cams_pmcoarse","bsc_pmcoarse")
expose_var  <- c("pm2.5","cams_pm25","pm10","cams_pm10","no2","cams_no2","o3","cams_o3",
                 "bsc_pm25","bscpm10","bscno2","pmcoarse","cams_pmcoarse","bsc_pmcoarse")
ctrtemp_all <- c("ma.temp","cb.temp")

# lags
MAX_LAG_APvar_vectors        <- c(2,2,2,2,2,2,3,3,2,2,2,3,3,3)
MAX_LAG_APvar_vectors_cardio <- c(1,1,2,2,1,1,2,2,2,1,1,2,2,2)
MAX_LAG_APvar_vectors_respir <- c(1,1,1,1,1,1,2,2,1,1,1,2,2,2)
MAX_LAG_T_vectors            <- rep(3, 14)

VAR_FUN_APvar  <- "lin"
VAR_knot_APvar <- c(0.25, 0.75)
studyperiod    <- "2003-2019"

regappar  <- list()
nl        <- 0

# years table (you seem not to use it later, so I’ll leave it)
periods <- data.frame(
  start = 2003:2013,
  end   = 2009:2019
)
nPER <- nrow(periods)

##### loop over age/sex/cause to merge stage 1 -------------------------------------------------
for (sAGE in c("allage","age000064","age065074","age075084","age085ppp")) {
  for (sSEX in c("allsex","men","women")) {
    for (sCAUSE in c("allcause","cardio","respir")) {
      
      # skip combos you don’t have
      if ((sAGE %in% c("age000064","age065074","age075084","age085ppp")) &&
          (sCAUSE %in% c("cardio","respir"))) {
        next
      }
      
      mortname <- paste0(sSEX, "_", sCAUSE, "_", sAGE, "_mort")
      
      # 2003–2013: CAMS not available
      exposed_range <- c(1,3,5,7,9,10,11,12,14)
      
      for (iSENSITIVITY_expose in exposed_range) {
        
        MAX_LAG_T <- MAX_LAG_T_vectors[iSENSITIVITY_expose]
        if (sCAUSE == "allcause") {
          MAX_LAG_APvar <- MAX_LAG_APvar_vectors[iSENSITIVITY_expose]
        } else if (sCAUSE == "cardio") {
          MAX_LAG_APvar <- MAX_LAG_APvar_vectors_cardio[iSENSITIVITY_expose]
        } else {
          MAX_LAG_APvar <- MAX_LAG_APvar_vectors_respir[iSENSITIVITY_expose]
        }
        
        # temperature control: O3 → cb.temp
        iSENSITIVITY_controltemp <- if (iSENSITIVITY_expose %in% c(7,8)) 2 else 1
        ctrtemp <- ctrtemp_all[iSENSITIVITY_controltemp]
        
        ## ---------- 1. overall (whole-period) stage 1 ----------------------
        fold_overall <- paste0(
          "/PROJECTES/ADAPTATION/proj/zchen/P20230906_ZC_population_weighted_AP/dataout/stage1/",
          ifelse(VAR_FUN_APvar == "lin", "", paste0(VAR_FUN_APvar,"/", paste(VAR_knot_APvar, collapse = "_"))),
          "/", studyperiod, "/",
          sSEX, "_", sAGE, "_", sCAUSE, "_SEN_", expose_var[iSENSITIVITY_expose], "_", ctrtemp, "/"
        )
        
        f_coeff_overall <- paste0(
          fold_overall, "Stage1_coeff_", expose_var[iSENSITIVITY_expose],
          "_SEN_", sSEX, "_", sAGE, "_", sCAUSE,
          "_LAGAP_", MAX_LAG_APvar, "_", ctrtemp, "_LAGTEMP_", MAX_LAG_T, ".Rdata"
        )
        
        if (!file.exists(f_coeff_overall)) {
          # nothing for this combo
          next
        }
        
        COEF_MODEL_APvar_overall <- readRDS(f_coeff_overall)
        VCOV_MODEL_APvar_overall <- readRDS(
          paste0(
            fold_overall, "Stage1_vcov_", expose_var[iSENSITIVITY_expose],
            "_SEN_", sSEX, "_", sAGE, "_", sCAUSE,
            "_LAGAP_", MAX_LAG_APvar, "_", ctrtemp, "_LAGTEMP_", MAX_LAG_T, ".Rdata"
          )
        )
        
        ## ---------- 2. time-varying stage 1 (your “timevarying/…”) --------
        fold_tv <- paste0(
          "/PROJECTES/ADAPTATION/proj/zchen/P20230906_ZC_population_weighted_AP/dataout/stage1/timevarying/",
          sSEX, "_", sAGE, "_", sCAUSE, "_SEN_",
          expose_var[iSENSITIVITY_expose], "_", ctrtemp, "/"
        )
        if (!file.exists(fold_tv)) dir.create(fold_tv, recursive = TRUE)
        
        COEF_MODEL_APvar_tv <- readRDS(
          paste0(
            fold_tv, "Stage1_coeff_", expose_var[iSENSITIVITY_expose],
            "_SEN_", sSEX, "_", sAGE, "_", sCAUSE,
            "_LAGAP_", MAX_LAG_APvar, "_", ctrtemp, "_LAGTEMP_", MAX_LAG_T, ".Rdata"
          )
        )
        VCOV_MODEL_APvar_tv <- readRDS(
          paste0(
            fold_tv, "Stage1_vcov_", expose_var[iSENSITIVITY_expose],
            "_SEN_", sSEX, "_", sAGE, "_", sCAUSE,
            "_LAGAP_", MAX_LAG_APvar, "_", ctrtemp, "_LAGTEMP_", MAX_LAG_T, ".Rdata"
          )
        )
        
        ## ---------- 3. extract regions from both objects -------------------
        ## overall object is usually: matrix (regions x 1) + list of vcovs
        vREG_overall <- rownames(COEF_MODEL_APvar_overall)
        ok_overall   <- sapply(VCOV_MODEL_APvar_overall, length) > 0
        vREG_overall <- vREG_overall[ok_overall]
        COEF_MODEL_APvar_overall <- COEF_MODEL_APvar_overall[vREG_overall, , drop = FALSE]
        VCOV_MODEL_APvar_overall <- VCOV_MODEL_APvar_overall[ok_overall]
        
        ## time-varying object is usually matrix (regions x 2) and a matrix/list of vcovs
        # make sure it is a matrix
        COEF_MODEL_APvar_tv <- as.matrix(COEF_MODEL_APvar_tv)
        vREG_tv             <- rownames(COEF_MODEL_APvar_tv)
        
        # main effect = column 1
        coef_tv_main <- COEF_MODEL_APvar_tv[, 1]
        # interaction = column 2 (if exists)
        coef_tv_int  <- if (ncol(COEF_MODEL_APvar_tv) >= 2) COEF_MODEL_APvar_tv[, 2] else rep(NA_real_, length(vREG_tv))
        
        # vcov: these were saved as list-of-lists per col, so be defensive
        # main vcov
        if (is.matrix(VCOV_MODEL_APvar_tv)) {
          # in case it was saved as matrix-of-lists
          vcov_tv_main <- VCOV_MODEL_APvar_tv[, 1]
          vcov_tv_int  <- if (ncol(VCOV_MODEL_APvar_tv) >= 2) VCOV_MODEL_APvar_tv[, 2] else rep(list(NULL), length(vREG_tv))
        } else {
          # if it is just a list (most likely): same vcov for the main part
          vcov_tv_main <- VCOV_MODEL_APvar_tv
          vcov_tv_int  <- VCOV_MODEL_APvar_tv
        }
        
        # keep only regions where we actually have a vcov
        ok_tv_main <- sapply(vcov_tv_main, length) > 0
        vREG_tv    <- vREG_tv[ok_tv_main]
        coef_tv_main <- coef_tv_main[ok_tv_main]
        coef_tv_int  <- coef_tv_int[ok_tv_main]
        vcov_tv_main <- vcov_tv_main[ok_tv_main]
        vcov_tv_int  <- vcov_tv_int[ok_tv_main]
        
        ## ---------- 4. bind three pieces for this AP/sex/age/cause ---------
        nl <- nl + 1
        regappar[[nl]] <-
          dplyr::bind_rows(
            # time-varying, main
            data.frame(
              REG        = vREG_tv,
              AP         = expose_var[iSENSITIVITY_expose],
              period     = "2011",                       # your label
              sAGE       = sAGE,
              sSEX       = sSEX,
              sCAUSE     = sCAUSE,
              MAX_LAG_T  = MAX_LAG_T,
              MAX_LAG_APvar = MAX_LAG_APvar,
              coef       = as.numeric(coef_tv_main),
              vcov       = unlist(vcov_tv_main)
            ),
            # time-varying interaction
            data.frame(
              REG        = vREG_tv,
              AP         = expose_var[iSENSITIVITY_expose],
              period     = "Time-varying Interaction",
              sAGE       = sAGE,
              sSEX       = sSEX,
              sCAUSE     = sCAUSE,
              MAX_LAG_T  = MAX_LAG_T,
              MAX_LAG_APvar = MAX_LAG_APvar,
              coef       = as.numeric(coef_tv_int),
              vcov       = unlist(vcov_tv_int)
            ),
            # whole-period
            data.frame(
              REG        = vREG_overall,
              AP         = expose_var[iSENSITIVITY_expose],
              period     = "Whole-period",
              sAGE       = sAGE,
              sSEX       = sSEX,
              sCAUSE     = sCAUSE,
              MAX_LAG_T  = MAX_LAG_T,
              MAX_LAG_APvar = MAX_LAG_APvar,
              coef       = as.numeric(COEF_MODEL_APvar_overall[, 1]),
              vcov       = unlist(VCOV_MODEL_APvar_overall)
            )
          )
      } # end exposure
    }
  }
}

# combine everything
regappar_df <- dplyr::bind_rows(regappar)

# add NUTS levels
regappar_df$NUTS2  <- substr(regappar_df$REG, 1, 4)
regappar_df$NUTS1  <- substr(regappar_df$REG, 1, 3)
regappar_df$NUTS0  <- substr(regappar_df$REG, 1, 2)
regappar_df$region <- region_obtain(regappar_df$NUTS0)

# --------------------------------------------------------------------
# join with population metadata (your code below stays the same)
# --------------------------------------------------------------------
load("/PROJECTES/ADAPTATION/proj/zchen/P20230906_ZC_population_weighted_AP/Zhao/population_cou_reg.Rdata")

elder_population <- pop_nuts3 %>%
  mutate(
    F_GE65_per_thousand = (F_Y65.74 + F_Y75.84 + F_Y_GE85) / F_TOTAL * 1000,
    M_GE65_per_thousand = (M_Y65.74 + M_Y75.84 + M_Y_GE85) / M_TOTAL * 1000,
    T_GE65_per_thousand = (T_Y65.74 + T_Y75.84 + T_Y_GE85) / T_TOTAL * 1000,
    F_GE75_per_thousand = (F_Y75.84 + F_Y_GE85) / F_TOTAL * 1000,
    M_GE75_per_thousand = (M_Y75.84 + M_Y_GE85) / M_TOTAL * 1000,
    T_GE75_per_thousand = (T_Y75.84 + T_Y_GE85) / T_TOTAL * 1000,
    F_GE85_per_thousand =  F_Y_GE85 / F_TOTAL * 1000,
    M_GE85_per_thousand =  M_Y_GE85 / M_TOTAL * 1000,
    T_GE85_per_thousand =  T_Y_GE85 / T_TOTAL * 1000
  ) %>%
  dplyr::select(
    geo,
    F_GE65_per_thousand, M_GE65_per_thousand, T_GE65_per_thousand,
    F_GE75_per_thousand, M_GE75_per_thousand, T_GE75_per_thousand,
    F_GE85_per_thousand, M_GE85_per_thousand, T_GE85_per_thousand
  )

# location_meta_df must already be in your environment
regappar_df <- regappar_df %>%
  dplyr::left_join(location_meta_df[, 1:25], by = c("REG" = "location")) %>%
  dplyr::left_join(elder_population,          by = c("REG" = "geo"))

#### spilt the merge dataset
metadf0=regappar_df[regappar_df$period=="Whole-period",]
metadf1=regappar_df[regappar_df$period=="Time-varying Interaction",]


######### stage 2
DATATABLE_DATA <- readRDS(
  file = "/PROJECTES/ADAPTATION/proj/zchen/P20230906_ZC_population_weighted_AP/dataout/DATATABLE_DATA_SEX_CAUSE_AGE_RH_20250923.Rdata"
)

var_df <- c(
  "pm25v0","pm2p5","pm10v0","pm10",
  "no2pred","no2",
  "o3_8hpred","o3_8h",
  "bscpm25","bscpm10","bscno2",
  "pmcoarse","cams_pmcoarse","bsc_pmcoarse"
)

expose_var <- c(
  "pm2.5","cams_pm25","pm10","cams_pm10",
  "no2","cams_no2",
  "o3","cams_o3",
  "bsc_pm25","bscpm10","bscno2",
  "pmcoarse","cams_pmcoarse","bsc_pmcoarse"
)

for (polsource in c("QML", "BSC")) {
  
  exposed_range <- dplyr::case_when(
    polsource == "QML"  ~ c(1, 12, 5, 7, 3),
    polsource == "BSC"  ~ c(9, 14, 11, 7, 10),
    polsource == "CAMS" ~ c(2, 13, 6, 8, 4)
  )
  
  for (sAGE in c("allage", "age000064", "age065074", "age075084", "age085ppp")) {
    for (sSEX in c("allsex", "men", "women")) {
      for (sCAUSE in c("allcause", "cardio", "respir")) {
        
        # skip invalid combos
        if ((sAGE %in% c("age000064", "age065074", "age075084", "age085ppp")) &&
            (sCAUSE %in% c("cardio", "respir"))) {
          next
        }
        
        mortname <- paste0(sSEX, "_", sCAUSE, "_", sAGE, "_mort")
        
        mort_df <- DATATABLE_DATA %>%
          dplyr::filter(year %in% 2003:2019) %>%
          dplyr::group_by(location) %>%
          dplyr::summarise(
            mort = sum(.data[[mortname]], na.rm = TRUE),
            .groups = "drop"
          )
        
        for (exp_ind in 1:5) {
          
          pvar <- expose_var[exposed_range[exp_ind]]
          
          meta_df <- metadf0 %>%
            dplyr::filter(
              AP    == pvar,
              sCAUSE == sCAUSE,
              sSEX   == sSEX,
              sAGE   == sAGE
            ) %>%
            dplyr::left_join(mort_df, by = c("REG" = "location"))
          
          # pick correct demographic covariate
          agecov <- switch(
            sAGE,
            "allage"    = "T_GE65_per_thousand",
            "age000064" = if (sSEX == "men") "M_GE65_per_thousand" else if (sSEX == "women") "F_GE65_per_thousand" else "T_GE65_per_thousand",
            "age065074" = if (sSEX == "men") "M_GE65_per_thousand" else if (sSEX == "women") "F_GE65_per_thousand" else "T_GE65_per_thousand",
            "age075084" = if (sSEX == "men") "M_GE75_per_thousand" else if (sSEX == "women") "F_GE75_per_thousand" else "T_GE75_per_thousand",
            "age085ppp" = if (sSEX == "men") "M_GE85_per_thousand" else if (sSEX == "women") "F_GE85_per_thousand" else "T_GE85_per_thousand",
            "T_GE65_per_thousand" # fallback
          )
          
          # pollutant means to put in the model
          var_df <- c(
            "pm25v0","pm2p5","pm10v0","pm10",
            "no2pred","no2",
            "o3_8hpred","o3_8h",
            "bscpm25","bscpm10","bscno2",
            "pmcoarse","cams_pmcoarse","bsc_pmcoarse"
          )
          
          formulax <- as.formula(
            paste0(
              "coef ~ ",
              var_df[c(1, 12, 5, 7, 3)[exp_ind]], "_mean",
              " + temp_mean + rh_mean + ",
              agecov,
              " + log(GDP_PurchasePower_Hab)"
            )
          )
          
          stage2res <- mixmeta::mixmeta(
            formulax,
            vcov,
            data    = meta_df,
            control = list(showiter = FALSE, igls.inititer = 10, maxiter = 1000),
            method  = "reml",
            random  = ~ 1 | REG
          )
          
          # save
          fold0 <- paste0(
            "/PROJECTES/ADAPTATION/proj/zchen/P20230906_ZC_population_weighted_AP/dataout/stage2/blup/metares/",
            polsource, "/"
          )
          if (!file_test("-d", fold0)) dir.create(fold0, recursive = TRUE)
          
          saveRDS(
            meta_df,
            file = paste0(
              fold0, "meta_df_", sAGE, "_", sSEX, "_", sCAUSE, "_", pvar, ".rds"
            )
          )
          saveRDS(
            stage2res,
            file = paste0(
              fold0, "stage2res_", sAGE, "_", sSEX, "_", sCAUSE, "_", pvar, ".rds"
            )
          )
        } # exp_ind
      }   # sCAUSE
    }     # sSEX
  }       # sAGE
}         # polsource

