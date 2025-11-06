################################################################################
# Calculation of the Location-Specific Associations
################################################################################

cat("\n= Calculation of the Location-Specific Associations =\n\n")

## ---------------------------------------------------------------------------
## 1. Global settings
## ---------------------------------------------------------------------------
# temperature exposure–response
VAR_FUN_T    <- "ns"
VAR_DEG_T    <- 2 ##for bs
# temperature knots (percentiles)
VAR_PRC_T    <- c(10, 75, 90) / 100

# air-pollution exposure–response
VAR_FUN_APvar <- "lin"   # default, can stay "lin" for most AP vars
LAG_FUN_APvar <- "integer"

# AP knots (only used when VAR_FUN_APvar != "lin")
VAR_knot_APvar <- c(0.25, 0.75)



# MIN lags for Temp
MIN_LAG_APvar <- 0
MIN_LAG_T     <- 0
if (MIN_LAG_T < 0) stop("ERROR: Invalid MIN_LAG_T")

# lag knots for temperature if not using integer
nlagknot    <- 3
LAG_KNOTS_T <- logknots(c(MIN_LAG_T, MAX_LAG_T), nlagknot)

# percentiles for predictions
PRED_PRC_T <- sort(unique(c(
  seq(0.0, 1.0, 0.1),
  seq(1.5, 5.0, 0.5),
  seq(6.0, 94.0, 1.0),
  seq(95.0, 98.5, 0.5),
  seq(99.0, 100.0, 0.1)
) / 100))
if (any(PRED_PRC_T < 0 | PRED_PRC_T > 1)) {
  stop("ERROR: Invalid percentile vector for predictions.")
}

# seasonality
DF_SEAS <- 7

# options: which AP and which control for temp
var_df <- c(
  "pm25v0", "pm2p5", "pm10v0", "pm10",
  "no2pred", "no2",
  "o3_8hpred", "o3",
  "bscpm25", "bscpm10", "bscno2",
  "pmcoarse", "cams_pmcoarse", "bsc_pmcoarse"
)
expose_var <- c(
  "pm2.5", "cams_pm25", "pm10", "cams_pm10",
  "no2", "cams_no2",
  "o3", "cams_o3",
  "bsc_pm25", "bscpm10", "bscno2",
  "pmcoarse", "cams_pmcoarse", "bsc_pmcoarse"
)

# temp control: 1 = MA temp, 2 = cb.temp
ctrtemp <- c("ma.temp", "cb.temp")

# max lags by pollutant & cause
MAX_LAG_APvar_vectors <-  c(2,2,2,2,2,2,3,3,2,2,2,3,3,3) ###lag of AP
MAX_LAG_APvar_vectors_cardio <-  c(1,1,2,2,1,1,2,2,2,1,1,2,2,2) ###lag of AP
MAX_LAG_APvar_vectors_respir <-  c(1,1,1,1,1,1,2,2,1,1,1,2,2,2) ###lag of AP
MAX_LAG_T_vectors <- c(3,3,3,3,3,3,3,3,3,3,3,3,3,3) #rep(21,14) ###lag of AP

## helper: build output folder
build_foldout <- function(studyperiod, sSEX, sAGE, sCAUSE,
                          expose_name, tempctrl,
                          var_fun = VAR_FUN_APvar,
                          var_knots = VAR_knot_APvar) {
  base <- "/PROJECTES/ADAPTATION/proj/zchen/P20230906_ZC_population_weighted_AP/dataout/stage1/"
  # for linear we don’t add subfolder with knots
  sub  <- if (var_fun == "lin") "" else paste0(var_fun, "/", paste(var_knots, collapse = "_"), "/")
  paste0(
    base,
    sub,
    studyperiod, "/",
    sSEX, "_", sAGE, "_", sCAUSE,
    "_SEN_", expose_name, "_", tempctrl, "/"
  )
}

## ---------------------------------------------------------------------------
## 2. Loop over scenarios
## ---------------------------------------------------------------------------
for (studyperiod in c("2003-2019")) {
  
  for (sAGE in c("allage", "age000064", "age065074", "age075084", "age085ppp")) {
    
    for (sSEX in c("allsex", "women", "men")) {
      
      for (sCAUSE in c("allcause", "cardio", "respir")) {
        
        # skip unavailable combos
        if ((sAGE %in% c("age000064","age065074","age075084","age085ppp")) &&
            (sCAUSE %in% c("cardio","respir"))) {
          next
        }
        
        mortname <- paste0(sSEX, "_", sCAUSE, "_", sAGE, "_mort")
        
        # cams not available early
        if (studyperiod %in% c("2003-2012", "2003-2019")) {
          exposed_range <- c(1,3,5,7,9,10,11,12,14)
        } else {
          exposed_range <- 1:14
        }
        
        for (iSENSITIVITY_expose in exposed_range) {
          
          # set max lags for this pollutant and cause
          MAX_LAG_T <- MAX_LAG_T_vectors[iSENSITIVITY_expose]
          if (sCAUSE == "allcause") {
            MAX_LAG_APvar <- MAX_LAG_APvar_vectors[iSENSITIVITY_expose]
          } else if (sCAUSE == "cardio") {
            MAX_LAG_APvar <- MAX_LAG_APvar_vectors_cardio[iSENSITIVITY_expose]
          } else {
            MAX_LAG_APvar <- MAX_LAG_APvar_vectors_respir[iSENSITIVITY_expose]
          }
          ##ozone use crossbasis for temperature as previous study
          if(iSENSITIVITY_expose %in% c(7,8)){iSENSITIVITY_controltemp = 2}  
          # temp lag function depends on choice
          if ((iSENSITIVITY_controltemp == 2) && (MAX_LAG_T <= 4)) {
            LAG_FUN_T <- "integer"
          } else {
            LAG_FUN_T <- "ns"
          }
          
          # AP lag function depends on length
          if (MAX_LAG_APvar <= 7) {
            LAG_FUN_APvar <- "integer"
          } else {
            LAG_FUN_APvar <- "ns"
            LAG_KNOTS_APvar <- logknots(c(MIN_LAG_APvar, MAX_LAG_APvar), nlagknot)
          }
          
          # check if we already ran this
          expose_name <- expose_var[iSENSITIVITY_expose]
          foldout     <- build_foldout(
            studyperiod,
            sSEX, sAGE, sCAUSE,
            expose_name,
            ctrtemp[iSENSITIVITY_controltemp]
          )
          
          out_file <- paste0(
            foldout,
            "Stage1_coeff_", expose_name,
            "_SEN_", sSEX, "_", sAGE, "_", sCAUSE,
            "_LAGAP_", MAX_LAG_APvar, "_",
            ctrtemp[iSENSITIVITY_controltemp],
            "_LAGTEMP_", MAX_LAG_T, ".Rdata"
          )
          if (file.exists(out_file)) next
          
          if (!file_test("-d", foldout)) dir.create(foldout, recursive = TRUE)
          
          ## -----------------------------------------------------------------
          ## 3. Allocate storage
          ## -----------------------------------------------------------------
          if (VAR_FUN_T == "ns") {
            COEF_MODEL_T <- matrix(NA, nREG, length(VAR_PRC_T) + 1, dimnames = list(vREG))
          } else if (VAR_FUN_T == "bs") {
            COEF_MODEL_T <- matrix(NA, nREG, length(VAR_PRC_T) + VAR_DEG_T, dimnames = list(vREG))
          } else {
            stop("Invalid VAR_FUN_T")
          }
          
          if (VAR_FUN_APvar == "lin") {
            COEF_MODEL_APvar <- matrix(NA, nREG, 1, dimnames = list(vREG))
          } else if (VAR_FUN_APvar == "ns") {
            COEF_MODEL_APvar <- matrix(NA, nREG, length(VAR_knot_APvar) + 1, dimnames = list(vREG))
          } else if (VAR_FUN_APvar == "bs") {
            COEF_MODEL_APvar <- matrix(NA, nREG, length(VAR_knot_APvar) + VAR_DEG_APvar, dimnames = list(vREG))
          } else {
            stop("Invalid VAR_FUN_APvar")
          }
          
          # lagged coef containers
          if (VAR_FUN_APvar == "lin") {
            if (LAG_FUN_APvar == "integer") {
              COEF_MODEL_APvar_lag <- matrix(NA, nREG, MAX_LAG_APvar + 1, dimnames = list(vREG))
            } else {
              COEF_MODEL_APvar_lag <- matrix(NA, nREG, nlagknot + 2, dimnames = list(vREG))
            }
            VCOV_MODEL_APvar_lag <- vector("list", nREG)
            names(VCOV_MODEL_APvar_lag) <- vREG
          } else {
            COEF_MODEL_APvar_lag <- vector("list", length(REF_PRC_AP))
            VCOV_MODEL_APvar_lag <- vector("list", length(REF_PRC_AP))
            names(COEF_MODEL_APvar_lag) <- sprintf("P%03g", 100 * REF_PRC_AP)
            names(VCOV_MODEL_APvar_lag) <- sprintf("P%03g", 100 * REF_PRC_AP)
            for (iP in seq_along(REF_PRC_AP)) {
              if (LAG_FUN_APvar == "integer") {
                COEF_MODEL_APvar_lag[[iP]] <- matrix(NA, nREG, MAX_LAG_APvar + 1, dimnames = list(vREG))
              } else {
                COEF_MODEL_APvar_lag[[iP]] <- matrix(NA, nREG, nlagknot + 2, dimnames = list(vREG))
              }
              VCOV_MODEL_APvar_lag[[iP]] <- vector("list", nREG)
              names(VCOV_MODEL_APvar_lag[[iP]]) <- vREG
            }
          }
          
          VCOV_MODEL_T             <- vector("list", nREG); names(VCOV_MODEL_T) <- vREG
          VCOV_MODEL_APvar         <- vector("list", nREG); names(VCOV_MODEL_APvar) <- vREG
          CROSS_PRED_REG_NOMETA_APvar <- vector("list", nREG); names(CROSS_PRED_REG_NOMETA_APvar) <- vREG
          vQAIC <- array(NA, dim = c(nREG), dimnames = list(vREG))
          GCV   <- vQAIC
          
          if (iSENSITIVITY_controltemp == 1) {
            cat("Adjusted ", 0, "-", MAX_LAG_T, " moving average Temperature\n")
          } else {
            cat("Adjusted Cross basis of 0 -", MAX_LAG_T, "Temperature\n")
          }
          cat("with", toupper(expose_name), "predictions\n")
          
          ## -----------------------------------------------------------------
          ## 4. Loop over regions
          ## -----------------------------------------------------------------
          for (iREG in seq_len(nREG)) {
            
            DATALIST_CALI[[iREG]]$dow <- weekdays(DATALIST_CALI[[iREG]]$date)
            
            # skip regions with no mortality or no AP
            if (sum(DATALIST_CALI[[iREG]][[mortname]], na.rm = TRUE) == 0) next
            if (sum(DATALIST_CALI[[iREG]][[var_df[iSENSITIVITY_expose]]], na.rm = TRUE) == 0) next
            
            cat("  Region", iREG, "/", nREG, ":", vREG_NAME[iREG], "(", vREG[iREG], ")\n")
            
            # subset years for this study period
            years   <- strsplit(studyperiod, "-")[[1]]
            yr_from <- as.integer(years[1])
            yr_to   <- as.integer(years[2])
            
            df_m <- DATALIST_CALI[[iREG]][DATALIST_CALI[[iREG]]$year %in% yr_from:yr_to, ]
            if (length(df_m$temp) < MAX_LAG_T) next
            
            # base seasonal formula
            FORMULA_SEA <- as.formula(
              paste0(mortname, " ~ ns(date, df = round(DF_SEAS * length(date) / 365.25))")
            )
            
            # temperature control
            if (iSENSITIVITY_controltemp == 1) {
              # moving average temp
              df_m$mtmean <- tsModel::runMean(df_m$temp, lags = 0:MAX_LAG_T)
              FORMULA_CRB <- as.formula(
                paste0(mortname,
                       " ~ ns(date, df = round(DF_SEAS * length(date) / 365.25))",
                       " + factor(dow) + CROSS_BASIS_APvar",
                       " + ns(mtmean, df = 6) + ns(rh, df = 3)")
              )
            } else {
              # cross-basis temp
              if (VAR_FUN_T == "ns" && LAG_FUN_T == "ns") {
                CROSS_BASIS_TEMP <- crossbasis(
                  df_m$temp, c(MIN_LAG_T, MAX_LAG_T),
                  argvar = list(
                    fun = VAR_FUN_T,
                    knots = quantile(df_m$temp, VAR_PRC_T, na.rm = TRUE),
                    Boundary.knots = range(df_m$temp, na.rm = TRUE)
                  ),
                  arglag = list(knots = LAG_KNOTS_T)
                )
              } else if (VAR_FUN_T == "ns" && LAG_FUN_T == "integer") {
                CROSS_BASIS_TEMP <- crossbasis(
                  df_m$temp, c(MIN_LAG_T, MAX_LAG_T),
                  argvar = list(
                    fun = VAR_FUN_T,
                    knots = quantile(df_m$temp, VAR_PRC_T, na.rm = TRUE),
                    Boundary.knots = range(df_m$temp, na.rm = TRUE)
                  ),
                  arglag = list(fun = "integer")
                )
              } else if (VAR_FUN_T == "bs") {
                CROSS_BASIS_TEMP <- crossbasis(
                  df_m$temp, c(MIN_LAG_T, MAX_LAG_T),
                  argvar = list(
                    fun = VAR_FUN_T,
                    degree = VAR_DEG_T,
                    knots = quantile(df_m$temp, VAR_PRC_T, na.rm = TRUE)
                  ),
                  arglag = list(knots = LAG_KNOTS_T)
                )
              } else {
                stop("Invalid temp VAR_FUN")
              }
              
              FORMULA_CRB <- as.formula(
                paste0(mortname,
                       " ~ ns(date, df = round(DF_SEAS * length(date) / 365.25))",
                       " + factor(dow) + CROSS_BASIS_APvar + CROSS_BASIS_TEMP",
                       " + ns(rh, df = 3)")
              )
            }
            
            # AP cross-basis
            ap_vec <- df_m[[var_df[iSENSITIVITY_expose]]]
            if (VAR_FUN_APvar == "lin" && LAG_FUN_APvar == "integer") {
              CROSS_BASIS_APvar <- crossbasis(
                ap_vec, c(MIN_LAG_APvar, MAX_LAG_APvar),
                argvar = list(fun = "lin"),
                arglag = list(fun = "integer")
              )
            } else if (VAR_FUN_APvar == "lin" && LAG_FUN_APvar == "ns") {
              CROSS_BASIS_APvar <- crossbasis(
                ap_vec, c(MIN_LAG_APvar, MAX_LAG_APvar),
                argvar = list(fun = "lin"),
                arglag = list(knots = LAG_KNOTS_APvar)
              )
            } else if (VAR_FUN_APvar == "ns" && LAG_FUN_APvar == "integer") {
              CROSS_BASIS_APvar <- crossbasis(
                ap_vec, c(MIN_LAG_APvar, MAX_LAG_APvar),
                argvar = list(
                  fun = "ns",
                  knots = quantile(ap_vec, VAR_knot_APvar, na.rm = TRUE),
                  Boundary.knots = range(ap_vec, na.rm = TRUE)
                ),
                arglag = list(fun = "integer")
              )
            } else if (VAR_FUN_APvar == "ns" && LAG_FUN_APvar == "ns") {
              CROSS_BASIS_APvar <- crossbasis(
                ap_vec, c(MIN_LAG_APvar, MAX_LAG_APvar),
                argvar = list(
                  fun = "ns",
                  knots = quantile(ap_vec, VAR_knot_APvar, na.rm = TRUE),
                  Boundary.knots = range(ap_vec, na.rm = TRUE)
                ),
                arglag = list(knots = LAG_KNOTS_APvar)
              )
            } else {
              stop("Invalid AP var/lag combination")
            }
            
            # fit model
            next_flag <- FALSE
            tryCatch({
              GLM_MODEL_CRB <- glm(
                formula = FORMULA_CRB,
                data    = df_m,
                family  = quasipoisson,
                na.action = "na.exclude"
              )
            }, error = function(e) {
              cat("ERROR:", conditionMessage(e), "\n")
              next_flag <<- TRUE
            })
            if (next_flag) next
            if (!GLM_MODEL_CRB$converged) next
            
            df_m$mort_pred_seas <- predict(GLM_MODEL_CRB, type = "response")
            vQAIC[iREG] <- QAIC(GLM_MODEL_CRB)
            GCV[iREG]   <- tryCatch(gcv_f(GLM_MODEL_CRB), error = function(e) NA)
            
            # crosspred AP (centred at 0)
            if (VAR_FUN_APvar == "lin") {
              CROSS_PRED_REG_NOMETA_APvar[[iREG]] <- suppressMessages(
                crosspred(CROSS_BASIS_APvar, GLM_MODEL_CRB,
                          model.link = "log", at = 10, bylag = 1, cen = 0)
              )
            } else {
              CROSS_PRED_REG_NOMETA_APvar[[iREG]] <- suppressMessages(
                crosspred(CROSS_BASIS_APvar, GLM_MODEL_CRB,
                          model.link = "log",
                          at = quantile(ap_vec, PRED_PRC_T, na.rm = TRUE),
                          bylag = 1, cen = 0)
              )
            }
            
            # reduced AP
            REDUCED_APvar <- crossreduce(
              CROSS_BASIS_APvar,
              GLM_MODEL_CRB,
              model.link = "log",
              cen = 0
            )
            COEF_MODEL_APvar[iREG, ] <- coef(REDUCED_APvar)
            VCOV_MODEL_APvar[[iREG]] <- vcov(REDUCED_APvar)
            
            # lag-specific reductions
            if (VAR_FUN_APvar == "lin") {
              REDUCED <- crossreduce(
                CROSS_BASIS_APvar,
                GLM_MODEL_CRB,
                type  = "var",
                value = 10,
                cen   = 0
              )
              COEF_MODEL_APvar_lag[iREG, ] <- coef(REDUCED)
              VCOV_MODEL_APvar_lag[[iREG]] <- vcov(REDUCED)
            } else {
              for (iPERC in seq_along(REF_PRC_AP)) {
                REDUCED <- crossreduce(
                  CROSS_BASIS_APvar,
                  GLM_MODEL_CRB,
                  type  = "var",
                  value = quantile(ap_vec, REF_PRC_AP[iPERC], na.rm = TRUE),
                  cen   = 0
                )
                COEF_MODEL_APvar_lag[[iPERC]][iREG, ] <- coef(REDUCED)
                VCOV_MODEL_APvar_lag[[iPERC]][[iREG]] <- vcov(REDUCED)
              }
            }
          } # end region loop
          
          cat("QAIC total:", round(sum(vQAIC, na.rm = TRUE)), "\n")
          cat("GCV  total:", round(sum(GCV,   na.rm = TRUE)), "\n")
          
          # save
          saveRDS(COEF_MODEL_APvar,        paste0(foldout, "Stage1_coeff_", expose_name, "_SEN_", sSEX, "_", sAGE, "_", sCAUSE, "_LAGAP_", MAX_LAG_APvar, "_", ctrtemp[iSENSITIVITY_controltemp], "_LAGTEMP_", MAX_LAG_T, ".Rdata"))
          saveRDS(VCOV_MODEL_APvar,        paste0(foldout, "Stage1_vcov_",  expose_name, "_SEN_", sSEX, "_", sAGE, "_", sCAUSE, "_LAGAP_", MAX_LAG_APvar, "_", ctrtemp[iSENSITIVITY_controltemp], "_LAGTEMP_", MAX_LAG_T, ".Rdata"))
          saveRDS(COEF_MODEL_APvar_lag,    paste0(foldout, "Stage1_lag_coeff_", expose_name, "_SEN_", sSEX, "_", sAGE, "_", sCAUSE, "_LAGAP_", MAX_LAG_APvar, "_", ctrtemp[iSENSITIVITY_controltemp], "_LAGTEMP_", MAX_LAG_T, ".Rdata"))
          saveRDS(VCOV_MODEL_APvar_lag,    paste0(foldout, "Stage1_lag_vcov_",  expose_name, "_SEN_", sSEX, "_", sAGE, "_", sCAUSE, "_LAGAP_", MAX_LAG_APvar, "_", ctrtemp[iSENSITIVITY_controltemp], "_LAGTEMP_", MAX_LAG_T, ".Rdata"))
          saveRDS(CROSS_PRED_REG_NOMETA_APvar,
                  paste0(foldout, "Stage1_crosspred_premeta_", expose_name, "_SEN_", sSEX, "_", sAGE, "_", sCAUSE, "_LAGAP_", MAX_LAG_APvar, "_", ctrtemp[iSENSITIVITY_controltemp], "_LAGTEMP_", MAX_LAG_T, ".Rdata"))
          
          TABLE_QAIC <- data.frame(
            NUTS = c("Europe", seq_len(nREG)),
            name = c("Europe", vREG_NAME),
            code = c("Europe", vREG),
            qaic = sprintf("%.0f", round(c(sum(vQAIC, na.rm = TRUE), vQAIC), 0)),
            gcv  = sprintf("%.0f", round(c(sum(GCV,   na.rm = TRUE), GCV), 0))
          )
          write.csv(
            TABLE_QAIC,
            paste0(foldout, "TABLE_QAIC_REG_", sSEX, "_", sAGE, "_", sCAUSE, "_LAGAP_", MAX_LAG_APvar, "_", ctrtemp[iSENSITIVITY_controltemp], "_LAGTEMP_", MAX_LAG_T, ".csv"),
            row.names = FALSE
          )
          
        } # end pollutant loop
      }     # end cause
    }       # end sex
  }         # end age
}           # end studyperiod
