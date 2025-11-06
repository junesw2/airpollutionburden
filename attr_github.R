#' Simulate Attributable Mortality under Different Air Pollution Scenarios
#'
#' This function simulates attributable numbers and fractions of mortality
#' under various exposure scenarios (e.g., policy standards or counterfactuals),
#' using Monte Carlo simulation based on second-stage model results.
#'
#' @details
#' The function aggregates exposures, applies exposure-response functions,
#' and computes attributable mortality numbers and rates across simulations.
#' It supports stratification by sex, cause, and age group, and can output
#' annual estimates or overall averages.
#'
#' @param stage2res Results from the second-stage model (e.g., meta-regression).
#' @param meta_df Metadata associated with the stage 2 results (regions, parameters, etc.).
#' @param DATALIST_PRED1 List of predicted exposure data for all regions and years.
#' @param iSENSITIVITY_expose Exposure sensitivity setting (used to define AP).
#' @param vRNGs Character vector of scenario names.
#'   Defaults to \code{c("Attribute", "EU standard old", "EU standard 2024", "WHO standard 2021")}.
#' @param syear Start year for analysis (default = 2003).
#' @param eyear End year for analysis (default = 2023).
#' @param sSEX Vector of sexes included in the analysis.
#' @param sCAUSE Vector of causes of death included in the analysis.
#' @param sAGE Vector of age groups included in the analysis.
#' @param eachyear Logical; if \code{TRUE}, results are reported for each year separately,
#'   otherwise aggregated (default = TRUE).
#' @param nsim Number of Monte Carlo simulations (default = 1000).
#' @param byrate Scaling factor for excess rate computation (default = 1e7).
#' @param seed Random seed for reproducibility (default = 12345).
#'
#' @return A \code{data.table} containing simulated attributable mortality counts,
#' rates, and fractions by region, year, sex, cause, and age group.
#' 
#' Requires packages: MASS, dlnm, foreach, iterators, itertools, data.table, dplyr.
simulate_attribution <- function(stage2res, meta_df, DATALIST_PRED1, iSENSITIVITY_expose, vRNGs = c("Attribute","EU standard old","EU standard 2024","WHO standard 2021"),syear=2003,eyear=2023,
                                 sSEX=sSEX,sCAUSE=sCAUSE,sAGE=sAGE,eachyear=TRUE,
                                 nsim = 1000,  byrate = 1e6,
                                 seed = 12345) {
  # ---- deps
  require(MASS)
  require(dlnm)
  require(foreach)
  require(iterators)
  require(itertools)
  require(data.table)
  require(dplyr)
  require(mixmeta)
  pop_mort <- readRDS(file = "/PROJECTES/ADAPTATION/proj/zchen/P20230906_ZC_population_weighted_AP/pop_mort_EA_regions_2003_2030_zhao.Rdata")
  
  # ---- basics
  
  na <- nrow(meta_df)
  mortname <- paste0(sSEX, "_", sCAUSE, "_", sAGE, "_mort")
  MIN_LAG_APvar <- 0
  MAX_LAG_T_vectors = c(3,3,3,3,3,3,3,3,3,3,3,3,3,3) #rep(21,14) ###lag of AP
  MAX_LAG_APvar_vectors =  c(2,2,2,2,2,2,3,3,2,2,2,3,3,3) ###lag of AP
  MAX_LAG_APvar_vectors_cardio =  c(2,2,2,2,1,1,3,3,2,2,1,2,2,2) ###lag of AP
  MAX_LAG_APvar_vectors_respir =  c(1,1,1,1,1,1,3,3,1,1,1,2,2,2) ###lag of AP
  MIN_LAG_APvar=0
  if(sCAUSE=="allcause"){#
    MAX_LAG_APvar = MAX_LAG_APvar_vectors[iSENSITIVITY_expose]
  } else if(sCAUSE=="cardio"){#
    MAX_LAG_APvar = MAX_LAG_APvar_vectors_cardio[iSENSITIVITY_expose]
  } else if(sCAUSE=="respir"){#
    MAX_LAG_APvar = MAX_LAG_APvar_vectors_respir[iSENSITIVITY_expose]}
  meta_df_nomiss <- meta_df
  
  numvars <- names(meta_df_nomiss)[sapply(meta_df_nomiss, is.numeric)]
  
  for (v in numvars) {
    nas <- is.na(meta_df_nomiss[[v]])
    if (any(nas)) {
      meta_df_nomiss[[v]][nas] <- mean(meta_df_nomiss[[v]], na.rm = TRUE)
    }
  }
  set.seed(seed)
  # 1) simulate meta-coefficients
  metacoefsim <- MASS::mvrnorm(nsim, coef(stage2res), vcov(stage2res))
  nc=1
  # 2) design matrix (expand for multivariate outcome with Kronecker on identity(nc))
  design_mat <- model.matrix(delete.response(terms(stage2res)), meta_df_nomiss) %x% diag(nc)
  
  # 3) fixed part of city-age coefficients for all sims
  fixsim <- design_mat %*% t(metacoefsim)  # (na*nc) x nsim
  nreg <- dim(fixsim)[1]
  # 4) PREDICTING ONLY THE RANDOM-EFFECT RESIDUALS
  random <- blup(stage2res, type="residual",vcov=T) 
  ## for non linear
  # rancoef <- sapply(random, "[[", "blup") |> t()
  # ranvcov <- sapply(random, function(x) vechMat(x$vcov)) |> t()
  rancoef <- random$blup |> t()
  ranvcov <- random$vcov |> t()
  # 5) Generate random part from
  
  ransim <- foreach(co = iter(rancoef, 'cell'), v = iter(ranvcov, 'cell'), 
                    .combine = rbind) %do% {
                      t(MASS::mvrnorm(nsim, mu = as.numeric(co), Sigma = as.numeric(v)))
                    }
  
  # 6) Add them and rearrange as array
  coefsim <- fixsim + ransim
  
  # ---- per-city attribution over all sims
  # iterator over rows & over sims of coefficient vectors
  row_it <- 1:dim(meta_df)[1]
  simres <- foreach::foreach(i = row_it, .combine = bind_rows) %do% {
    ires <- meta_df[i, , drop = FALSE]
    ipop_mort <- pop_mort[pop_mort$REG==ires$REG, ]
    sim  <- coefsim[i, ]
    attr_reg_df <- DATALIST_PRED1[[ires$REG]] 
    attr_reg_df$cams_pmcoarse <- attr_reg_df$pm10 - attr_reg_df$pm2p5
    
    pollutant <- var_df[iSENSITIVITY_expose]
    attr_reg_df <- attr_reg_df %>% filter(!is.na(attr_reg_df[[pollutant]]) & !is.na(attr_reg_df[[mortname]]))
    if (length(attr_reg_df[[pollutant]]) == 0) return(data.table())
    attr_reg_df[[pollutant]][attr_reg_df[[pollutant]] < 0] <- 0
    
    # Region-specific cutoffs from 2019
    ref2019 <- attr_reg_df[attr_reg_df$year == 2019, pollutant, drop = TRUE]
    L5P_REG  <- quantile(ref2019, 0.99, na.rm = TRUE) * 0.95
    L10P_REG <- quantile(ref2019, 0.99, na.rm = TRUE) * 0.90
    L15P_REG <- quantile(ref2019, 0.99, na.rm = TRUE) * 0.85
    Min_REG  <- min(attr_reg_df[[pollutant]], na.rm = TRUE)
    
    # Standards
    stds <- list(
      pm25 = list(EU_old = NA, EU_2024 = 25, WHO_2005 = 25, WHO_2021 = 15),
      pmcoarse = list(EU_old = NA, EU_2024 = 20, WHO_2005 = 25, WHO_2021 = 30),
      pm10 = list(EU_old = 50, EU_2024 = 45, WHO_2005 = 50, WHO_2021 = 45),
      no2 = list(EU_old = NA, EU_2024 = 50, WHO_2005 = NA, WHO_2021 = 25),
      o3 = list(EU_old = 120, EU_2024 = 100, WHO_2005 = 100, WHO_2021 = 100)
    )
    
    # Map pollutant name
    poll_class <- case_when(
      pollutant %in% c("pm25v0","pm2p5","bscpm25") ~ "pm25",
      pollutant %in% c("pmcoarse","cams_pmcoarse","bsc_pmcoarse") ~ "pmcoarse",
      pollutant %in% c("pm10v0","pm10","bscpm10") ~ "pm10",
      pollutant %in% c("no2pred","no2","bscno2") ~ "no2",
      pollutant %in% c("o3_8hpred","o3_8h") ~ "o3",
      TRUE ~ "other"
    )
    if (poll_class == "other") stop("Invalid pollutant: ", pollutant)
    
    # Basis
    VAR_FUN_APvar <- "lin"
    ONEBASIS_APvar <- onebasis(attr_reg_df[[pollutant]], fun = VAR_FUN_APvar)
    
    # Centering baseline
    cen_default <- list(
      "EU standard old"    = stds[[poll_class]]$EU_old,
      "EU standard 2024"   = stds[[poll_class]]$EU_2024,
      "WHO standard 2005"  = stds[[poll_class]]$WHO_2005,
      "WHO standard 2021"  = stds[[poll_class]]$WHO_2021,
      "Attribute"          = ifelse(poll_class == "o3", 70, 0),
      "5% lower than 2019" = L5P_REG,
      "10% lower than 2019"= L10P_REG,
      "15% lower than 2019"= L15P_REG,
      "Regional Minimum"   = Min_REG
    )
    
    # Center exposures
    bvarcen <- scale(ONEBASIS_APvar, 
                     center = onebasis(cen_default[[vRNGs[1]]], fun = VAR_FUN_APvar), 
                     scale = FALSE)
    
    # Lagged mortality
    LAGGED_MORT_MATRIX <- tsModel::Lag(attr_reg_df[[mortname]], 
                                       -seq(MIN_LAG_APvar, MAX_LAG_APvar))
    LAGGED_MORT_VECTOR <- rowMeans(LAGGED_MORT_MATRIX, na.rm = FALSE)
    
    # Daily AN per sim
    anday_mat <- (1 - exp(- bvarcen %*% sim)) * LAGGED_MORT_VECTOR
    
    # Period definition
    if (!eachyear) {
      vPER <- paste0(syear,"-",eyear)
      pername <- paste0(syear,"_",eyear)
    } else {
      vPER <- paste0(syear:eyear,"-",syear:eyear)
      pername <- "Eachyear"
    }
    sPER <- as.numeric(substr(gsub("-", "", vPER), 1, 4))
    ePER <- as.numeric(substr(gsub("-", "", vPER), 5, 8))
    nPER <- length(vPER)
    
    outlist <- list()
    
    for (iPER in seq_len(nPER)) {
      vTIM <- which(attr_reg_df$year >= sPER[iPER] & attr_reg_df$year <= ePER[iPER])
      dyear <- 365.25
      if (poll_class == "o3") {
        vTIM <- which(attr_reg_df$year >= sPER[iPER] & attr_reg_df$year <= ePER[iPER] &
                        attr_reg_df$temp >= median(attr_reg_df$temp, na.rm = TRUE))
        dyear <- 365.25 / 2
      }
      
      for (rng in vRNGs) {
        thr <- switch(rng,
                      "EU standard old"    = attr_reg_df[[pollutant]][vTIM] > stds[[poll_class]]$EU_old,
                      "EU standard 2024"   = attr_reg_df[[pollutant]][vTIM] > stds[[poll_class]]$EU_2024,
                      "WHO standard 2005"  = attr_reg_df[[pollutant]][vTIM] > stds[[poll_class]]$WHO_2005,
                      "WHO standard 2021"  = attr_reg_df[[pollutant]][vTIM] > stds[[poll_class]]$WHO_2021,
                      "Attribute"          = attr_reg_df[[pollutant]][vTIM] > 0,
                      "5% lower than 2019" = attr_reg_df[[pollutant]][vTIM] > L5P_REG,
                      "10% lower than 2019"= attr_reg_df[[pollutant]][vTIM] > L10P_REG,
                      "15% lower than 2019"= attr_reg_df[[pollutant]][vTIM] > L15P_REG,
                      "Regional Minimum"   = attr_reg_df[[pollutant]][vTIM] > Min_REG,
                      stop("Invalid AP Range: ", rng)
        )
        
        heatind <- vTIM[thr]
        
        outlist[[paste(iPER, rng)]] <- data.table::data.table(
          REG     = ires$REG,
          AP      = pollutant,
          pop     = ipop_mort$pop[ipop_mort$year==sPER[iPER]],
          cou     = ires$NUTS0,
          geozone = ires$region,
          sim     = seq_len(nsim),
          year    = sPER[iPER],
          period  = vPER[iPER],
          type    = rng,
          an_total = colSums(anday_mat[heatind, , drop = FALSE], na.rm = TRUE),
          mort    = sum(attr_reg_df[[mortname]][vTIM], na.rm = TRUE)
        )
      }
    }
    data.table::rbindlist(outlist)
  }
  
  # excess death rates per city-age-sim
  simres <- simres |>
    dplyr::mutate(dplyr::across(dplyr::starts_with("an_"),
                                ~ byrate * .x / pop,
                                .names = "excessrate_{.col}")) |>
    dplyr::mutate(dplyr::across(dplyr::starts_with("an_"),
                                ~  .x / mort *100,
                                .names = "fraction_{.col}")) |>
    dplyr::rename_with(~ gsub("_an_", "_", .x, fixed = TRUE))
  
  # CIs for city-age
  attr_nuts3_df_out <- simres |>
    dplyr::group_by(REG, AP, cou, year, type) |>
    dplyr::summarise(
      dplyr::across(
        c(dplyr::starts_with("an_"),
          dplyr::starts_with("excessrate"),
          dplyr::starts_with("fraction")),
        list(
          value = ~ pmax(stats::quantile(.x, 0.5, na.rm = TRUE), 0),
          low   = ~ stats::quantile(.x, 0.025, na.rm = TRUE),
          high  = ~ stats::quantile(.x, 0.975, na.rm = TRUE)
        )
      ),
      pop  = sum(pop, na.rm = TRUE)/nsim,
      mort = sum(mort, na.rm = TRUE)/nsim,
      .groups = "drop"
    )
  # ---- Region x age aggregation (per sim)
  cousim <- simres[
    , c(
      lapply(.SD, sum, na.rm = TRUE),
      .(pop = sum(pop, na.rm = TRUE),
        mort = sum(mort, na.rm = TRUE))
    ),
    by = .(cou, AP, year, type, sim),
    .SDcols = patterns("^an_")
  ][
    , (paste0("excessrate_", grep("^an_", names(.SD), value = TRUE))) := lapply(.SD, \(x) byrate * x / pop), 
    by = .(cou, AP, year, type, sim),
    .SDcols = patterns("^an_")
  ][
    , (paste0("fraction_", grep("^an_", names(.SD), value = TRUE))) := 
      lapply(.SD, \(x) x / mort * 100), 
    by = .(cou, AP, year, type, sim),
    .SDcols = patterns("^an_")
  ]
  
  attr_cou_df_out <- cousim |>
    dplyr::group_by(cou, AP,  year, type) |>
    dplyr::summarise(
      dplyr::across(dplyr::starts_with("an_") | dplyr::starts_with("excessrate") |  dplyr::starts_with("fraction"),
                    list(value = ~ max(stats::quantile(.x, 0.5, na.rm = TRUE),0),
                         low = ~ stats::quantile(.x, 0.025, na.rm = TRUE),
                         high = ~ abs(stats::quantile(.x, 0.975, na.rm = TRUE)))),
      pop = sum(pop, na.rm = TRUE)/nsim,
      mort = sum(mort, na.rm = TRUE)/nsim,
      .groups = "drop"
    )
  
  # ---- Region  aggregation (per sim)
  regionsim <-  simres[
    , c(
      lapply(.SD, sum, na.rm = TRUE),
      .(pop = sum(pop, na.rm = TRUE),
        mort = sum(mort, na.rm = TRUE))
    ),
    by = .(geozone, AP, year, type, sim),
    .SDcols = patterns("^an_")
  ][
    , (paste0("excessrate_", grep("^an_", names(.SD), value = TRUE))) := lapply(.SD, \(x) byrate * x / pop), 
    by = .(geozone,AP, year, type, sim),
    .SDcols = patterns("^an_")
  ][
    , (paste0("fraction_", grep("^an_", names(.SD), value = TRUE))) := 
      lapply(.SD, \(x) x / mort * 100), 
    by = .(geozone, AP, year, type, sim),
    .SDcols = patterns("^an_")
  ]
  
  
  attr_reg_df_out <- regionsim |>
    dplyr::group_by(geozone, AP,  year, type) |>
    dplyr::summarise(
      dplyr::across(dplyr::starts_with("an_") | dplyr::starts_with("excessrate") |  dplyr::starts_with("fraction"),
                    list(value = ~ max(stats::quantile(.x, 0.5, na.rm = TRUE),0),
                         low = ~ stats::quantile(.x, 0.025, na.rm = TRUE),
                         high = ~ abs(stats::quantile(.x, 0.975, na.rm = TRUE)))),
      pop = sum(pop, na.rm = TRUE)/nsim,
      mort = sum(mort, na.rm = TRUE)/nsim,
      .groups = "drop"
    )
  
  # ---- Region  aggregation (per sim)
  eusim <- simres[
    , c(
      lapply(.SD, sum, na.rm = TRUE),
      .(pop = sum(pop, na.rm = TRUE),
        mort = sum(mort, na.rm = TRUE))
    ),
    by = .( AP, year, type, sim),
    .SDcols = patterns("^an_")
  ][
    , (paste0("excessrate_", grep("^an_", names(.SD), value = TRUE))) := lapply(.SD, \(x) byrate * x / pop), 
    by = .(AP, year, type, sim),
    .SDcols = patterns("^an_")
  ][
    , (paste0("fraction_", grep("^an_", names(.SD), value = TRUE))) := 
      lapply(.SD, \(x) x / mort * 100), 
    by = .( AP, year, type, sim),
    .SDcols = patterns("^an_")
  ]
  
  
  attr_eu_df_out <- eusim |>
    dplyr::group_by(AP,  year, type) |>
    dplyr::summarise(
      dplyr::across(dplyr::starts_with("an_") | dplyr::starts_with("excessrate") |  dplyr::starts_with("fraction"),
                    list(value = ~ max(stats::quantile(.x, 0.5, na.rm = TRUE),0),
                         low = ~ stats::quantile(.x, 0.025, na.rm = TRUE),
                         high = ~ abs(stats::quantile(.x, 0.975, na.rm = TRUE)))),
      pop = sum(pop, na.rm = TRUE)/nsim,
      mort = sum(mort, na.rm = TRUE)/nsim,
      .groups = "drop"
    )
  
  ###save the ATTR results
  foldout = paste0( "/PROJECTES/ADAPTATION/proj/zchen/P20230906_ZC_population_weighted_AP/dataout/stage2/blup/table1/policysave/", sSEX, "_", sAGE, "_", sCAUSE, "_SEN_",polsource,"/update/" );
  if( !file_test( "-d", foldout ) ){ dir.create( foldout, recursive = TRUE ); }
  print("");
  print("= Saving the Data =");
  # Raw Data
  save( simres,
        attr_nuts3_df_out,
        attr_cou_df_out,
        attr_reg_df_out,
        attr_eu_df_out,
        file = paste0( foldout, "Stage2_blupattrl_",pername,"_",syear,"_",eyear,"_loc_cou_reg_eu_prev_diflagmort_",pvar,".Rdata" ) );
  
  
  
}