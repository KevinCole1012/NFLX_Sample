# TODO: Add comment
# 
###############################################################################

Toolbox.calculatePortfolioExposures <- function(
		data,
		weightVars, 
		exposureVariables, 
		groupVariable
) {	
	if(is.null(groupVariable)) {
		weights <- data[weightVars]
		dataSub <- data %>% 
				dplyr::select_at(exposureVariables)
		
		results <- Toolbox.weightedSumFast(dataSub, weights, TRUE)
	} else {
		if(length(weightVars) != 1) {
			Toolbox.logError("weightVar should have a length of 1 when using groupVariable")
			stop()		
		}
		
		dataT <- Toolbox.pivotWide(
						data %>%
								dplyr::select_at(c(groupVariable, weightVars)) %>%
								dplyr::filter_at(groupVariable, dplyr::any_vars(!is.na(.))) %>%
								dplyr::mutate(id_ = dplyr::row_number()),
						"id_",
						groupVariable, 
						weightVars,
						fillValue = 0
				) %>%
				dplyr::arrange(id_)
		
		groupValues <- as.character(unique(data[[groupVariable]]))
		groupValues <- groupValues[!is.na(groupValues)]
		
		weights <- dataT[groupValues]
		dataSub <- data %>%
				dplyr::select_at(c(groupVariable, exposureVariables)) %>%
				dplyr::filter_at(groupVariable, dplyr::any_vars(!is.na(.))) %>%
				dplyr::select_at(c(exposureVariables)) %>%
				dplyr::mutate(id_ = dplyr::row_number()) %>%
				dplyr::arrange(id_) %>%
				dplyr::select(-id_)
		
		results <- Toolbox.weightedSumFast(dataSub, weights, TRUE) %>%
				tibble::rownames_to_column(groupVariable)
	}
	
	return(results)
}

#' Calculate weight within each quantile of the given variables.
Toolbox.calculatePortfolioQuantileWeights <- function(
		data, 
		weightVariable, 
		quantileVariables,
		numQuantiles
) {
	varMap <- c(
			"weight_" = weightVariable, 
			quantileVariables
	)

	quantiles <- data %>%
			dplyr::ungroup() %>%
			dplyr::select(dplyr::all_of(varMap)) %>%
			dplyr::mutate(
					dplyr::across(
							dplyr::all_of(quantileVariables),
							~ dplyr::ntile(., numQuantiles)
					)
			)

	# NOTE: Obs with missing quantiles are dropped, this is done because SQL doesn't allow a missing (NULL) value for a primary key
	
	quantileWeights <- list()
	
	for(quantileVar in quantileVariables) {		
		quantileWeights[[quantileVar]] <- quantiles %>%
				dplyr::rename_at(quantileVar, ~ "quantile") %>%
				dplyr::filter(is.finite(quantile)) %>%
				dplyr::group_by(quantile) %>%
				dplyr::summarize(weight = sum(weight_, na.rm = TRUE))
	}
	
	results <- dplyr::bind_rows(quantileWeights, .id = "variable_id")
	
	return(results)
}

Toolbox.calculatePortfolioSummaryParams <- function(
	numQuantiles = 5,
	calcExposures = TRUE,
	calcQuantiles = TRUE,
	calcGroupExposures = TRUE,
	calcSideExposures = TRUE,
	calcSideQuantiles = TRUE,
	calcSideGroupExposures = TRUE,
	calcRiskSummary = TRUE,
	calcMCAR = TRUE,
	calcBeta = TRUE,
	calcDiagTC = TRUE,
	calcFullTC = TRUE	
) {
	results <- list(
			numQuantiles = numQuantiles,
			calcExposures = calcExposures,
			calcQuantiles = calcQuantiles,
			calcGroupExposures = calcGroupExposures,
			calcSideExposures = calcSideExposures,
			calcSideQuantiles = calcSideQuantiles,
			calcSideGroupExposures = calcSideGroupExposures,
			calcRiskSummary = calcRiskSummary,
			calcMCAR = calcMCAR,
			calcBeta = calcBeta,
			calcDiagTC = calcDiagTC,
			calcFullTC = calcFullTC	
	)
	
	return(results)
}

#' Calculate summary statistics for a portfolio include exposures, number of names, total risk, etc.
#' @param data data.frame containing security-level exposures and weights
#' @param identifierVariables Variables used to identify securities
#' @param weightVar The variable that has the portfolio weight
#' @param benchmarkVar The variable that has the benchmark weight
#' @param marketVar The variable that has the market weight
#' @param covData List containing covariance matrix for each risk model
#' @param riskFactorData data.frame containing metadata for each risk factor
#' @param groupVariables Variables that will group calculations (i.e. period_dt) 
Toolbox.calculatePortfolioSummary <- function(
		data, 
		covData, 
		riskFactorData, 
		identifierVariables,  
		weightVar, 
		benchmarkVar, 
		marketVar, 
		sharesVar, 
		priceVar, 
		equityVar,
		signalVariables,
		exposureVariables,
		groupVariables,
		groupExposureVariables,
		advVar,
		riskVariablePrefix = NULL,
		quantileVariables = NULL,
		calcParams = NULL
) {	
	if(is.null(calcParams)) {
		calcParams <- Toolbox.calculatePortfolioSummaryParams()
	}

	# Assign variables in the current environment from calcParams
	
	list2env(calcParams, environment())
	
	# TODO: If risk factors are missing then they are silently dropped causing miscalculated risk estimates.  This
	#		condition should be checked, but the challenge is that it would force inclusion of all risk factors 
	#		which can greatly increase size of dataset
	
	Toolbox.logFine("weightVar: %s", weightVar)
	Toolbox.logFine("benchmarkVar: %s", benchmarkVar)
	Toolbox.logFine("marketVar: %s", marketVar)
	Toolbox.logFine("sharesVar: %s", sharesVar)
	Toolbox.logFine("priceVar: %s", priceVar)
	Toolbox.logFine("equityVar: %s", equityVar)
	Toolbox.logFine("signalVariables: %s", signalVariables)
	Toolbox.logFine("exposureVariables: %s", exposureVariables)	
	Toolbox.logFine("groupVariables: %s", groupVariables)	
	Toolbox.logFine("groupExposureVariables: %s", groupExposureVariables)	
	Toolbox.logFine("quantileVariables: %s", quantileVariables)	
	Toolbox.logFine("advVar: %s", advVar)	
	
	results <- list()
	
	# Use only equity securities.  This is mainly a problem with the country_* variables, specifically country_USA 

	varMap <- c(
			"weight_" = weightVar, 
			"equity_" = equityVar
	)
	
	if(!is.null(benchmarkVar)) {
		varMap <- c(varMap, "weight_benchmark_" = benchmarkVar)
	}
	
	if(!is.null(marketVar)) {
		varMap <- c(varMap, "weight_market_" = marketVar)
	}
	
	if(!is.null(sharesVar)) {
		varMap <- c(varMap, "shares_" = sharesVar)
	}
	
	if(!is.null(priceVar)) {
		varMap <- c(varMap, "price_" = priceVar)
	}
	
	if(!is.null(advVar)) {
		varMap <- c(varMap, "adv_" = advVar)
	}
	
	Toolbox.logFine("varMap: %s", varMap)

	data <- data %>%
			dplyr::rename(dplyr::all_of(varMap)) %>%
			dplyr::ungroup()

	dataEquity <- data %>%
			dplyr::filter(equity_) %>%
			dplyr::mutate(
					weight_long_ = pmax(weight_, 0),
					weight_short_ = pmin(weight_, 0)
			)

	if(calcExposures) {
		# TODO: Re-do this after supporting multiple weight variables in Toolbox.calculatePortfolioExposures()
	
		Toolbox.logDebug("Calculate portfolio exposures")
		
		exposureVariablesSub <- dplyr::intersect(exposureVariables, names(dataEquity))
		
		missingExposureVariables <- dplyr::setdiff(exposureVariables, names(dataEquity))
		
		if(length(missingExposureVariables) > 0) {
			Toolbox.logError("missingExposureVariables: %s", missingExposureVariables)
		}
	
		exposures <- Toolbox.calculatePortfolioExposures(
						dataEquity, 
						"weight_", 
						exposureVariablesSub, 
						NULL
				) %>%
				tidyr::gather("variable_id", "managed")
	
		if(calcSideExposures) {
			Toolbox.logFine("Calculate long exposures")
			
			exposuresLong <- Toolbox.calculatePortfolioExposures(
							dataEquity, 
							"weight_long_", 
							exposureVariablesSub, 
							NULL
					) %>%
					tidyr::gather("variable_id", "long")
			
			Toolbox.logFine("Calculate short exposures")
			
			exposuresShort <- Toolbox.calculatePortfolioExposures(
							dataEquity, 
							"weight_short_", 
							exposureVariablesSub, 
							NULL
					) %>%
					tidyr::gather("variable_id", "short")
	
			exposures <- exposures %>%
					dplyr::left_join(exposuresLong, by = "variable_id")  %>%
					dplyr::left_join(exposuresShort, by = "variable_id")
		}
		
		if(!is.null(benchmarkVar)) {
			Toolbox.logFine("Calculate benchmark exposures")
			
			benchmarkExposures <- Toolbox.calculatePortfolioExposures(
							dataEquity, 
							"weight_benchmark_", 
							exposureVariablesSub, 
							NULL
					) %>%
					tidyr::gather("variable_id", "benchmark")
			
			exposures <- exposures %>%
					dplyr::left_join(benchmarkExposures, by = "variable_id") %>%
					dplyr::mutate(active = managed - benchmark)
		}
		
		results$exposures <- exposures
	}

	if(calcGroupExposures && !is.null(groupVariables) && !is.null(groupExposureVariables)) {
		# TODO: Re-do this after supporting multiple weight variables in Toolbox.calculatePortfolioExposures()
		
		Toolbox.logDebug("Calculate group portfolio exposures")
		
		groupExposureVariablesSub <- dplyr::intersect(groupExposureVariables, names(dataEquity))
		
		missingGroupExposureVariables <- dplyr::setdiff(groupExposureVariables, names(dataEquity))
		
		if(length(missingGroupExposureVariables) > 0) {
			Toolbox.logError("missingGroupExposureVariables: %s", missingGroupExposureVariables)
		}
		
		groupExposures <- lapply(
				Toolbox.createIdentityList(groupVariables),
				\(groupVariable) {
					Toolbox.logFine("groupVariable: %s", groupVariable)
					
					exposures <- Toolbox.calculatePortfolioExposures(
									data = dataEquity, 
									weightVars = "weight_", 
									exposureVariables = c(groupExposureVariablesSub, "equity_"), 
									groupVariable = groupVariable
							) %>%
							dplyr::rename(dplyr::all_of(c(group_id = groupVariable))) %>%
							dplyr::group_by(group_id) %>%
							Toolbox.pivotNarrow("variable_id", "managed", groupExposureVariablesSub) %>%
							dplyr::rename(weight_managed = equity_)
				}
			) %>%
			dplyr::bind_rows(.id = "group_type")
	
		if(!is.null(benchmarkVar)) {
			benchmarkGroupExposures <- lapply(
							Toolbox.createIdentityList(groupVariables),
							\(groupVariable) {
								Toolbox.logFine("groupVariable: %s", groupVariable)
								
								exposures <- Toolbox.calculatePortfolioExposures(
												data = dataEquity, 
												weightVars = "weight_benchmark_", 
												exposureVariables = c(groupExposureVariablesSub, "equity_"), 
												groupVariable = groupVariable
										) %>%
										dplyr::rename(dplyr::all_of(c(group_id = groupVariable))) %>%
										dplyr::group_by(group_id) %>%
										Toolbox.pivotNarrow("variable_id", "benchmark", groupExposureVariablesSub) %>%
										dplyr::rename(weight_benchmark = equity_)
							}
					) %>%
					dplyr::bind_rows(.id = "group_type")
			
			groupExposures <- groupExposures %>%
					dplyr::left_join(
							benchmarkGroupExposures, 
							by = c("variable_id", "group_type", "group_id")
					) %>%
					dplyr::mutate(
							active = managed - benchmark,
							weight_active = weight_managed - weight_benchmark
					)
		}
	
		results$groupExposures <- groupExposures
	}
	
	# Quantile analysis

	if(calcQuantiles && !is.null(quantileVariables)) {
		Toolbox.logDebug("Calculate quantile weights")
		
		# TODO: In the following the quantiles are re-calculated for each weight variable.  This 
		#		could be made more efficient by pre-computing quantiles and then getting weights within them.
		
		quantileWeights <- Toolbox.calculatePortfolioQuantileWeights(
						dataEquity, 
						"weight_", 
						quantileVariables,
						numQuantiles
				) %>%
				dplyr::rename(managed = weight)
		
		if(calcSideQuantiles) {
			Toolbox.logFine("Calculate long quantile weights")
			
			quantileWeightsLong <- Toolbox.calculatePortfolioQuantileWeights(
							dataEquity, 
							"weight_long_", 
							quantileVariables,
							numQuantiles
					) %>%
					dplyr::rename(long = weight)
			
			Toolbox.logFine("Calculate short quantile weights")
			
			quantileWeightsShort <- Toolbox.calculatePortfolioQuantileWeights(
							dataEquity, 
							"weight_short_", 
							quantileVariables,
							numQuantiles
					) %>%
					dplyr::rename(short = weight)
			
			quantileWeights <- quantileWeights %>%
					dplyr::left_join(quantileWeightsLong, by = c("variable_id", "quantile"))  %>%
					dplyr::left_join(quantileWeightsShort, by = c("variable_id", "quantile"))
		}
		
		if(!is.null(benchmarkVar)) {
			Toolbox.logFine("Calculate benchmark quantile weights")
			
			quantileWeightsBenchmark <- Toolbox.calculatePortfolioQuantileWeights(
							dataEquity, 
							"weight_benchmark_", 
							quantileVariables,
							numQuantiles
					) %>%
					dplyr::rename(benchmark = weight)
			
			quantileWeights <- quantileWeights %>%
					dplyr::left_join(quantileWeightsBenchmark, by = c("variable_id", "quantile")) %>%
					dplyr::mutate(active = managed - benchmark)
		}
		
		results$quantileWeights <- quantileWeights 		
	}
	
	# Security-level analysis

	Toolbox.logFine("Calculate security-level analysis")

	varMap2 <- c(
			identifierVariables,
			"weight" = "weight_", 
			"price" = "price_", 
			"equity" = "equity_"
	)
	
	if(!is.null(benchmarkVar)) {
		varMap2 <- c(varMap2, "benchmark" = "weight_benchmark_")
	}

	if(!is.null(sharesVar)) {
		varMap2 <- c(varMap2, "shares" = "shares_")
	}
	
	if(!is.null(advVar)) {
		varMap2 <- c(varMap2, "adv" = "adv_")
	}
	
	results$portfolio <- data %>%
			dplyr::ungroup() %>%
			dplyr::select_at(varMap2)

	if(is.null(benchmarkVar)) {
		results$portfolio <- results$portfolio %>%
				dplyr::mutate(
						benchmark = 0,
						active = weight
				)
	} else {
		results$portfolio <- results$portfolio %>%
				dplyr::mutate(active = weight - benchmark)
	}
	
	if(is.null(sharesVar)) {
		results$portfolio <- results$portfolio %>%
				dplyr::mutate(shares = NA_real_)		
	}
	
	if(is.null(priceVar)) {
		results$portfolio <- results$portfolio %>%
				dplyr::mutate(price = NA_real_)		
	}
	
	if(is.null(equityVar)) {
		results$portfolio <- results$portfolio %>%
				dplyr::mutate(equity = TRUE)		
	}
	
	if(is.null(advVar)) {
		results$portfolio <- results$portfolio %>%
				dplyr::mutate(adv_perc = NA_real_)
	} else {
		results$portfolio <- results$portfolio %>%
				dplyr::mutate(
						adv_perc = ifelse(adv == 0, NA_real_, shares / adv),
						adv_perc_abs = abs(adv_perc)
				)
	}
	
	# Calculate summary statistics
	
	Toolbox.logDebug("Calculate summary")
		
	results$summary <- results$portfolio %>%
			dplyr::summarize(
					value =  sum(shares * price, na.rm = TRUE),
					cash = sum(weight * ifelse(equity, 0, 1), na.rm = TRUE),
					long = sum(weight * ifelse(equity & (weight > 0), 1, 0), na.rm = TRUE),
					short = sum(weight * ifelse(equity & (weight < 0), 1, 0), na.rm = TRUE),
					
					n = sum(equity & weight != 0, na.rm = TRUE),
					n_long = sum(equity & (weight > 0), na.rm = TRUE),
					n_short = sum(equity & (weight < 0), na.rm = TRUE),
					
					n_effective = Toolbox.effectiveN(ifelse(equity, weight, 0)),
					n_effective_long = Toolbox.effectiveN(ifelse(equity & weight > 0, weight, 0)),
					n_effective_short = Toolbox.effectiveN(ifelse(equity & weight < 0, weight, 0)),
					
					n_effective_active = Toolbox.effectiveN(ifelse(equity, active, 0)),
					
					max = Toolbox.safeMax(weight * ifelse(equity, 1, 0)),
					min = Toolbox.safeMin(weight * ifelse(equity, 1, 0)),
					
					max_active = Toolbox.safeMax(active * ifelse(equity, 1, 0)),
					min_active = Toolbox.safeMin(active * ifelse(equity, 1, 0)),
					
					adv_perc_abs_max = Toolbox.safeMax(adv_perc_abs),
					adv_perc_abs_median = median(ifelse(adv_perc_abs == 0, NA_real_, adv_perc_abs), na.rm = TRUE)
			)

	# Risk calculations

	Toolbox.logDebug("Calculate risk")
	
	wgtM <- data[["weight_"]]
	wgtEquityM <- dataEquity[["weight_"]]
	
	if(!is.null(benchmarkVar)) {
		benchmarkM <- data[["weight_benchmark_"]]		
		benchmarkEquityM <- dataEquity[["weight_benchmark_"]]		
	} else {
		benchmarkM <- rep(0, length(wgtM))
		benchmarkEquityM <- rep(0, length(wgtEquityM))
	}
	
	if(!is.null(marketVar)) {
		marketM <- data[["weight_market_"]]		
		marketEquityM <- dataEquity[["weight_market_"]]		
	} else {
		marketM <- rep(0, length(wgtM))
		marketEquityM <- rep(0, length(wgtEquityM))
	}

	allRiskModelDecomp <- list()
	allRiskModelMCAR <- list()
	allRiskModelBeta <- list()
	allRiskModelDiagTC <- list()
	allRiskModelFullTC <- list()
	
	if(calcRiskSummary || calcMCAR || calcBeta || calcFullTC || calcDiagTC) {
		for(riskModelID in names(covData)) {
			Toolbox.logDebug("riskModelID: %s", riskModelID)
			
			riskFactorNames <- riskFactorData[[riskModelID]]$name
			riskFactorVars <- paste(riskVariablePrefix, riskModelID, riskFactorNames, sep = "_")
			
			riskExp <- data %>%
					dplyr::ungroup() %>%
					dplyr::select(dplyr::any_of(riskFactorVars))
	
			riskExpEquity <- dataEquity %>%
					dplyr::ungroup() %>%
					dplyr::select(dplyr::any_of(riskFactorVars))
					
			names(riskExp) <- sub(paste0("^", riskVariablePrefix, "_", riskModelID, "_"), "", names(riskExp)) 
			names(riskExpEquity) <- names(riskExp) 
			
			sriskM <- riskExp$srisk			
			sriskEquityM <- riskExpEquity$srisk			
			
			allCurrencyDecomp <- list()
			allCurrencyMCAR <- list()
			allCurrencyBeta <- list()
			allCurrencyFullTC <- list()
			
			for(currency in names(covData[[riskModelID]])) {
				Toolbox.logDebug("currency: %s", currency)
				
				cov <- covData[[riskModelID]][[currency]]
							
				expM <- Toolbox.exposuresAsMatrix(riskExp, cov)	
				expEquityM <- Toolbox.exposuresAsMatrix(riskExpEquity, cov)	
				
				covM <- Toolbox.covarianceAsMatrix(cov)
				
				if(calcRiskSummary) {
					riskDecomp <- list()
		
					Toolbox.logFine("Calculate managed risk")
					
					riskDecomp$managed <- Toolbox.riskDecomp(wgtM, expM, sriskM, covM)    
		
					totalRisk <- riskDecomp$managed$trisk				
					
					if(!is.null(benchmarkVar)) {
						Toolbox.logFine("Calculate benchmark risk")
						
						riskDecomp$benchmark <- Toolbox.riskDecomp(benchmarkM, expM, sriskM, covM)
						
						Toolbox.logFine("Calculate active risk")
						
						riskDecomp$active <- Toolbox.riskDecomp(wgtM - benchmarkM, expM, sriskM, covM)
						
						activeRisk <- riskDecomp$active$trisk				
					} else {
						activeRisk <- totalRisk
					}
		
					if(!is.null(marketVar)) {
						Toolbox.logFine("Calculate market risk")
						
						riskDecomp$market <- Toolbox.riskDecomp(marketM, expM, sriskM, covM)
					}			
					
					allCurrencyDecomp[[currency]] <- dplyr::bind_rows(riskDecomp, .id = "type")
				}
							
				# MCAR / MCTR
	
				if(calcMCAR) {
					mcar <- Toolbox.marginalRisk(wgtM - benchmarkM, expM, sriskM, covM)
					mctr <- Toolbox.marginalRisk(wgtM, expM, sriskM, covM)									
					
					allCurrencyMCAR[[currency]] <- 
							cbind(
									data[identifierVariables],
									data.frame(
											mcar = mcar,
											mctr = mcar,
											mcar_rel = mcar * (wgtM - benchmarkM) / activeRisk,
											mctr_rel = mctr * wgtM / totalRisk
									)
							)
				}
	
				# Predicted beta

				if(calcBeta) {
					# TODO: Add beta decomposed by factor group
	
					allCurrencyBeta[[currency]] <- data.frame(
							managed = Toolbox.portfolioBeta(wgtM, marketM, expM, sriskM, covM),
							long = Toolbox.portfolioBeta(pmax(wgtM, 0), marketM, expM, sriskM, covM),
							short = Toolbox.portfolioBeta(-pmin(wgtM, 0), marketM, expM, sriskM, covM),
							benchmark = Toolbox.portfolioBeta(benchmarkM, marketM, expM, sriskM, covM)
						) %>%
						dplyr::mutate(active = managed - benchmark)
				}

				# Full TC
			
				if(calcFullTC) {
					Toolbox.logDebug("Calculate Full Covariance TC")
								
					allSignalFullTC <- list()
					
					for(signalVariable in signalVariables) {
						Toolbox.logDebug("signalVariable: %s", signalVariable)

						signalEquityM <- dataEquity[[signalVariable]]
						
						tc <- Toolbox.calcFullTC(
								signalEquityM, 
								wgtEquityM - benchmarkEquityM, 
								expEquityM, 
								sriskEquityM, 
								covM
						)
						
						allSignalFullTC[[signalVariable]] <- data.frame(tc = tc) 			
					}
					
					allCurrencyFullTC[[currency]] <- dplyr::bind_rows(allSignalFullTC, .id = "signal_variable_id")
				}
			}
			
			if(calcRiskSummary) {
				allRiskModelDecomp[[riskModelID]] <- dplyr::bind_rows(allCurrencyDecomp, .id = "currency")
			}
			
			if(calcMCAR) {
				allRiskModelMCAR[[riskModelID]] <- dplyr::bind_rows(allCurrencyMCAR, .id = "currency")
			}
			
			if(calcBeta) {
				allRiskModelBeta[[riskModelID]] <- dplyr::bind_rows(allCurrencyBeta, .id = "currency")
			}
			
			if(calcFullTC) {
				allRiskModelFullTC[[riskModelID]] <- dplyr::bind_rows(allCurrencyFullTC, .id = "currency")
			}
			
			# Diagonal TC
			
			if(calcDiagTC) {
				Toolbox.logDebug("Calculate Diagonal TC")
				
				allSignalDiagTC <- list()
				
				for(signalVariable in signalVariables) {
					Toolbox.logDebug("signalVariable: %s", signalVariable)
					
					signalEquityM <- dataEquity[[signalVariable]]
					
					tc <- Toolbox.calcDiagTC(signalEquityM, wgtEquityM - benchmarkEquityM, sriskEquityM)
					
					allSignalDiagTC[[signalVariable]] <- data.frame(tc = tc) 			
				}
		
				allRiskModelDiagTC[[riskModelID]] <- dplyr::bind_rows(allSignalDiagTC, .id = "signal_variable_id")
			}
		}
	}
		
	# Flatten the risk summary data
	
	results$risk <- list()
	
	if(calcRiskSummary) {
		results$risk$summary <- dplyr::bind_rows(allRiskModelDecomp, .id = "risk_model_id")
	}
	
	if(calcMCAR) {
		results$risk$mcar <- dplyr::bind_rows(allRiskModelMCAR, .id = "risk_model_id")
	}
	
	if(calcBeta) {
		results$risk$beta <- dplyr::bind_rows(allRiskModelBeta, .id = "risk_model_id")
	}
	
	if(calcDiagTC || calcFullTC) {
		results$tc <- list()
		
		if(calcDiagTC) {
			results$tc$diagonal <- dplyr::bind_rows(allRiskModelDiagTC, .id = "risk_model_id")
		}
		
		if(calcFullTC) {
			results$tc$full = dplyr::bind_rows(allRiskModelFullTC, .id = "risk_model_id")
		}
	}

	return(results)
}

#' Calculate portfolio-level return adjusted for various costs
#' @param data data.frame containing security-level exposures and weights
#' @param identifierVariables Variables used to identify securities
#' @param weightVar The variable that has the portfolio weight
#' @param benchmarkVar The variable that has the benchmark weight
#' @param tradeVarVar The variable that has the change in weight due to trading
#' @param equityVar The variable that determines whether a security is equity (and subject to cost adjustment calculation)
#' @param groupVariables Variables that will group calculations (i.e. period_dt) 
Toolbox.calculatePortfolioReturnSummary <- function(
		data, 
		identifierVariables,  
		weightVar, 
		benchmarkVar, 
		marketVar, 
		tradeVar, 
		equityVar,  
		returnVariables, 
		transactionsCost, 
		stockLoanCost, 
		leverageCost
) {
	Toolbox.logFine("weightVar: %s", weightVar)
	Toolbox.logFine("benchmarkVar: %s", benchmarkVar)
	Toolbox.logFine("marketVar: %s", marketVar)
	Toolbox.logFine("tradeVar: %s", tradeVar)
	Toolbox.logFine("returnVariables: %s", returnVariables)
	Toolbox.logFine("transactionsCost: %s", transactionsCost)
	Toolbox.logFine("stockLoanCost: %s", stockLoanCost)
	Toolbox.logFine("leverageCost: %s", leverageCost)
	
	Toolbox.logDebug("Calculate portfolio returns")

	varMap <- c("weight_" = weightVar, "equity_" = equityVar)
	
	if(!is.null(tradeVar)) {
		varMap <- c(varMap, "trade_" = tradeVar)
	}

	Toolbox.logFine("varMap: %s", varMap)
	
	dataEquity <- data %>%
			dplyr::rename(dplyr::all_of(varMap)) %>%
			dplyr::filter(equity_) %>%
			dplyr::mutate(
					weight_long_ = pmax(weight_, 0),
					weight_short_ = pmin(weight_, 0)
			)
	
	# TODO: Re-do this after supporting multiple weight variables in Toolbox.calculatePortfolioExposures()

	Toolbox.logFine("Calculate portfolio returns")
	
	portReturns <- Toolbox.calculatePortfolioExposures(
					dataEquity, 
					"weight_", 
					c(returnVariables), 
					NULL
			) %>%
			tidyr::gather("variable_id", "portfolio")

	Toolbox.logFine("Calculate long returns")
	
	longReturns <- Toolbox.calculatePortfolioExposures(
					dataEquity, 
					"weight_long_", 
					c(returnVariables), 
					NULL
			) %>%
			tidyr::gather("variable_id", "long_portfolio")

	Toolbox.logFine("Calculate short returns")
	
	shortReturns <- Toolbox.calculatePortfolioExposures(
					dataEquity, 
					"weight_short_", 
					c(returnVariables), 
					NULL
			) %>%
			tidyr::gather("variable_id", "short_portfolio")
	
	# TODO: Make robust to missing/empty stockLoanCost, leverageCost
	
	Toolbox.logFine("Calculate long/short costs")

	longShortCosts <- dataEquity %>%
			dplyr::select_at(c("weight_", "weight_long_", "weight_short_")) %>%
			dplyr::summarize(					
					long_weight = sum(weight_long_, na.rm = TRUE),								
					short_weight = -sum(weight_short_, na.rm = TRUE),								
					short_cost = short_weight * stockLoanCost,
					leverage = max(long_weight - 1, 0),
					leverage_cost = leverage * leverageCost
			)

	portReturns <- portReturns %>%
			Toolbox.smartLeftJoin(longReturns, by = "variable_id") %>%
			Toolbox.smartLeftJoin(shortReturns, by = "variable_id") %>%
			Toolbox.smartLeftJoin(longShortCosts, by = NULL)		
	
	if(!is.null(transactionsCost) && !is.null(tradeVar)) {
		Toolbox.logFine("Calculate trade costs")
		
		tradeCosts <- dataEquity %>%
				dplyr::select_at("trade_") %>%
				dplyr::summarize(
						trade_weight = sum(abs(trade_), na.rm = TRUE),								
						trade_cost = trade_weight * transactionsCost
				)
		
		portReturns <- portReturns %>%
				Toolbox.smartLeftJoin(tradeCosts, by = NULL) %>%
				dplyr::mutate(total_cost = short_cost + leverage_cost + trade_cost)		
	} 
	else {
		portReturns <- portReturns %>%
				dplyr::mutate(total_cost = short_cost + leverage_cost)
	}
	
	portReturns <- portReturns %>%
			dplyr::mutate(
					portfolio_adj = portfolio - total_cost,
					short_portfolio_adj = short_portfolio - short_cost,
					long_portfolio_adj = long_portfolio - leverage_cost
			)
	
	if(!is.null(benchmarkVar)) {
		Toolbox.logDebug("Calculate benchmark returns")
		
		benchmarkReturns <- Toolbox.calculatePortfolioExposures(
						dataEquity, 
						benchmarkVar, 
						returnVariables, 
						NULL
				) %>%
				tidyr::gather("variable_id", "benchmark")
		
		portReturns <- portReturns %>%
				dplyr::left_join(benchmarkReturns, by = "variable_id") %>%
				dplyr::mutate(
						active = portfolio - benchmark, 
						active_adj = portfolio_adj - benchmark
				)
	} else {
		portReturns <- portReturns %>%
				dplyr::mutate(
						benchmark = 0, 
						active = portfolio,
						active_adj = portfolio_adj - benchmark
				)
	}

	if(!is.null(marketVar)) {
		Toolbox.logDebug("Calculate market returns")
		
		marketReturns <- Toolbox.calculatePortfolioExposures(
						dataEquity, 
						marketVar, 
						returnVariables,
						NULL
				) %>%
				tidyr::gather("variable_id", "market")
		
		portReturns <- portReturns %>%
				dplyr::left_join(marketReturns, by = "variable_id")
	} else {
		portReturns <- portReturns %>%
				dplyr::mutate(
						market = 0
				)
	}
	
	results <- list(
			returns = portReturns
	)
	
	return(results)
}

