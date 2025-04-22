# TODO: Add comment
# 
###############################################################################

suppressWarnings(suppressMessages({								
					library(furrr, warn.conflicts = FALSE, quietly = TRUE)
					library(progressr, warn.conflicts = FALSE, quietly = TRUE)
				}))

Toolbox.histBetasCalculateIndiv <- function(
		data, 
		groupVariable,
		returnVariable, 
		indepVariables,
		baseIndepVariables,
		numPeriods,
		minimumObservations,
		calcSplitBetas = FALSE,
		useFast = TRUE
) {
	allIndepVariables <- c(indepVariables, baseIndepVariables)
	
	# Subset to ensure that all variable are finite (non-NA/Inf/NaN)
	
	dataSub <- data %>%
			dplyr::filter_at(
					dplyr::all_of(c(returnVariable, allIndepVariables)), 
					dplyr::all_vars(is.finite(.))
			)
	
	if(nrow(dataSub) > max(
			minimumObservations, 
			1 + length(allIndepVariables),
			numPeriods
	)) {
		dataMat <- as.matrix(dataSub[allIndepVariables])

		if(calcSplitBetas) {
			x <- cbind(1, dataMat, pmax(dataMat, 0, na.rm = TRUE))
			colnames(x) <- c('(Intercept)', allIndepVariables, paste0(allIndepVariables, "_split_"))
		} else {
			x <- cbind(1, dataMat)
			colnames(x) <- c('(Intercept)', allIndepVariables)			
		}
		
		y <- as.matrix(dataSub[returnVariable])
				
		Toolbox.logFine("Calculate rolling regression")

		# TODO: Implement exponentially weighted regression here
		
		fit <- Toolbox.calcRollingRegressionFast(
				x, 
				y, 
				numPeriods, 
				fullStats = TRUE,
				groups = dataSub[[groupVariable]], 
				rollingWindow = TRUE,
				minimumObservations = minimumObservations,
				partial = FALSE,
				useFast = useFast
		)
		
		if(!is.null(fit)) {
			Toolbox.logFine("Calculate betas")

			Toolbox.logFine("Construct data.frame")
			
			results <- data.frame(
					sigma = fit$sigmas, 
					r_squared = fit$r.squareds
			)
			
			results[[groupVariable]] <- fit$groups
			
			# NOTE: If there is only one valid group then the coefs are returned as a vector, not a matrix
			
			if(is.matrix(fit$coefs)) {
				if(calcSplitBetas) {
					if(length(indepVariables) == 1) {						
						results$beta_neg <- fit$coefs[, indepVariables]
						results$beta_pos <- fit$coefs[, paste0(indepVariables, "_split_")] + results$beta_neg						
					} else {
						results$beta_neg <- rowSums(fit$coefs[, indepVariables])
						results$beta_pos <- rowSums(fit$coefs[, paste0(indepVariables, "_split_")]) + results$beta_neg
					}
					
					if(!is.null(baseIndepVariables)) {
						if(length(baseIndepVariables) == 1) {						
							results$beta_base_neg <- fit$coefs[, baseIndepVariables]
							results$beta_base_pos <- fit$coefs[, paste0(baseIndepVariables, "_split_")] + results$beta_base_neg							
						} else {
							results$beta_base_neg <- rowSums(fit$coefs[, baseIndepVariables])
							results$beta_base_pos <- rowSums(fit$coefs[, paste0(baseIndepVariables, "_split_")]) + results$beta_base_neg
						}
					}
				} else {
					if(length(indepVariables) == 1) {
						results$beta <- fit$coefs[, indepVariables]
					} else {
						results$beta <- rowSums(fit$coefs[, indepVariables])						
					}
					
					if(!is.null(baseIndepVariables)) {
						if(length(baseIndepVariables) == 1) {
							results$beta_base <- fit$coefs[, baseIndepVariables]
						}
						else {
							results$beta_base <- rowSums(fit$coefs[, baseIndepVariables])							
						}
							
					}					
				}
			} else {
				if(calcSplitBetas) {
					results$beta_neg <- sum(fit$coefs[indepVariables])
					results$beta_pos <- sum(fit$coefs[paste0(indepVariables, "_split_")]) + results$beta_neg
					
					if(!is.null(baseIndepVariables)) {
						results$beta_base_neg <- sum(fit$coefs[baseIndepVariables])
						results$beta_base_pos <- sum(fit$coefs[paste0(baseIndepVariables, "_split_")]) + results$beta_base_neg
					}
				} else {
					results$beta <- sum(fit$coefs[indepVariables])
				
					if(!is.null(baseIndepVariables)) {
						results$beta_base <- sum(fit$coefs[baseIndepVariables])
					}
				}
			}			
		} else {
			Toolbox.logFine("Null fit")
			
			results <- NULL			
		}
	} else {
		Toolbox.logFine("Insufficient observations")
	
		results <- NULL
	}
		
	return(results)
}

Toolbox.histBetasIndepVars <- function(returnVariable, numLeads, numLags) {
	varNames <- returnVariable
	
	if(numLeads > 0) {
		varNames <- c(
				varNames,
				paste0(returnVariable, paste0("_lead_", seq(numLeads)), "_")
		)
	}
	
	if(numLags > 0) {
		varNames <- c(
				varNames,
				paste0(returnVariable, paste0("_lag_", seq(numLags)), "_")
		)
	}
		
	return(varNames)
}

Toolbox.histBetasCalculate <- function(
		data, 
		dateVariable,
		returnVariable, 
		indexReturnVariables, 
		indexReturnResidVariables, 
		baseIndexReturnVariable, 
		numPeriods,
		minimumObservations,
		periodVariable = NULL,
		periodicity = "day",
		dateOffset = 0,
		winsorizeReturns = c(.01, .99),
		calcSplitBetas = FALSE,
		numLeads = numLeads,
		numLags = numLags,
		useFast = TRUE,
		progressor = NULL
) {
	# Eliminate all of the initial rows that have zero returns because Toolbox.calcRollingRegressionFast() doesn't play well with them

	firstNonZeroReturn <- which(data[[returnVariable]] != 0)[1]
	
	if(is.na(firstNonZeroReturn)) {
		return(NULL)
	} else if(firstNonZeroReturn > 1) {
		data <- data %>%
				dplyr::slice(- seq(1:(firstNonZeroReturn - 1)))		
	}	
	
	if(!is.null(winsorizeReturns)) {
		data <- data %>%
				dplyr::mutate_at(dplyr::all_of(returnVariable), ~ Toolbox.winsorize(., winsorizeReturns[1], winsorizeReturns[2]))
	}
	
	if(is.null(periodVariable)) {		
		data <- data %>%
				dplyr::arrange(dplyr::all_of(dateVariable)) %>%
				dplyr::mutate(
						group_date_ = factor(lubridate::floor_date(eval(as.name(dateVariable)) + numLeads + dateOffset, periodicity)),
						group_ = as.integer(group_date_)
				)
		
		groupMapping <- data %>%
				dplyr::distinct(group_date_, group_)
		
		groupVariable <- "group_"
	} else {
		groupVariable <- periodVariable
	}

	# Non-residual betas

	Toolbox.logFine("Estimate non-residual betas")
	
	allBetas <- list()
	allSplitBetas <- list()
	
	for(i in seq_along(indexReturnVariables)) {		
		indexReturnVariable <- indexReturnVariables[i]

		Toolbox.logFine("indexReturnVariable: %s", indexReturnVariable)
		
		indepVariables <- Toolbox.histBetasIndepVars(indexReturnVariable, numLeads, numLags)
				
		allBetas[[indexReturnVariable]] <- Toolbox.histBetasCalculateIndiv(
				data, 
				groupVariable,
				returnVariable, 
				indepVariables,
				NULL,
				numPeriods,
				minimumObservations,
				calcSplitBetas = calcSplitBetas,
				useFast = useFast
		)
	}

	# Residual betas
	
	if(!is.null(baseIndexReturnVariable)) {
		Toolbox.logFine("Estimate residual betas")
		
		for(i in seq_along(indexReturnResidVariables)) {
			indexReturnResidVariable <- indexReturnResidVariables[i]
			
			indexReturnResidVariableEff <- paste0(indexReturnResidVariable, "_resid")

			indepVariables <- Toolbox.histBetasIndepVars(indexReturnResidVariableEff, numLeads, numLags)
			baseIndepVariables <- Toolbox.histBetasIndepVars(baseIndexReturnVariable, numLeads, numLags)

			allBetas[[indexReturnResidVariableEff]] <- Toolbox.histBetasCalculateIndiv(
					data, 
					groupVariable,
					returnVariable, 
					indepVariables,
					baseIndepVariables,
					numPeriods,
					minimumObservations,
					calcSplitBetas = calcSplitBetas,
					useFast = useFast
			)
		}		
	}

	if(length(allBetas) > 0) {	
		Toolbox.logFine("Merge estimated betas")	
		
		results <- dplyr::bind_rows(allBetas, .id = "index") %>%
				dplyr::distinct_at(c("index", groupVariable), .keep_all = TRUE)
	
		if(is.null(periodVariable)) {
			Toolbox.logFine("Map group to date")
			
			results <- results %>%
					dplyr::left_join(groupMapping, by = groupVariable) %>%
					dplyr::mutate(date_ = as.POSIXct(group_date_) + lubridate::period(1, periodicity)) %>%
					dplyr::rename_at("date_", ~ { dateVariable } ) %>%
					dplyr::select(-group_, -group_date_) 
		}
	} else {
		results <- NULL
	}
	
	if(!is.null(progressor)) {
		progressor()
	}
	
	return(results)
}

#' Construct lead/lag variables for all variables listed in indexVariables
Toolbox.histBetasCalculateLeadLag <- function(
		data, 
		dateVariable, 
		indexVariables,
		numLeads,
		numLags		
) {
	if(numLeads > 0) {
		for(lead in seq(numLeads)) {
			Toolbox.logFine("lead: %s", lead)
			
			data <- data %>%
					dplyr::arrange(dplyr::all_of(dateVariable)) %>%
					dplyr::mutate(
							dplyr::across(
									dplyr::all_of(indexVariables), 
									list(
											lead_ = ~ dplyr::lead(., lead)
									),
									.names = "{.col}_{.fn}"
							)
					) %>%
					dplyr::rename_at(
							paste0(indexVariables, "_lead_"), 
							~ paste0(., lead, "_")
					)		
		}
	}
	
	if(numLags) {
		for(lag in seq(numLags)) {
			Toolbox.logFine("lag: %s", lag)
			
			data <- data %>%
					dplyr::arrange(dplyr::all_of(dateVariable)) %>%
					dplyr::mutate(
							dplyr::across(
											dplyr::all_of(indexVariables), 
											list(
													lag_ = ~ dplyr::lag(., lag)
											),
											.names = "{.col}_{.fn}"
							)
					) %>%
					dplyr::rename_at(
							paste0(indexVariables, "_lag_"), 
							~ paste0(., lag, "_")
					)		
		}
	}	
	
	return(data)
}

#' @oaram returnsData data.frame containing security returns data 
#' @oaram indexReturnsData data.frame contain index returns data
#' @param numPeriods The number of periods (i.e. months) to include in estimation (see periodicity argument)
#' @param indexReturnLocalVariables The variable names of index returns to use in calculating local betas 
#' @param indexReturnResidLocalVariables The variable names of index returns to use when calculating local residual betas 
#' @param indexReturnUSDVariables The variable names of index returns to use in calculating USD betas 
#' @param indexReturnResidUSDVariables The variable names of index returns to use when calculating USD residual betas
#' @param baseIndexReturnVariable The variable name of index returns to include in regression when calculating residual betas
#' @param minimumObservations The minimum number of days to require when estimating betas
#' @param dateVariable The variable giving the return date (same for both security and index returns)
#' @param identifierVariables The variable names for security identifiers
#' @param returnLocalVariable The variable name for local returns
#' @param returnUSDVariable The variable name for usd returns
#' @param winsorizeReturns Winsorization levels for returns (in the time-series for each period separately).  Default is c(.01, .99). 
#' @param logReturns Whether to calculate betas on log returns (both securities and indices) 
#' @param periodicity The periodicity for estimated betas. One of year, month, week
#' @param dateOffset The number of days to be applied as an offset in order to calculate the period corresponding to the beta.  Typically 0 (for daily data), may be 7 for weekly data.
#' @param winsorizeBetas Winsorization levels for betas in th cross-section.  Default is c(.05, .95). 
#' @param calcSplitBetas Whether to estimate betas for positive and negative returns the factor (e.g. market)
#' @param numLeads The number of leads to apply to RHS return variable
#' @param numLags The number of lags to apply to RHS return variable
Toolbox.histBetasCalculateAll <- function(
		returnsData,
		indexReturnsData,
		numPeriods,
		indexReturnLocalVariables, 
		indexReturnResidLocalVariables, 
		indexReturnUSDVariables, 
		indexReturnResidUSDVariables, 
		baseIndexReturnVariable, 
		minimumObservations, 
		dateVariable = "date",
		identifierVariables = c("identifier", "exchange_country"),
		returnLocalVariable = "totalreturn_loc",
		returnUSDVariable = "totalreturn_usd",
		winsorizeReturns = c(.01, .99),
		logReturns = FALSE,
		periodicity = "month",
		dateOffset = 0,
		winsorizeBetas = c(.05, .95),
		calcSplitBetas = FALSE,
		numLeads = 1,
		numLags = 1,
		useFast = TRUE,
		parallel = FALSE
) {
	
	# TODO: Generalize split betas to "transformed" betas where a function is given to perform transformation.   
	#		- A post-transformation function would also need to be given to roll up pos/neg betas from components
	#		- This would facilitate ^ 2 and ^ 3 betas

	# TODO: Support exponentially weighted betas (via halfLife parameter)
	
	# Calculate lead and lag of index variables
	
	Toolbox.logInfo("Calculate lead and lag of index variables")
	
	# Convert to character in case it is a factor, otherwise a variable is appended as it's factor level, not value
	
	if(!is.null(indexReturnLocalVariables)) 
		indexReturnLocalVariables <- as.character(indexReturnLocalVariables) 
	if(!is.null(indexReturnUSDVariables)) 
		indexReturnUSDVariables <- as.character(indexReturnUSDVariables) 
	if(!is.null(indexReturnResidLocalVariables)) 
		indexReturnResidLocalVariables <- as.character(indexReturnResidLocalVariables) 
	if(!is.null(indexReturnResidUSDVariables)) 
		indexReturnResidUSDVariables <- as.character(indexReturnResidUSDVariables) 
	if(!is.null(baseIndexReturnVariable)) 
		baseIndexReturnVariable <- as.character(baseIndexReturnVariable)
	
	# Ensure that all index return variables have no missing values and are log-transformed as necessary
	
	indexReturnVars <- unique(c(
					indexReturnLocalVariables, 
					indexReturnUSDVariables, 
					indexReturnResidLocalVariables, 
					indexReturnResidUSDVariables, 
					baseIndexReturnVariable
			))
	
	Toolbox.logDebug("indexReturnVars: %s", indexReturnVars)

	indexReturnsData <- indexReturnsData %>%
			dplyr::ungroup() %>%
			dplyr::select(dplyr::all_of(c(indexReturnVars, dateVariable))) %>%
			dplyr::arrange(dplyr::all_of(dateVariable)) %>%
			dplyr::mutate(
					group_date_ = factor(lubridate::floor_date(eval(as.name(dateVariable)) + numLeads + dateOffset, periodicity)),
					group_ = as.integer(group_date_)
			) 	

	groupMapping <- indexReturnsData %>%
			dplyr::distinct(group_date_, group_)
	
	# Log transform, if required
	
	if(logReturns) {
		Toolbox.logDebug("Log-transforming index returns")
		
		indexReturnsData <- indexReturnsData %>%
				dplyr::mutate_at(indexReturnVars, ~ log(1 + .))		
	}
	
	# Create lagged values
	# Assign group_ integer based on the value of periodicity.

	Toolbox.logDebug("Calculate lead/lag index returns")

	# Only lead/lag the non-residual indexes

	indexReturnNonResidVars <- unique(c(
					indexReturnLocalVariables, 
					indexReturnUSDVariables, 
					baseIndexReturnVariable
			))
	
	# Calculate the leads/lags of the index returns

	indexReturnsData <- Toolbox.histBetasCalculateLeadLag(
			indexReturnsData, 
			dateVariable, 
			indexReturnNonResidVars,
			numLeads = numLeads,
			numLags = numLags
	)

	# Purge residual index returns of the base index returns (include lead and lag values) and
	# append to indexReturnsData
	
	indexReturnResidVariables <- unique(c(indexReturnResidLocalVariables, indexReturnResidUSDVariables))
	
	if(!is.null(baseIndexReturnVariable) && length(indexReturnResidVariables) > 0) {
		Toolbox.logInfo("Purge residual index return variables of base index returns")	
		
		indepVariables <- Toolbox.histBetasIndepVars(baseIndexReturnVariable, numLeads, numLags)

		Toolbox.logDebug("indepVariables: %s", indepVariables)
		
		for(i in seq_along(indexReturnResidVariables)) {
			indexReturnResidVariable <- indexReturnResidVariables[i]

			Toolbox.logDebug("indexReturnResidVariable: %s", indexReturnResidVariable)

			indexReturnsDataSub <- indexReturnsData %>%
					dplyr::filter_at(dplyr::all_of(c(indexReturnResidVariable, indepVariables)), dplyr::all_vars(is.finite(.)))
			
			x <- cbind(1, as.matrix(indexReturnsDataSub[indepVariables]))
			y <- as.matrix(indexReturnsDataSub[indexReturnResidVariable])			
			
			colnames(x) <- c('(Intercept)', indepVariables)
			
			# TODO: Validate the approach of using rolling regression to get residual index returns.
		
			Toolbox.logWarn("Validate the approach of using rolling regression to get residual index returns.")
			
			fit <- Toolbox.calcRollingRegressionFast(
					x, 
					y, 
					numPeriods, 
					fullStats = TRUE,
					groups = indexReturnsDataSub$group_, 
					rollingWindow = TRUE,
					minimumObservations = minimumObservations,
					partial = TRUE,
					useFast = useFast
			)
			
			# Calculate residual returns using coeffs estimated in each respective group
			# This is joined into indexReturnsDataSub in case data are dropped due to missing values
		
			indexReturnsDataSub$residual_ <- (y - rowSums(x * fit$coefs))[,1]
									
			indexReturnsData <- indexReturnsData %>%
					dplyr::left_join(indexReturnsDataSub %>%
									dplyr::select(dplyr::all_of(c(dateVariable, "residual_"))) %>%
									dplyr::rename_at("residual_", ~ { paste0(indexReturnResidVariable, "_resid") } ), 
							by = dateVariable)
		}
		
		indexReturnsData <- Toolbox.histBetasCalculateLeadLag(
				indexReturnsData, 
				dateVariable, 
				paste0(indexReturnResidVariables, "_resid"),
				numLeads = numLeads,
				numLags = numLags
		)
	}
	
	Toolbox.logDebug("Construct model data")
	
	modelData <- returnsData %>%
			dplyr::select(dplyr::all_of(c(identifierVariables, dateVariable, returnLocalVariable, returnUSDVariable))) %>%
			dplyr::filter_at(c(returnLocalVariable, returnUSDVariable), dplyr::all_vars(is.finite(.))) %>%
			dplyr::inner_join(indexReturnsData, by = dateVariable)

	if(logReturns) {
		Toolbox.logDebug("Log-transforming security returns")
		
		modelData <- modelData %>%
				dplyr::mutate_at(c(returnLocalVariable, returnUSDVariable), ~ log(1 + .))
	}

	Toolbox.logDebug("Nest model data")
	
	modelDataNest <- modelData %>%
			dplyr::group_by_at(identifierVariables) %>%
			tidyr::nest()	
 
	# NOTE: using furrr doesn't appear to speed things up at all.  IN fact it seems to hang while 
	#			initializing the first worker, not making any progress.
	
	if(parallel) {
		Toolbox.logDebug("Using furrr::future_map()")
		
		mapFunc <- furrr::future_map
	} else {
		Toolbox.logDebug("Using purrr::map()")
		
		mapFunc <- purrr::map
	}
	
	Toolbox.logDebug("nrow(modelDataNest):: %s", nrow(modelDataNest))

	# TODO: These (USD and local calculations) should be wrapped into a single function which is executed in mapFunc() so that 
	# 		in the parallel situation we only have one call to mapFunc() which should be faster given the need to
	#		initialize workers as a fixed cost.
	
	if(!is.null(returnUSDVariable) && !is.null(indexReturnUSDVariables)) {
		p <- progressr::progressor(steps = nrow(modelDataNest), message = "Calculating USD Betas")
		
		Toolbox.logDebug("Calculate USD betas")
		
		allBetasUSD <- modelDataNest %>%
				dplyr::mutate(
						beta = mapFunc(
								data, 
								~ Toolbox.histBetasCalculate(
										.,
										dateVariable,
										returnUSDVariable, 
										indexReturnUSDVariables, 
										indexReturnResidUSDVariables, 
										baseIndexReturnVariable,
										numPeriods,
										minimumObservations,
										periodVariable = "group_",
										periodicity = periodicity,
										dateOffset = dateOffset,
										winsorizeReturns = winsorizeReturns,
										calcSplitBetas = calcSplitBetas,
										numLeads = numLeads,
										numLags = numLags,
										useFast = useFast,
										progressor = p
							)
						)) %>%
				dplyr::select(-data) %>%
				tidyr::unnest(beta) %>%
				dplyr::mutate(currency = "usd")
	} else {
		allBetasUSD <- NULL
	}

	if(!is.null(returnLocalVariable) && !is.null(indexReturnLocalVariables)) {
		p <- progressr::progressor(steps = nrow(modelDataNest), message = "Calculating Local Betas")
		
		Toolbox.logDebug("Calculate Local betas")
						
		allBetasLocal <- modelDataNest %>%
				dplyr::mutate(
						beta = mapFunc(data, 
								~ Toolbox.histBetasCalculate(
										.,
										dateVariable,
										returnLocalVariable, 
										indexReturnLocalVariables, 
										indexReturnResidLocalVariables, 
										baseIndexReturnVariable,
										numPeriods,
										minimumObservations,
										periodVariable = "group_",
										periodicity = periodicity,
										dateOffset = dateOffset,
										winsorizeReturns = winsorizeReturns,
										calcSplitBetas = calcSplitBetas,
										numLeads = numLeads,
										numLags = numLags,
										useFast = useFast,
										progressor = p)
						)) %>%
				dplyr::select(-data) %>%
				tidyr::unnest(beta) %>%
				dplyr::mutate(currency = "loc")
	} else {
		allBetasLocal <- NULL
	}
	
	Toolbox.logInfo("Combine betas from all currencies")	
	
	results <- dplyr::bind_rows(allBetasUSD, allBetasLocal)
	
	if(!is.null(winsorizeBetas)) {		
		Toolbox.logInfo("Winsorize betas")	

		if(calcSplitBetas) {
			results <- results %>%
					dplyr::group_by(index, currency, group_) %>%
					dplyr::mutate(
							beta_neg_win = Toolbox.winsorize(beta_neg, winsorizeBetas[1], winsorizeBetas[2]),
							beta_pos_win = Toolbox.winsorize(beta_pos, winsorizeBetas[1], winsorizeBetas[2])
					)
			
			if(!is.null(baseIndexReturnVariable) && length(c(indexReturnResidLocalVariables, indexReturnResidUSDVariables)) > 0) {
				results <- results %>%
						dplyr::group_by(index, currency, group_) %>%
						dplyr::mutate(
								beta_base_neg_win = Toolbox.winsorize(beta_base_neg, winsorizeBetas[1], winsorizeBetas[2]),
								beta_base_pos_win = Toolbox.winsorize(beta_base_pos, winsorizeBetas[1], winsorizeBetas[2])
						)			
			}
		} else {
			results <- results %>%
					dplyr::group_by(index, currency, group_) %>%
					dplyr::mutate(
							beta_win = Toolbox.winsorize(beta, winsorizeBetas[1], winsorizeBetas[2])				
					)
			
			if(!is.null(baseIndexReturnVariable) && length(c(indexReturnResidLocalVariables, indexReturnResidUSDVariables)) > 0) {
				results <- results %>%
						dplyr::group_by(index, currency, group_) %>%
						dplyr::mutate(
								beta_base_win = Toolbox.winsorize(beta_base, winsorizeBetas[1], winsorizeBetas[2])				
						)			
			}			
		}
	}

	Toolbox.logFine("Map group to date")
	
	results <- results %>%
			dplyr::ungroup() %>%
			dplyr::left_join(groupMapping, by = "group_") %>%
			dplyr::mutate(date_ = as.POSIXct(group_date_) + lubridate::period(1, periodicity)) %>%
			dplyr::rename_at("date_", ~ { dateVariable } ) %>%
			dplyr::select(-group_, -group_date_)
	
	return(results)
}


