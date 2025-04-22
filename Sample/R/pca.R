# Functions to support the calculation and analysis of PCA-based models
# 
###############################################################################

#' Construct list containing model data for PCA estimation.   This also includes the
#' calculation of the covariance matrix for the transposed security returns.  This
#' is split off from Toolbox.pcaCalculatePCAModel() to improve performance across multiple 
#' invocations (eliminating the need to re-calculate the cov matrix)
Toolbox.pcaConstructModelData <- function(
		dataRet, 
		identifierVariables, 
		dateVariable, 
		returnVariableUSD,
		returnVariableLocal,
		marketFactor = TRUE,
		fxFactors = TRUE,
		centerReturns = FALSE, 
		scaleReturns = FALSE
) {
	Toolbox.logDebug("Get valid return subset")
	
	# Eliminate observations with missing returns
	
	returnVariables <- c(returnVariableUSD, returnVariableLocal)
	
	dataRetSub <- dataRet %>%
			dplyr::filter_at(dplyr::all_of(returnVariables), dplyr::all_vars(is.finite(.))) %>%
			dplyr::distinct(!!!dplyr::syms(c(identifierVariables, dateVariable)), .keep_all = TRUE)

	if(fxFactors && !is.null(returnVariableLocal)) {
		Toolbox.logDebug("Calculate currency factor returns")

		fxReturns <- dataRetSub %>%
				dplyr::group_by_at(c(dateVariable, "exchange_country")) %>%
				dplyr::summarize(
						return = Toolbox.fuzz(
								mean((1 + eval(as.name(returnVariableUSD))) / ( 1 + eval(as.name(returnVariableLocal))) - 1)
						)
				)
				
		returnVariableEff <- returnVariableLocal
	} else {
		returnVariableEff <- returnVariableUSD		
		fxReturns <- NULL
	}
	
	if(marketFactor) {
		Toolbox.logDebug("Calculate market factor")
		
		# Calculate returns to (equal weight) market and subtract from security returns in each period
		
		marketReturns <- dataRetSub %>%
				dplyr::group_by_at(dateVariable) %>%
				dplyr::summarize(market_ = mean(eval(as.name(returnVariableEff))))
		
		dataRetSub <- dataRetSub %>%
				dplyr::left_join(marketReturns, by = dateVariable) %>%
				dplyr::mutate(return__ = eval(as.name(returnVariableEff)) - market_) %>%
				dplyr::select(- market_)
		
		returnVariableEff <- "return__"		
	} else {
		marketReturns <- NULL
	}	
	
	Toolbox.logDebug("Transpose returns")
	
	dataRetSubT <- dataRetSub %>%
			dplyr::select_at(c(identifierVariables, dateVariable, returnVariableEff)) %>%
			dplyr::arrange_at(dateVariable) %>%
			tidyr::spread(dateVariable, returnVariableEff)
	
	Toolbox.logDebug("Construct identifier mapping")
	
	identifierMap <- dataRetSub %>%
			dplyr::ungroup() %>%
			dplyr::select(dplyr::all_of(identifierVariables)) %>%
			tidyr::unite(
					identifier_full_, 
					dplyr::all_of(identifierVariables), sep = "|", remove = FALSE
			) %>%
			dplyr::distinct(identifier_full_, .keep_all = TRUE)
	
	Toolbox.logFine("Convert to matrix")
	
	dataRetSubTM <- scale(
			as.matrix(
					dataRetSubT %>% 
							dplyr::select(-dplyr::all_of(identifierVariables))
			), 
			center = centerReturns, 
			scale = scaleReturns
	)
	
	if(is.unsorted(colnames(dataRetSubTM))) {
		Toolbox.logError("Dates in dataRetSubTM not in ascending order")
		stop("Dates in dataRetSubTM not in ascending order")
	}
	
	returnCenter <- attr(dataRetSubTM, "scaled:center")
	returnScale <- attr(dataRetSubTM, "scaled:scale")
	
	rownames(dataRetSubTM) <- apply( dataRetSubT[ ,identifierVariables], 1, paste, collapse = "|")
	
	Toolbox.logDebug("Calculate covariance")
	
	# dof: Count of non-missing returns
	# cov: The covariance matrix (actually a design or correlation matrix depending on values of center and scale)
	
	dataRetSubTMFinite <- is.finite(dataRetSubTM)
	dof <- t(dataRetSubTMFinite) %*% dataRetSubTMFinite
	
	# Set missing returns to 0
	
	dataRetSubTM[!is.finite(dataRetSubTM)] <- 0
	
	# NOTE: Force use of BLAS because we know that there aren't any non-finite values and this will be fasted
	
	oldMatProdOpt <- getOption("matprod")
	options(matprod = "blas")
	cov <- crossprod(dataRetSubTM) / (dof - 1)
	options(matprod = oldMatProdOpt)
	
	# Set missing values to 0 (shouldn't be possible)
	
	cov[!is.finite(cov)] <- 0
	
	if(any(colnames(dataRetSubTM) != colnames(cov)) || any(colnames(dataRetSubTM) != rownames(cov))) {
		Toolbox.logError("Column names dataRetSubTM do not match col/row names of cov")
		stop("Column names dataRetSubTM do not match col/row names of cov")
	}
	
	allDates <- sort(unique(dataRetSub[[dateVariable]]))
	
	results <- list(
		marketReturns = marketReturns,
		fxReturns = fxReturns,
		
		dataRetSubTM = dataRetSubTM,
		dataRetSubTMFinite = dataRetSubTMFinite,
		
		returnCenter = returnCenter,
		returnScale = returnScale,

		cov = cov,
		
		identifierMap = identifierMap,
		allDates = allDates
	)
	
	return(results)
}

#' Calculate PCA risk model
#' @param window The maximum number of days in estimation window 
#' @param minWindow The minimum number of days in a window.  1 means that we start from the beginning (i.e. partial)
#' @param annualizationFactor The scale factor applied to annualized statistics.  (daily scale = 252, monthly scale = 12, quarterly scale = 4)
#' @param periodicity takes ("month", "week", "day").  Determines how often the PCA analysis is computed.  Dramatically saves on time to only compute "monthly" vs. "daily"
#' @param requiredWeight The minium required total weight (as a fraction of of a security to have it's loadings calculated 
Toolbox.pcaCalculatePCAModel <- function(
		modelData,
		startDate = NULL,
		endDate = NULL,
		numFactors = NULL, 
		halfLife = NULL, 
		window = NULL, 
		minWindow = NULL,
		marketFactor = TRUE,
		fxFactors = TRUE,
		centerLoadings = TRUE, 
		scaleLoadings = TRUE,
		annualizationFactor = 1, 
		periodicity = "month",
		requiredWeight = .75,
		requiredWeightWindow = NULL,
		residualReturns = FALSE,
		deflippify = TRUE
) {
	
	# TODO: Replicate (or improve on) SAS implementation's checks for missing values
	# TODO: Allow specification of additional "fixed" factors that will be used to purge returns and then whose factor returns are also returned (along with covariances)
	# TODO: Extract the code through calculation of covariance matrix to a separate function so that it can be re-used across multiple invocations of this function
	
	Toolbox.logDebug("annualizationFactor: %s", annualizationFactor)

	# Extract model data variables

	marketReturns <- modelData$marketReturns
	
	if(fxFactors) {
		fxReturnsT <- modelData$fxReturns %>%
				Toolbox.pivotWide(
						"datadate", 
						"exchange_country", 
						"return", 
						fillValue = 0,
						prefix = "curr_"
				)
	}
	
	dataRetSubTM <- modelData$dataRetSubTM
	dataRetSubTMFinite <- modelData$dataRetSubTMFinite
	cov <- modelData$cov
	identifierMap <- modelData$identifierMap
	allDates <- modelData$allDates
		
	# TODO: Make sure market and fx returns align with modelData dates 
	
	# If no window given then assume maximum length
	
	if(is.null(window)) {		
		window <- length(allDates)
	}

	Toolbox.logDebug("window: %s", window)
	
	# If minWindow is null then assume minimum window size is window
	
	if(is.null(minWindow)) {
		minWindow <- window
	}

	Toolbox.logDebug("minWindow: %s", minWindow)
	
	if(is.null(requiredWeightWindow)) {
		requiredWeightWindow <- minWindow
	}

	requiredWeightWindow <- min(requiredWeightWindow, window)
	
	Toolbox.logDebug("requiredWeightWindow: %s", requiredWeightWindow)
	
	# Get all dates for which model should be estimated

	# Note, early dates are excluded based on whether a partial estimation window is allowed

	allComputationDates <- as.POSIXct(sort(unique(cut(tail(allDates, -minWindow), periodicity)))[-1])
	
	if(!is.null(startDate)) {
		allComputationDates <- allComputationDates[allComputationDates >= startDate]
	}
	
	if(!is.null(endDate)) {
		allComputationDates <- allComputationDates[allComputationDates <= endDate]
	}
	
	# If no numFactors then assume window since that would be the maximum

	if(is.null(numFactors)) {
		numFactors <- window
	}
	
	Toolbox.logDebug("numFactors: %s", numFactors)
	
	firstModelIndex <- 1
	
	Toolbox.logDebug("firstModelDateIndex: %s", firstModelIndex)
	
	lastModelIndex <- length(allComputationDates)
	
	Toolbox.logDebug("lastModelDateIndex: %s", lastModelIndex)
	
	allModels <- Toolbox.createList(lastModelIndex - firstModelIndex + 1)
	
	factLoadingsStdLast <- NULL
	
	for(modelIndex in seq(firstModelIndex, lastModelIndex)) {
		Toolbox.logFine("modelIndex: %s", modelIndex)		
		
		currDate <- allComputationDates[modelIndex]

		Toolbox.logDebug("currDate: %s", currDate)
		
		# Identify the last date based on the currdate provided by the modelIndex

		lastDateIndex <- max(which(allDates < currDate))		
		
		Toolbox.logFine("lastDateIndex: %s", lastDateIndex)
		
		fromDateIndex <- max(1, lastDateIndex - window + 1) 
		
		Toolbox.logFine("fromDateIndex: %s", fromDateIndex)		
		
		subDates <- allDates[fromDateIndex:lastDateIndex]
		
		if(length(subDates) > 1) {
			Toolbox.logFine("min(subDates): %s", min(subDates))		
			Toolbox.logFine("max(subDates): %s", max(subDates))		
		
			if(marketFactor) {
				marketReturnsSub <- marketReturns$market_[fromDateIndex:lastDateIndex]
			}
			
			if(fxFactors) {
				fxReturnsTSub <- fxReturnsT[fromDateIndex:lastDateIndex,]
			}
			
			covSub <- cov[fromDateIndex:lastDateIndex, fromDateIndex:lastDateIndex]
			
			if(Toolbox.safeIsFinite(halfLife)) {
				Toolbox.logFine("Apply half-life to covSub values")
				
				dateDiffs <- (as.numeric(currDate) - as.numeric(subDates)) / (24 * 60 * 60)
				weights <- exp( - dateDiffs * log(2) / halfLife)
				weights <- length(weights) * weights / sum(weights)
				
				covSub <- covSub * (weights %*% t(weights))
			} else {
				weights <- rep(1, length(subDates))
			}

			# Force to sum to N
					
			numFactorsEff <- min(numFactors, lastDateIndex - fromDateIndex + 1)
			
			Toolbox.logFine("numFactorsEff: %s", numFactorsEff)

			Toolbox.logFine("Calculate eigenvector decomposition")
			
			eigenResults <- eigen(covSub)
			
			factReturns <- eigenResults$vectors[,seq(1:numFactorsEff)]
			
			# factorNames is a "factor" to ensure that they are plotted in the correct order
			
			factorNames <- paste0("PC", seq(1, numFactorsEff))
			factorNames <- factor(factorNames, levels = factorNames)
			
			colnames(factReturns) <- factorNames
						
			factStdDev <- sqrt(eigenResults$values[seq(1:numFactorsEff)])
			
			factFraction <- eigenResults$values[seq(1:numFactorsEff)] / sum(eigenResults$values)
			
			# Calculate security loadings on factors
			
			Toolbox.logFine("Calculate security loadings on factors")
			
			dataRetSubTMFiniteSub <- dataRetSubTMFinite[, max(1, lastDateIndex - requiredWeightWindow + 1):lastDateIndex]
			weightFiniteSub <- tail(weights, requiredWeightWindow)
						
			# The fraction of the weight applied to returns for each security that corresponds to finite (non-missing)
			# returns.  If finite returns have > .75 of weight then these securities have valid loadings.  In exponentially-weighted
			# models this has the nice effect of valuing more recent returns more than older returns.

			finiteRetWeights <- apply(weightFiniteSub * dataRetSubTMFiniteSub, 1, sum) / sum(weightFiniteSub)

			dataRetSubTMSub <- dataRetSubTM[finiteRetWeights > requiredWeight, fromDateIndex:lastDateIndex]
			
			# NOTE: Missing (non-finite) returns are set to 0 above which has the effect of dropping them
			# 		from the regression
			
			# Weight the security returns

			dataRetSubTMSubWtd <- t(t(dataRetSubTMSub) * weights) 
						
			# NOTE: This is done directly (instead of using lm()) for performance reasons
	
			Toolbox.logFine("Calculating ginv()")
	
			factLoadings <- t(ginv(t(factReturns) %*% factReturns) %*% t(factReturns) %*% t(dataRetSubTMSubWtd))
			
			colnames(factLoadings) <- factorNames
			
			# Standardize loadings
	
			Toolbox.logFine("Standardize loadings")		
	
			factLoadingsStd <- scale(factLoadings, center = centerLoadings, scale = scaleLoadings)
			
			factLoadingsCenter <- attr(factLoadingsStd, "scaled:center")
			factLoadingsScale <- attr(factLoadingsStd, "scaled:scale")
						
			if(deflippify && !is.null(factLoadingsStdLast)) {
				Toolbox.logFine("Deflippify factor loadings")		

				# Detect whether factor loadings have "flipped" relative to the prior period's and, if so
				# then flip then back.  This is an issue because the sign of factors are not
				# uniquely defined
		
				factLoadingsStdT <- as.data.frame(factLoadingsStd) %>%
						tibble::rownames_to_column("identifier_full_") %>%
						dplyr::group_by(identifier_full_) %>%
						Toolbox.pivotNarrow("factor_id", "exposure")
				
				factLoadingsStdLastT <- as.data.frame(factLoadingsStdLast) %>%
						tibble::rownames_to_column("identifier_full_") %>%
						dplyr::group_by(identifier_full_) %>%
						Toolbox.pivotNarrow("factor_id", "exposure_last")

				# TODO: Have a threshold for abs(corr) so that it only flips back if the corr is high enough in maginitude
				#			that it is a _flip_ rather than a change in factor identity
	
				factLoadingsCorr <- factLoadingsStdT %>%
						dplyr::inner_join(factLoadingsStdLastT, by = c("identifier_full_", "factor_id")) %>%
						dplyr::group_by(factor_id) %>%
						dplyr::summarize(
								corr = cor(exposure, exposure_last, use = "pairwise.complete"),
								corr_sign = ifelse(corr == 0, 1, sign(corr))
						)
				
				# Missing corr_sign can happen if the # of factors change from period to period.  This
				# can happen when # factors > minWindow in early periods
	
				factLoadingsStdAdj <- factLoadingsStdT %>%
						dplyr::left_join(factLoadingsCorr, by = "factor_id") %>%
						dplyr::mutate(corr_sign = dplyr::coalesce(corr_sign, 1)) %>%
						dplyr::mutate(exposure_adj = exposure * corr_sign) %>%
						Toolbox.pivotWide("identifier_full_", "factor_id", "exposure_adj", fillValue = NA) %>%
						tibble::column_to_rownames("identifier_full_") %>%
						as.matrix
			} else {
				factLoadingsStdAdj <- factLoadingsStd
			}
			
			factLoadingsStdLast <- factLoadingsStdAdj
			
			# Calculate factor returns based on loadings
			# These may vary from factReturns due to standardization, security coverage, etc. but mostly
			# standardization of the loadings

			Toolbox.logFine("Calculate effective factor returns")		

			factReturnsEff <- t(ginv(t(factLoadingsStdAdj) %*% factLoadingsStdAdj) %*% t(factLoadingsStdAdj) %*% dataRetSubTMSubWtd)
									
			colnames(factReturnsEff) <- factorNames
						
			# Calculate residual returns and specific risk

			Toolbox.logFine("Calculate residual returns")		

			residReturns <- t(dataRetSubTMSub) - factReturnsEff %*% t(factLoadingsStdAdj)

			Toolbox.logFine("Calculate specific risk")		
			
			specRisk <- apply(
					residReturns, 
					2, 
					function(x) { 
						c(
								sd(x, na.rm = TRUE) * sqrt(annualizationFactor),
								sqrt(Hmisc::wtd.var(x, weights = weights, normwt = TRUE, na.rm = TRUE) * annualizationFactor)
						)
					}
			)
			
			rownames(specRisk) <- c("srisk", "srisk_weighted")

			# Add market returns to factor returns (if necessary)

			if(marketFactor) {
				Toolbox.logFine("Add market factor")		
				
				factorNamesEff <- c("market", as.character(factorNames))
				factorNamesEff <- factor(factorNamesEff, levels = factorNamesEff)
				
				factReturnsEff <- cbind(
						marketReturnsSub * weights,
						factReturnsEff
				)
				
				colnames(factReturnsEff) <- factorNamesEff
				
				factLoadingsStdAdj <- cbind(
						1, 
						factLoadingsStdAdj
				)
				
				colnames(factLoadingsStdAdj) <- factorNamesEff
			} else {
				factorNamesEff <- factorNames
			}

			if(fxFactors) {
				Toolbox.logFine("Add FX factors")		
				
				fxFactorNames <- colnames(fxReturnsTSub)[-1]
				
				factorNamesEff <- c(as.character(factorNamesEff), fxFactorNames)
				factorNamesEff <- factor(factorNamesEff, levels = factorNamesEff)
				
				factReturnsEff <- cbind(
						factReturnsEff,
						fxReturnsTSub[-1] * weights						
				)
				
				colnames(factReturnsEff) <- factorNamesEff
				
				fxFactorLoadings <- data.frame(identifier_full_ = rownames(factLoadingsStdAdj)) %>%
						dplyr::left_join(identifierMap, by = "identifier_full_") %>%
						Toolbox.constructDummyVariables(
								variable = "exchange_country", 
								prefix = "curr_", 
								required = gsub("^curr_", "", fxFactorNames)
						) %>%
						dplyr::select_at(fxFactorNames)
				
				factLoadingsStdAdj <- cbind(
						factLoadingsStdAdj,
						fxFactorLoadings
				)
				
				colnames(factLoadingsStdAdj) <- factorNamesEff
			}
			
			
			# Calculate annualized covariance matrix
	
			Toolbox.logFine("Calculate factor covariance matrix")		
	
			# TODO: Confirm that the "/ sqrt(weights)" works correctly for the PC# factors
			
			factCov <- cov(factReturnsEff / sqrt(weights), use = "pairwise.complete") * annualizationFactor

			Toolbox.logFine("Construct factFractionDF")		
			
			factFractionDF <- data.frame(
					factor = factorNames, 
					fraction = factFraction, 
					stringsAsFactors = FALSE
			)

			Toolbox.logFine("Construct factStdDevDF")		
			
			factStdDevDF <- data.frame(
					factor = factorNames, 
					std_dev = factStdDev, 
					stringsAsFactors = FALSE
			)

			Toolbox.logFine("Construct factLoadingsStdDF")		
			
			factLoadingsStdDF <- data.frame(
							factLoadingsStdAdj, 
							stringsAsFactors = FALSE
					) %>%
					tibble::rownames_to_column("identifier_full_") %>%
					dplyr::left_join(identifierMap, by = "identifier_full_") %>%
					dplyr::select(-identifier_full_)

			Toolbox.logFine("Construct factReturnsDF")		
			
			factReturnsDF <- data.frame(
					factReturns, 
					date = subDates, 
					stringsAsFactors = FALSE
			)

			Toolbox.logFine("Construct factReturnsEffDF")		
			
			factReturnsEffDF <- data.frame(
					factReturnsEff, 
					date = subDates, 
					stringsAsFactors = FALSE
			)

			# Residual returns which may be of use if PCA is being used to "purge" return of PCA factor exposures

			Toolbox.logFine("Construct residRetsDF")		

			if(residualReturns) {
				residRetsDF <- Toolbox.pivotNarrow(
									data.frame(
										date = subDates,
										residReturns, 
										stringsAsFactors = FALSE,
										check.names = FALSE
									) %>% dplyr::group_by(date),
									"identifier_full_", 
									"return_res"
								) %>%
								dplyr::left_join(identifierMap, by = "identifier_full_") %>%
								dplyr::select(-identifier_full_)
			} else {
				residRetsDF <- NULL
			}

			Toolbox.logFine("Construct specRiskDF")		
					
			specRiskDF <- data.frame(t(specRisk)) %>%
						tibble::rownames_to_column("identifier_full_") %>%
						dplyr::left_join(identifierMap, by = "identifier_full_") %>%
						dplyr::select(-identifier_full_)					

			Toolbox.logFine("Construct model results")		
					
			allModels[[modelIndex]] <- list(
					fraction = factFractionDF,
					stdDev = factStdDevDF,

					weights = weights,
					
					loadings = factLoadingsStdDF,
					loadingsCenter = factLoadingsCenter, 
					loadingsScale = factLoadingsScale, 
					
					factReturnsRaw = factReturnsDF, 
					factReturns = factReturnsEffDF,
					
					covRet = covSub,
					
					cov = factCov,
					
					returnRes = residRetsDF,
					specRisk = specRiskDF
			)
		} else {
			allModels[[modelIndex]] <- NULL					
		}
	}

	Toolbox.logDebug("Estimation finished")		
	
	results <- list(
			models = allModels,
			allDates = allComputationDates[seq(firstModelIndex, lastModelIndex)]
	) 
	
	return(results)
}

#' Flatten the output from Toolbox.pcaCalculatePCAModel() to facilitate analysis of
#' loadings and factor returns
Toolbox.pcaCombinePCAModel <- function(
		pcaModel, 
		wide = TRUE
) {
	datesDF <- data.frame(period_dt = pcaModel$allDates, period_index = seq_along(pcaModel$allDates))
	
	Toolbox.logDebug("Combine loadings")
	
	loadingsDF <- dplyr::bind_rows(
					lapply(
							pcaModel$models, 
							function(x) {
								results <- x$loadings
								
								if(!wide) {
									results <- results %>%
											dplyr::group_by(identifier, exchange_country) %>%
											Toolbox.pivotNarrow(
													"factor",
													"value"
											)
								}
								
								return(results)
							}
					),
					.id = "period_index"
			) %>%
			dplyr::mutate(period_index = as.integer(period_index)) %>%
			dplyr::left_join(datesDF, by = "period_index") %>%
			dplyr::select(-period_index)

	Toolbox.logDebug("Combine specRisk")
	
	specRiskDF <- dplyr::bind_rows(
					lapply(
							pcaModel$models, 
							function(x) {
								return(x$specRisk)
							}
					),
					.id = "period_index"
			) %>%
			dplyr::mutate(period_index = as.integer(period_index)) %>%
			dplyr::left_join(datesDF, by = "period_index") %>%
			dplyr::select(-period_index)
	
	Toolbox.logDebug("Combine factReturns")
	
	factReturnsDF <- dplyr::bind_rows(
					lapply(
							pcaModel$models, 
							function(x) {
								results <- x$factReturns
								
								if(!wide) {
									results <- results %>%
											dplyr::group_by(date) %>%
											Toolbox.pivotNarrow(
													"factor",
													"return"
											)
								}
								
								return(results)								
							}
					),
					.id = "period_index"
			) %>%
			dplyr::mutate(period_index = as.integer(period_index)) %>%
			dplyr::left_join(datesDF, by = "period_index") %>%
			dplyr::select(-period_index)

	Toolbox.logDebug("Combine cov")
	
	covDF <- dplyr::bind_rows(
					lapply(
							pcaModel$models, 
							function(x) {
								results <- data.frame(x$cov) %>%
										dplyr::add_rownames("name_row")

								if(!wide) {
									results <- results %>%
											dplyr::group_by(name_row) %>%
											Toolbox.pivotNarrow(
													"name_col",
													"cov"
											)
								}
								
								return(results)								
							}
					),
					.id = "period_index"
			) %>%
			dplyr::mutate(period_index = as.integer(period_index)) %>%
			dplyr::left_join(datesDF, by = "period_index") %>%
			dplyr::select(-period_index)
	
	Toolbox.logDebug("Combine faction")
	
	fractionDF <- dplyr::bind_rows(
					lapply(
							pcaModel$models, 
							function(x) {
								return(x$fraction)
							}
					),
					.id = "period_index"
			) %>%
			dplyr::mutate(period_index = as.integer(period_index)) %>%
			dplyr::left_join(datesDF, by = "period_index") %>%
			dplyr::select(-period_index)
	
	Toolbox.logDebug("Combine stdDev")
	
	stdDevDF <- dplyr::bind_rows(
					lapply(
							pcaModel$models, 
							function(x) {
								return(x$stdDev)
							}
					),
					.id = "period_index"
			) %>%
			dplyr::mutate(period_index = as.integer(period_index)) %>%
			dplyr::left_join(datesDF, by = "period_index") %>%
			dplyr::select(-period_index)
	
	Toolbox.logDebug("Combine loadingStats")
	
	loadingsStatsDF <- dplyr::bind_rows(
					lapply(
							pcaModel$models, 
							function(x) {
								results <- data.frame(
										center = x$loadingsCenter,
										scale = x$loadingsScale,
										factor = names(x$loadingsCenter)
								)
								
								return(results)
							}
					),
					.id = "period_index"
			) %>%
			dplyr::mutate(period_index = as.integer(period_index)) %>%
			dplyr::left_join(datesDF, by = "period_index") %>%
			dplyr::select(-period_index)
	
	Toolbox.logDebug("Construct metadata")
	
	riskFactorMetadata <- data.frame(
					name = c("srisk_weighted", head(colnames(loadingsDF), -3))
			) %>%
			dplyr::mutate(
					type = dplyr::case_when(
							name == "srisk_weighted" ~ "specific",
							name == "market" ~ "market",
							grepl("^PC.*", name) ~ "common",
							grepl("^curr_.*", name) ~ "currency"
					)
			)
	
	results <- list(
			loadings = loadingsDF,
			specRisk = specRiskDF,
			factReturns = factReturnsDF,
			cov = covDF,
			fraction = fractionDF,
			loadingsStats = loadingsStatsDF,
			riskFactorMetadata = riskFactorMetadata
	)
	
	return(results)
}
