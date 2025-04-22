# TODO: Add comment
# 
# Author: kcole
###############################################################################

suppressWarnings(suppressMessages({								
	library(car, warn.conflicts = FALSE, quietly = TRUE)
    library(segmented, warn.conflicts = FALSE, quietly = TRUE)
}))

#' Compute robust standard errors for the results from lm()
#' Stolen from: https://gist.github.com/tony91782/954316
#' @param model The model from lm()
#' @param type One of hc0:hc4
Toolbox.calculateRobustStdErr <- function(model, type=c("hc3", "hc0", "hc1", "hc2", "hc4"), ...){	
	type <- match.arg(type)
	V <- hccm(model, type=type)
	sumry <- summary(model)
	table <- coef(sumry)
	table[,2] <- sqrt(diag(V))
	table[,3] <- table[,1]/table[,2]
	table[,4] <- 2*pt(abs(table[,3]), df.residual(model), lower.tail=FALSE)
	
	sumry$coefficients <- table
	p <- nrow(table)
	hyp <- cbind(0, diag(p - 1))
	sumry$fstatistic[1] <- linearHypothesis(model, hyp,white.adjust=type)[2,"F"]
	sumry$type <- type
	
	return(sumry)
}

#' Compute the R-squared for the given model fit object.  This is necessary in the 
#' case of WLS since the R-squared computation that R performs is on the weighted 
#' data and it is generally much larger than the R-squared computed on the unweighted, 
#' which is usually how it is being used.  That is, as a measure of the amount of
#' variation explained.
#' See http://gseacademic.harvard.edu/~willetjo/pdf%20files/Willett_Singer_AS88.pdf
#' See https://stats.stackexchange.com/questions/83826/is-a-weighted-r2-in-robust-linear-model-meaningful-for-goodness-of-fit-analys
Toolbox.calculateRSquared <- function(fit) {  
	SSe <- sum((fit$residuals)^2)  
	observed <- fit$residuals+fit$fitted.values
	SSt <- sum((observed-mean(observed))^2)  
	value <- 1 - SSe / SSt
	return(value)
}  

#' Compute the VIF for all regressors in the given lm fit
Toolbox.calculateVIF <- function(lmFit) {
	# Run vif() within safely() to eliminate problems with "aliased" factors.  These are factors that are
	# colinear.
	
	vif <- data.frame(as.list(safely(car::vif)(lmFit)$result))
	
	if(is.null(vif)) {
		return(NULL)
	} else {
		return(tidyr::gather(vif, "term", "vif"))
	}
}

#' Compute the relative importance of a rowwise collection of linear model summaries.  Note that the
#' sum of the lmg() measures is equal to the total R-squared of the regression
#' @param lmSummaries The linear model summaries
#' @param object The object containing the summaries
Toolbox.calculateRelativeImportance <- function(lmSummaries, object) {
	# FIXIT: Properly handle the "object" argument instead of hard-coding "fit"
	return(Toolbox.applyRowWiseDF(lmSummaries, fit,  
			function(x) {
				Toolbox.logDebug("Start")
				
				if(is.null(x)) {
					df <- data.frame()
				}	
				else {
					lmg <- calc.relimp(x, type = c("lmg", "first", "last", "betasq"))
					df <- data.frame(lmg = lmg$lmg, first = lmg$first, last = lmg$last, betasq = lmg$betasq) %>% 
									tibble::rownames_to_column("term")
				}
				
				Toolbox.logDebug("End")
				
				return(df)				
			}))
}

#' Calculate the residuals from the given fit as applied to the data.  This is useful when
#' the fit is performed over a subset of the given dataframe.  Coefficients and summary
#' statistics (e.g. R2) are also returned in a standard format
#' @param fit The lm() result
#' @param data The data to which the regression is applied
#' @param identrifierVariables Variables used to identify observations in data
#' @param variables The names of the dependent variables
#' @param suffix Suffix to be applied to the variables to designate residuals
Toolbox.calcResidCoeffSummary <- function(
		fit, 
		data, 
		identifierVariables, 
		variables, 
		suffix
) {
	# TODO: Support identifierVariables
	
	if(length(variables) > 1) {
		# Left join prediction to source data

		Toolbox.logFine("Calculate prediction")

		# NOTE: The standard approach of using stats::predict() doesn't work when there are factor
		#		levels present in the out-of-sample data into which the model is being imputed. 
		#		Therefore, I need to roll a custom prediction, which isn't difficult but
		#		it needs to be careful to align the variables.
		
		vars <- rownames(coef(fit))

		# NOTE: Subset the model matrix to only include those variables that were estimated in the 
		#   	fitted universe.  This can happen, for example, when there are more levels to a factor
		#   	variable in the ex-fitted universe.   This essentially forces the assumed coeficient for
		#   	those variables to be 0.
		
		currentNAAction <- options('na.action')$na.action
		options(na.action = 'na.include')		
		modelMatrix <- model.matrix(formula(fit), data)[,vars]
		options(na.action = currentNAAction)
		
		# Force null coefficients to be 0, this arises often when there are linearly dependent dummies
		
		coefs <- coef(fit)
		coefs[is.na(coefs)] <- 0
		
		predictM <- modelMatrix %*% coefs

		Toolbox.logFine("Convert prediction to dataframe")
		
		predictDF <- data.frame(predictM)
		names(predictDF) <- variables 
		predictDF <- cbind(data.frame(index_ = data$index_, stringsAsFactors = FALSE), predictDF)
		
		# Calculate residuals by skinny-fying the two datasets, calculating prediction from
		# variable and then fattening the dataset
		
		Toolbox.logFine("Pivot narrow")
	
		# TODO: This may be more efficient if calculated at the matrix level, by column
	
		# TODO: replace reshape2::dcast() with Toolbox.pivot*()

		residuals <- reshape2::melt(data, id.vars = "index_", measure.vars = variables, value.name = "value") %>%
				dplyr::left_join(
						reshape2::melt(predictDF, id.vars = "index_", measure.vars = variables, value.name = "predicted"), 
						by = c("index_", "variable")) %>%
				dplyr::mutate(residual_ = value - predicted, variable = paste0(variable, suffix)) %>%
				reshape2::dcast(index_ ~ variable, value.var = "residual_")
		
		# glance() and tidy() don't support multiple LHS so we need to use summary() and then call 
		# those on each element of the summary list

		Toolbox.logFine("Extract coeff and summary")
	
		fitSummary <- summary(fit)
		names(fitSummary) <- variables
		
		coeffs <- dplyr::bind_rows(lapply(fitSummary, function(x) { broom::tidy(x)}), .id = "variable_id") %>%
				dplyr::rename(std_error = std.error, t_stat = statistic, p_value = p.value)   
		
		summaries <- dplyr::bind_rows(lapply(fitSummary, function(x) { broom::glance(x)}), .id = "variable_id") %>%
				dplyr::rename(r_squared = r.squared, adj_r_squared = adj.r.squared, p_value = p.value)
	}
	else {
		Toolbox.logFine("Calculate prediction")

		# NOTE: The standard approach of using stats::predict() doesn't work when there are factor
		#		levels present in the out-of-sample data into which the model is being imputed. 
		#		Therefore, I need to roll a custom prediction, which isn't difficult but
		#		it needs to be careful to align the variables.
	
		vars <- names(coef(fit))
		
		# NOTE: Subset the model matrix to only include those variables that were estimated in the 
		#   	fitted universe.  This can happen, for example, when there are more levels to a factor
		#   	variable in the ex-fitted universe.   This essentially forces the assumed coeficient for
		#   	those variables to be 0.
		
		currentNAAction <- options('na.action')$na.action
		options(na.action = 'na.include')		
		modelMatrix <- model.matrix(formula(fit), data)[,vars]
		options(na.action = currentNAAction)
		
		# Force null coefficients to be 0. This arises often when there are linearly dependent dummies
		
		coefs <- coef(fit)
		coefs[is.na(coefs)] <- 0
		
		prediction <- modelMatrix %*% coefs 
		
		residuals <- data.frame(
				data$index_, 
				data[[variables]] - prediction, 
				stringsAsFactors = FALSE
		)
		
		names(residuals) <- c("index_", paste0(variables, suffix))

		Toolbox.logFine("Extract coeff and summary")
		
		coeffs <- broom::tidy(fit) %>% dplyr::mutate(variable_id = variables) %>%
				dplyr::rename(std_error = std.error, t_stat = statistic, p_value = p.value)   
		
		summaries <- broom::glance(fit) %>% 
				dplyr::mutate(variable_id = variables) %>%
				dplyr::rename(r_squared = r.squared, adj_r_squared = adj.r.squared, p_value = p.value, 
						log_lik = logLik, aic = AIC, bic = BIC, df_residual = df.residual)
	}
	
	results <- list(
			residuals = residuals, 
			coeffs = coeffs, 
			summaries = summaries
	)
	
	return(results)	
}

Toolbox.neutralizeVariablesEach <- function(
		data, 
		identifierVariables,
		variables,
		neutVariables, 
		weightVariable,
		suffix,
		joint,
		universeVariables,
		intercept) {
	if(!is.null(universeVariables)) {
		# First fit purging regression to be a subset that includes those in any of the given universes	
		universeData <- data %>%
				dplyr::filter_at(dplyr::all_of(universeVariables), dplyr::any_vars(. == TRUE))			
	} else {
		universeData <- data
	}

	Toolbox.logFine("joint: %s", joint)
	
	if(joint) {
		universeDataSub <- universeData %>%
				# Eliminate groups for which all variables have NA for all values
				dplyr::filter_at(c(variables, weightVariable), dplyr::all_vars(any(is.finite(.))))

		if(nrow(universeDataSub) > 0) {			
			# Construct linear model
			
			# Include multiple dependent variables by listing them in cbind()
			# See https://stackoverflow.com/questions/42909238/how-to-specify-formula-in-linear-model-with-100-dependent-variables-without-havi
			
			form <- as.formula(glue::glue("cbind({paste0(variables, collapse = ',')}) ~ {paste0(neutVariables, collapse = ' + ')} {ifelse(intercept, '', '-1')}"))
			
			# Estimation regression model
	
			Toolbox.logFine("Estimate model")
			
			if(is.null(weightVariable)) {
				fit <- lm(
						form, 
						universeData, 
						na.action = "na.exclude", 
						model = FALSE, 
						x = FALSE, 
						y = FALSE
				)
			}
			else {
				fit <- lm(
						form, 
						universeData, 
						na.action = "na.exclude", 
						weights = eval(as.name(weightVariable)), 
						model = FALSE, 
						x = FALSE, 
						y = FALSE
				)				
			}
	
			Toolbox.logFine("Calculate residuals")
			
			results <- Toolbox.calcResidCoeffSummary(
					fit, 
					data, 
					identifierVariables, 
					variables, 
					suffix
			)
		} else {
			results <- NULL
		}
		
		return(results)
	} else {
		allResiduals <- data %>%
				dplyr::select_at(identifierVariables)
		
		allCoeffs <- list()
		allSummaries <- list()
		
		for(variable in variables) {
			Toolbox.logFine("variable: %s", variable)
			
			universeDataSub <- universeData %>%
					# Eliminate groups for which all variables have NA for all values
					dplyr::filter_at(c(variable, weightVariable), dplyr::all_vars(any(is.finite(.))))

			if(nrow(universeDataSub) > 0) {
				Toolbox.logFine("Neutralizing")
				
				form <- as.formula(glue::glue("{variable} ~ {paste0(neutVariables, collapse = ' + ')} {ifelse(intercept, '', '-1')}"))
			
				# Estimate regression model
				
				# TODO: Improve performance by using lm.fit() or something similar
			
				if(is.null(weightVariable)) {
					fit <- lm(
							form, 
							universeDataSub, 
							na.action = "na.exclude", 
							model = FALSE, 
							x = FALSE, 
							y = FALSE
					)
				}
				else {
					fit <- lm(
							form, 
							universeDataSub, 
							na.action = "na.exclude", 
							weights = eval(as.name(weightVariable)), 
							model = FALSE, 
							x = FALSE, 
							y = FALSE
					)				
				}
							
				fitResults <- Toolbox.calcResidCoeffSummary(
						fit, 
						data, 
						identifierVariables, 
						variable, 
						suffix
				)
				
				allResiduals <- allResiduals %>%
						dplyr::full_join(fitResults$residuals, by = identifierVariables)
				allCoeffs[[variable]] <- fitResults$coeffs
				allSummaries[[variable]] <- fitResults$summaries
			}
		}
		
		results <- list(
				residuals = allResiduals,
				coeffs = dplyr::bind_rows(allCoeffs),
				summaries = dplyr::bind_rows(allSummaries)
		)
		
		return(results)
	}	
}	


#' Neutralize the given variables by regressing them on neutVariables and returning the residuals.
#' Note that if all observations within a group for any variable are NA then all observations for that
#' group will be excluded (residuals all NA)
#' @param modelData data.frame containing the data
#' @param groupVariables Variables to be used to form groups
#' @param variables The variables to neutralize
#' @param neutVariables The variables to use in neutralizing
#' @param weightVariable The variable to use in weighting the neutralization regression
#' @param suffix The suffix to apply to the neutralized variables
#' @param joint Whether to neutralize variables jointly
#' @param universeVariables logical variables that determine the subset over which to neutralize.  The neutralization is extrapolated to securities outside the neutralization subset.
#' @param intercept Whether to include an intercept in the neutralizing regression
#' @param zeroMissings Whether missing values in neutVariables should be set to zero
#' @param addNAFactorLevel Whether factor neutVariables should have "<NA>" added as a level when it is missing 
Toolbox.neutralizeVariables <- function(
		modelData, 
		groupVariables, 
		variables, 
		neutVariables, 
		weightVariable = NULL, 
		suffix = "_res", 
		joint = TRUE, 
		universeVariables = NULL, 
		intercept = TRUE,
		zeroMissings = TRUE,
		addNAFactorLevel = TRUE
) {
	# Replace missing values of neutVariables with 0
	# Create index_ variable to facilitate joining residuals
	
	Toolbox.logDebug("groupVariables: %s", groupVariables)
	Toolbox.logDebug("variables: %s", variables)
	Toolbox.logDebug("neutVariables: %s", neutVariables)
	Toolbox.logDebug("weightVariable: %s", weightVariable)
	Toolbox.logDebug("suffix: %s", suffix)
	Toolbox.logDebug("joint: %s", joint)
	Toolbox.logDebug("universeVariables: %s", universeVariables)

	uniqueVariables <- unique(variables)
	
	Toolbox.logDebug("uniqueVariables: %s", uniqueVariables)
	
	# Add _index as a key field so that residuals can be joined in later step
	
	modelData <- modelData %>%
			dplyr::group_by_at(groupVariables) %>%
			dplyr::mutate(index_ = row_number())
	
	modelDataSub <- modelData %>%
			dplyr::ungroup() %>%
			dplyr::select(dplyr::all_of(c("index_", groupVariables, uniqueVariables, neutVariables, weightVariable, universeVariables)))

	isFactorVar <- function(x) is.factor(x) | is.character(x)		
	factorVariables <- names(modelDataSub)[sapply(modelDataSub, isFactorVar)]
	
	if(addNAFactorLevel) {
		# TODO: This will not likely work well when a factor variable has non-character levels
		
		factorNeutVariables <- dplyr::intersect(
				factorVariables,
				neutVariables
		)
		
		Toolbox.logDebug("factorNeutVariables: %s", factorNeutVariables)
		
		if(length(factorNeutVariables) > 0) {
			modelDataSub <- modelDataSub %>%
					dplyr::ungroup() %>%
					dplyr::mutate_at(
							factorNeutVariables,
							~ factor(ifelse(is.na(.), "<NA>", as.character(.)))
					)
		}		
	}
	
	if(zeroMissings) {
		nonFactorNeutVariables <- dplyr::setdiff(
				neutVariables,
				factorVariables
		)
		
		Toolbox.logDebug("nonFactorNeutVariables: %s", nonFactorNeutVariables)		
		
		modelDataSub <- modelDataSub %>%
				dplyr::ungroup() %>%
				dplyr::mutate_at(nonFactorNeutVariables, list(~ tidyr::replace_na(., 0)))				
	}
	
	allNeuts <- modelDataSub %>%
			dplyr::group_by_at(groupVariables) %>%
			tidyr::nest() %>%
			dplyr::mutate(
					neut = purrr::map(
							data, 
							~ Toolbox.neutralizeVariablesEach(
									., 
									"index_",
									uniqueVariables, 
									neutVariables, 
									weightVariable, 
									suffix, 
									joint, 
									universeVariables,
									intercept
							)
					),
					residuals = purrr::map(neut, ~ .$residuals),
					coeffs = purrr::map(neut, ~ .$coeffs),
					summaries = purrr::map(neut, ~ .$summaries)
			)
	
	allResiduals <- allNeuts %>%
			tidyr::unnest(residuals) %>%
			dplyr::select(-data, -neut, -coeffs, -summaries)
	
	allCoeffs <- allNeuts %>%
			tidyr::unnest(coeffs) %>%
			dplyr::select(-data, -neut, -residuals, -summaries)
	
	allSummaries <- allNeuts %>%
			tidyr::unnest(summaries) %>%
			dplyr::select(-data, -neut, -residuals, -coeffs)			

	Toolbox.logDebug("Join residuals")
	
	modelData <- modelData %>%			
			dplyr::left_join(allResiduals, by = c(groupVariables, "index_")) %>%
			dplyr::select(-index_) %>%
			dplyr::ungroup()

	results <- list(
			data = modelData,
			coeffs = allCoeffs,
			summaries = allSummaries
	)
	
	Toolbox.logDebug("Complete")
	
	return(results)
} 

#' Compute univariate regressions for one or more independent variables
#' @param data Model data.frame
#' @param groupVariables The variables that should be used to group regresssions (i.e. "period_dt"
#' @param dependentVariable The dependent variable in the regression
#' @param independentVariables The independent variables in the regression
#' @param weightVariable The weight variable
Toolbox.computeUnivariateRegressions <- function(
		data, 
		groupVariables,  
		dependentVariable, 
		independentVariables,  
		weightVariable = NULL, 
		includeIntercept = TRUE
) {	
	allCoeffs <- list()
	
	safeLM <- function(formula, data, weightVariable) {
		# TODO: Using lm.fit may be faster, but without the t-stats, etc.
		
		if(is.null(weightVariable)) {
			fit <- lm(formula, data = data, na.action = "na.exclude", model = FALSE)
		}
		else {
			# TODO: Not sure why I need to do this weight_ hack here
			
			fit <- lm(formula, data = data %>% dplyr::mutate(weight_ = eval(as.name(weightVariable))),
					weights = weight_, na.action = "na.exclude", model = FALSE)			
		}
		
		return(fit)
	}		
		
	dataNested <- data %>%
			dplyr::select(dplyr::all_of(c(groupVariables, dependentVariable, independentVariables, weightVariable))) %>%
			Toolbox.smartGroupBy(groupVariables) %>%
			dplyr::filter_at(
					c(dependentVariable, independentVariables, weightVariable), 
					dplyr::all_vars(any(is.finite(.)))
			) %>%
			tidyr::nest()

	# TODO: Improve this by inverting the loop to be "inside" safeLM
	
	for(indepVar in independentVariables) {
		Toolbox.logFine("indepVar: %s", indepVar)

		modelFormula <- as.formula(glue::glue("{dependentVariable} ~ {ifelse(includeIntercept, '', '0 + ')}{indepVar}"))

		Toolbox.logFine("modelFormula: %s", modelFormula)
		
		allCoeffs[[indepVar]] <- dataNested  %>%
				dplyr::mutate(
						fit = purrr::map(data, ~ safeLM(modelFormula, data = ., weightVariable)),
						tidy = purrr::map(fit, broom::tidy)
				) %>%
				tidyr::unnest(tidy) %>%
				dplyr::select(-data, -fit) %>%
				dplyr::rename(p_value = p.value, std_error = std.error)
	}
	
	return(dplyr::bind_rows(allCoeffs, .id = "variable_id"))
}

#' estimate spline model given a linear model and knots
#' the preferred model is univariate analysis
#' the searches for a solution iteratively over the starting breakpoints         
#' produces a list which contains the typical lm modeling goodies, as well as knot-specific 
#' @param model The model from lm(), such as: lm(y ~ x, data = df)
#' @param numKnots the number of knots to search over, breakpoints start with percentiles
Toolbox.calculateSplineFit <- function(model, numKnots = 1, control = seg.control()) {
    
    # Extract the x variable from the lm model object

    if(numKnots==0){
        return(model)
    } else {
		results <- segmented::segmented(
				model, 
				seg.Z = as.formula(paste0("~",names(model$model)[-1])), 
				npsi=numKnots, 
				control = control
		)
		
        return(results)
    }
    
}

#' select the optimal number of splines given a max number of knots and an improvement in r-squared
#' @seealso Toolbox.calculateSplineFit()
#' @param model The model from lm(), such as: lm(y ~ x, data = df)
#' @param numKnotsMax the maximum number of knots to search over
#' @param rsqImprov the minimum step-up in r-squared
Toolbox.calculateSplineFitOptimal <- function(model, numKnotsMax = 5, rsqImprov = .01, control = seg.control()) {    
    #instantiate df of r-squareds
    rSqDf <- data.frame(num_knots = numeric(), adj_rsq = numeric())
    
    #first row is the linear model
    rSqDf[1,] <- c(0,summary(model)$adj.r.squared)
    
    #additional rows are spline models with knots
    for(numKnots in seq(numKnotsMax)){ 
		Toolbox.logDebug("numKnots: %s", numKnots)
		
		# TODO: Cache the model so that we can return it later without re-estimating
		
		rSq <- (Toolbox.calculateSplineFit(model, numKnots = numKnots, control) %>% summary())["adj.r.squared"]

		Toolbox.logDebug("rSq: %s", rSq)
		
        rSqDf[(numKnots+1),] <- c(numKnots, rSq)        
    }

    #calculate changes in rsq
	
    rSqDf <- rSqDf %>%
        dplyr::mutate(chg_adj_rsq = adj_rsq - lag(adj_rsq)) %>%
        dplyr::mutate(chg_adj_rsq = ifelse(is.na(chg_adj_rsq), adj_rsq, chg_adj_rsq))
    
    #define the optimal number of knots, which satisfies the rsq hurdle
	
    optKnots <- max((rSqDf %>% dplyr::filter(chg_adj_rsq > rsqImprov))$num_knots)
    
	Toolbox.logDebug("optKnots: %s", optKnots, control)
	
	# Return the appropriate spline
	
    results <- Toolbox.calculateSplineFit(model, optKnots)
	
	return(results)
}
