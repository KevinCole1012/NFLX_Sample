# Perform industry/country decomposition of security-level exposure/return data
# 
###############################################################################

#' Compute industry/country decomposition for the given data.frame
Toolbox.calcCountryIndustryDecomposition <- function(
		decompData, 
		decompVars = c("beta_loc", "beta_usa"), 
		returnVars = c("totalreturn_loc", "totalreturn_usd"),
		weightVar = "weight", 
		identifierVars = NULL,
		industryVar = "risk_industry",
		countryVar = "risk_country"
) {
	# NOTE: See %IMM_ComputeGroupVar() for more details on the estimation of "group" exposures
	
	Toolbox.logFine("decompVars: %s", decompVars)
	Toolbox.logFine("returnVars: %s", returnVars)
	Toolbox.logFine("weightVariable: %s", weightVar)
	
	Toolbox.logFine("industryVar: %s", industryVar)
	Toolbox.logFine("countryVar: %s", countryVar)
	
	if(is.null(identifierVars)) {
		identifierVars <- c("identifier", "exchange_country")
	}
	
	# Force to be a character variable since factors don't play well
	
	decompVars <- as.character(decompVars)
	
	if(is.null(weightVar)) {
		weightVar <- "one_"
		
		decompData <- decompData %>%
				dplyr::mutate(one_ = 1)
	}

	# Eliminate data that has any missing decomposition variables, countries or industries
	# Replace missing returns with 0
	# Force weights to sum to 1
	
	# NOTE: If any of decompVars is missing then the observation is dropped.   This can be a problem with sparse variables.

	varMap <- c(
			identifierVars,
			decompVars,
			returnVars,
			"risk_industry_" = industryVar, 
			"risk_country_" = countryVar,
			"weight__"
	)
	
	decompDataND <- decompData %>%
			dplyr::ungroup() %>%
			dplyr::distinct_at(identifierVars, .keep_all = TRUE) %>%
			dplyr::filter_at(
					dplyr::all_of(c(decompVars, countryVar, industryVar, weightVar)), 
					dplyr::all_vars(!is.na(.))
			) %>%
			dplyr::mutate(weight__ = eval(as.name(weightVar)) / sum(eval(as.name(weightVar)), na.rm = TRUE)) %>%
			dplyr::select_at(varMap)	
		
	if(!is.null(returnVars)) {
		Toolbox.logFine("Set missing returns to 0")
		
		decompDataND <- decompDataND %>%
				dplyr::mutate(
						dplyr::across(
								dplyr::all_of(returnVars), 
								~ tidyr::replace_na(., 0)
						)
				)
	}
	
	if(dim(decompDataND)[1] > 0) {	
		Toolbox.logFine("Construct countryWeights")
		
		countryWeights <- decompDataND %>%
				dplyr::group_by(risk_country_) %>%
				dplyr::summarize(weight_country = sum(weight__, na.rm = TRUE)) %>%
				dplyr::mutate(
						weight_country = weight_country / sum(weight_country), 
						weight_country_inv = 1 / weight_country
				) %>%
				dplyr::filter(weight_country > 0)

		Toolbox.logFine("Construct industryWeights")
		
		industryWeights <- decompDataND %>%
				dplyr::group_by(risk_industry_) %>%
				dplyr::summarize(weight_industry = sum(weight__, na.rm = TRUE)) %>%
				dplyr::mutate(
						weight_industry = weight_industry / sum(weight_industry), 
						weight_industry_inv = 1 / weight_industry
				) %>%
				dplyr::filter(weight_industry > 0)
		
		Toolbox.logFine("Construct decompData")
		
		decompDataND <- decompDataND %>%
				dplyr::left_join(industryWeights, by = "risk_industry_") %>%
				dplyr::left_join(countryWeights, by = "risk_country_") %>%
				dplyr::mutate(
						risk_country_ = as.factor(risk_country_), 
						risk_industry_ = as.factor(risk_industry_)
				)
		
		Toolbox.logFine("Construct 'dummies'")	
		
		decompDataND <- decompDataND %>%
				dplyr::left_join(
						decompDataND %>%
								dplyr::select_at(c(identifierVars, "risk_industry_", "weight_industry_inv")) %>%
								dplyr::mutate(risk_industry_ = paste0("industry_", risk_industry_)) %>%
								tidyr::spread(risk_industry_, weight_industry_inv, fill = 0),
						by = identifierVars
				) %>%
				dplyr::left_join(
						decompDataND %>%
								dplyr::select_at(c(identifierVars, "risk_country_", "weight_country_inv")) %>%
								dplyr::mutate(risk_country_ = paste0("country_", risk_country_)) %>%
								tidyr::spread(risk_country_, weight_country_inv, fill = 0),
						by = identifierVars
				) %>%
				dplyr::mutate(one = 1) 	
		
		Toolbox.logFine("Build matrices")	
		
		riskIndustries <- unique(paste0("industry_", industryWeights$risk_industry_))
		riskCountries <- unique(paste0("country_", countryWeights$risk_country_))
		
		X <- as.matrix(subset(decompDataND, select = c("one", riskCountries, riskIndustries)))
		
		R <- rep(c(0, 1, 0, 0, 0, 1), 
				c(
						1, length(riskCountries), length(riskIndustries), 
						1, length(riskCountries), length(riskIndustries)
				)
		)
		
		dim(R) <- c(length(R) / 2, 2)
		
		allVars <- c(decompVars, returnVars)
		
		betaReturn <- as.matrix(subset(decompDataND, select = allVars))
		
		Toolbox.logFine("Build XWX")	
		
		W <- decompDataND$weight__
		dim(W) <- c(1, length(W))
		
		XW <- sweep(t(X), MARGIN = 2, W, `*`)
		XWX <- XW %*% X
		
		zero22 <- rep(0, 4)
		dim(zero22) <- c(2, 2)
		
		zero24 <- rep(0, 2 * length(allVars))
		dim(zero24) <- c(2, length(allVars))
		
		Toolbox.logFine("Build A and B")	
		
		A <- cbind(rbind(XWX, t(R)), rbind(R, zero22))
		
		if(is.finite(determinant(A)$modulus)) {
			B <- rbind(XW %*% betaReturn, zero24)
			
			Toolbox.logFine("Compute decomposition")	
			
			# NOTE: Use MASS::ginv() instead of solve() because in some periods (early) the country and industry dummies can
			#         be linearly dependent due to sparse coverage of the variable to be decomposed
			
			betaReturnDecomp <- ((MASS::ginv(A) %*% B) / c(1, countryWeights$weight_country, industryWeights$weight_industry, 1, 1))[1 : (1 + length(riskCountries) + length(riskIndustries)),]
			
			betaDF <- data.frame(
					name = c("market", riskCountries, riskIndustries), 
					type = c("market", rep("country", length(riskCountries)), rep("industry", length(riskIndustries))),
					weight = c(1, countryWeights$weight_country, industryWeights$weight_industry),
					betaReturnDecomp
			)
			
			if(is.null(colnames(betaReturnDecomp))) {
				# Hack to handle case of single model variable
				names(betaDF)[4] <- allVars
			}
			
			# Convert X, whose values are the relative weight of an industry/country, into a matrix of 0/1 values
			# and multiply by beta decomposition to get "predicted" beta.  Subtract this from the betas to get the residual beta.
			
			dummyX <- ifelse(X > 0, 1, 0)
			
			residuals <- betaReturn - dummyX %*% betaReturnDecomp
			
			betaResid <- data.frame(
					decompDataND %>% 
							dplyr::select_at(c(identifierVars, "risk_industry_", "risk_country_")) %>%
							dplyr::rename_at("risk_industry_", ~ industryVar) %>%
							dplyr::rename_at("risk_country_", ~ countryVar),
					residuals
			)
			
			rSquared <- 1 - colSums(residuals ^ 2, na.rm = TRUE) / colSums(sweep(betaReturn, MARGIN = 2, STATS = colMeans(betaReturn, na.rm = TRUE)) ^ 2, na.rm = TRUE) 
			
			decompSummary <- data.frame(rSquared) %>%
					dplyr::rename(r_squared = rSquared) %>%
					tibble::rownames_to_column("variable_id")
			
			result <- list(
					decomposition = betaDF, 
					residual = betaResid, 
					summary = decompSummary
			)
		}
		else {
			Toolbox.logError("A is singular")	
			
			# A is singular and therefore solve(A) will fail, so return a missing data.frame
			
			result = list(
					decomposition = data.frame(), 
					residual = data.frame(), 
					summary = data.frame()
			)			
		}
	}
	else {
		result = list(
				decomposition = data.frame(), 
				residual = data.frame(), 
				summary = data.frame()
		)
	}
	
	return(result)
}

#' Decompose betas by industry and country in each period.
Toolbox.calcCountryIndustryDecompByPeriod <- function(
		decompData, 
		decompVars = c("beta_loc", "beta_usa"), 
		returnVars = c("totalreturn_loc", "totalreturn_usd"),
		weightVar = "weight", 
		identifierVars = NULL,
		dateVar = "period_dt",
		industryVar = "risk_industry",
		countryVar = "risk_country"
) {
	Toolbox.logFine("decompVars: %s", decompVars)
	Toolbox.logFine("returnVars: %s", returnVars)
	Toolbox.logFine("weightVariable: %s", weightVar)
	
	Toolbox.logFine("dateVar: %s", dateVar)
	Toolbox.logFine("industryVar: %s", industryVar)
	Toolbox.logFine("countryVar: %s", countryVar)
	
	Toolbox.logDebug("Calculate decomposition")
		
	betaDecomp <- decompData %>%
			dplyr::group_by_at(dateVar) %>%
			tidyr::nest() %>%
			dplyr::mutate(
					decompAll = purrr::map(
							data,
							~ Toolbox.calcCountryIndustryDecomposition(
									., 
									decompVars, 
									returnVars, 
									weightVar, 
									identifierVars,
									industryVar,
									countryVar
							)
					),
					decomposition = lapply(decompAll, function(x) x$decomposition),
					residual = lapply(decompAll, function(x) x$residual),
					summary = lapply(decompAll, function(x) x$summary)
			)
	
	Toolbox.logDebug("Extract decomposition")
	
	decomposition <- betaDecomp %>%
			dplyr::select_at(c(dateVar, "decomposition")) %>%
			tidyr::unnest(cols = decomposition) %>%
			dplyr::ungroup()
	
	residual <- betaDecomp %>%
			dplyr::select_at(c(dateVar, "residual")) %>%
			tidyr::unnest(cols = residual) %>%
			dplyr::ungroup()
	
	summary <- betaDecomp %>%
			dplyr::select_at(c(dateVar, "summary")) %>%
			tidyr::unnest(cols = summary) %>%
			dplyr::ungroup()
	
	results <- list(
			decomposition = decomposition, 
			residual = residual, 
			summary = summary
	)
	
	return(results)
}

#' Merge the results of Toolbox.calcCountryIndustryDecompByPeriod() into a single "side" data.frame 
Toolbox.mergeCountryIndustryDecomp <- function(
		decompData, 
		decompVars, 
		returnVars,
		weightVar = "weight",
		identifierVars = NULL,
		dateVar = "period_dt",
		industryVar = "risk_industry",
		countryVar = "risk_country",
		wide = TRUE
) {	
	if(wide) {
		decompDataMarket <- decompData$decomposition %>%
				dplyr::filter(type == "market") %>%
				dplyr::select(-name, -type) %>%
				dplyr::rename_at(dplyr::all_of(c(weightVar, decompVars, returnVars)), 
						list( ~paste0(., "_market")))
		
		decompDataIndustry <- decompData$decomposition %>%
				dplyr::filter(type == "industry") %>%
				dplyr::mutate(risk_industry_ = substr(name, 10, nchar(as.character(name)))) %>%
				dplyr::rename_at("risk_industry_", ~ industryVar) %>%
				dplyr::select(-name, -type) %>%
				dplyr::rename_at(dplyr::all_of(c(weightVar, decompVars, returnVars)), 
						list( ~paste0(., "_industry")))
		
		decompDataCountry <- decompData$decomposition %>%
				dplyr::filter(type == "country") %>%
				dplyr::mutate(risk_country_ = substr(name, 9, nchar(as.character(name)))) %>%
				dplyr::rename_at("risk_country_", ~ countryVar) %>%
				dplyr::select(-name, -type) %>%
				dplyr::rename_at(dplyr::all_of(c(weightVar, decompVars, returnVars)), 
						list( ~paste0(., "_country")))
		
		decompDataResidual <- decompData$residual %>%
				dplyr::rename_at(dplyr::all_of(c(decompVars, returnVars)), 
						list( ~paste0(., "_residual")))

		if(is.null(dateVar)) {
			decompDataMerged <- decompDataResidual %>%
					dplyr::cross_join(decompDataMarket) %>%
					dplyr::left_join(decompDataIndustry, by = c(industryVar)) %>%
					dplyr::left_join(decompDataCountry, by = c(countryVar))			
		} else {
			decompDataMerged <- decompDataResidual %>%
					dplyr::left_join(decompDataMarket, by = dateVar) %>%
					dplyr::left_join(decompDataIndustry, by = c(dateVar, industryVar)) %>%
					dplyr::left_join(decompDataCountry, by = c(dateVar, countryVar))			
		}
	} else {
		# NOTE: All of the variables here are called *_value to avoid a collision between
		#		that column and the industryVar or countryVar
		
		# TODO: Support weightVar
		
		decompDataResidualT <- decompData$residual %>%
				dplyr::group_by_at(c(dateVar, identifierVars, industryVar, countryVar)) %>%
				Toolbox.pivotNarrow("variable_id", "residual_decomp")

		decompDataMarketT <- decompData$decomposition %>%
				dplyr::filter(type == "market") %>%
				dplyr::select(-type, -weight, -name) %>%
				dplyr::group_by_at(dateVar) %>%
				Toolbox.pivotNarrow("variable_id", "market_decomp")
				
		decompDataIndustryT <- decompData$decomposition %>%
				dplyr::filter(type == "industry") %>%
				dplyr::mutate(risk_industry_ = substr(name, 10, nchar(as.character(name)))) %>%
				dplyr::rename_at("risk_industry_", ~ industryVar) %>%
				dplyr::select(-type, -weight, -name) %>%
				dplyr::group_by_at(c(dateVar, industryVar)) %>%
				Toolbox.pivotNarrow("variable_id", "industry_decomp")
		
		decompDataCountryT <- decompData$decomposition %>%
				dplyr::filter(type == "country") %>%
				dplyr::mutate(risk_country_ = substr(name, 9, nchar(as.character(name)))) %>%
				dplyr::rename_at("risk_country_", ~ countryVar) %>%
				dplyr::select(-type, -weight, -name) %>%
				dplyr::group_by_at(c(dateVar, countryVar)) %>%
				Toolbox.pivotNarrow("variable_id", "country_decomp")
		
		decompDataMerged <- decompDataResidualT %>%
				dplyr::left_join(decompDataMarketT, by = c(dateVar, "variable_id")) %>%
				dplyr::left_join(decompDataIndustryT, by = c(dateVar, industryVar, "variable_id")) %>%
				dplyr::left_join(decompDataCountryT, by = c(dateVar, countryVar, "variable_id"))
	}
	
	return(decompDataMerged)	
}
