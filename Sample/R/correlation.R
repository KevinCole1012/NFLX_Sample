# Functions for computing correlations
# 
###############################################################################

#' Enhances the standard cor() function to check for the existance of a minimum number of finite observations
Toolbox.cor <- function(
		x, 
		y, 
		use = "everything",
		method = c("pearson", "kendall", "spearman"), 
		minObservations = 0
) {
	if(pmin(sum(is.finite(x) & is.finite(y)), na.rm = TRUE) >= minObservations) {
		sdX <- stats::sd(x, na.rm = TRUE)
		sdY <- stats::sd(y, na.rm = TRUE)
			
		if(is.na(sdX) || is.na(sdY) || sdX == 0 || sdY == 0) {
			return(NA_real_)
		} else {
			return(cor(x, y, use, method))
		}
	} else {
		return(NA_real_)
	}
}

#' Calculate weighted correlation of the two given variables
Toolbox.corrWeighted <- function(
		x,
		y,
		wgt
) {
	if(is.vector(x) && is.vector(y) && is.vector(wgt)) {
		# Subset to values that are non-NA for all x, y and wgt
		nonNA <- is.finite(x) & is.finite(y) & is.finite(wgt) & (dplyr::coalesce(wgt, 0) > 0)
		
		results <- stats::cov.wt(cbind(x[nonNA], y[nonNA]), wgt[nonNA], cor = TRUE)
		
		return(results$cor[1,2])		
	} else {
		Toolbox.logError("x, y and wgt need to be vectors")
		stop("x, y and wgt need to be vectors")
	}	
}

#' Compute correlations of the given variables by groupingVariable
#' @param data The data.frame containing the variables to correlated
#' @param groupingVariables The variables that should be used to create groups
#' @param includeVariables The variables to be correlated (vector of characters)
#' @param ... Parameters passed to cor() function
Toolbox.calcGroupedCorrelation <- function(
		data, 
		groupingVariables, 
		variables = NULL, 
		withVariables = NULL, 
		wide = FALSE, 
		excludeDiag = TRUE, 
		...
) {	
	if(is.null(variables)) {
		numericVariables <- names(data)[unlist(lapply(data, is.numeric))] 
		variables <- dplyr::setdiff(numericVariables, groupingVariables)
	}

	if(is.null(withVariables)) {
		groupedCorr <- data %>%
				Toolbox.smartGroupBy(groupingVariables) %>%
				dplyr::do(
						data.frame(cor(dplyr::select_at(., variables), ...)) %>%
								tibble::rownames_to_column("variable_")
				)
	} else {
		groupedCorr <- data %>%
				Toolbox.smartGroupBy(groupingVariables) %>%
				dplyr::do(
						data.frame(cor(dplyr::select_at(., variables), dplyr::select_at(., withVariables), ...)) %>%
								tibble::rownames_to_column("variable_")
				)
	}
	
	if(!wide) {
		# Skinny data is preferred so use gather
		
		groupedCorr <- groupedCorr %>%
				Toolbox.smartGroupBy(c(groupingVariables, "variable_")) %>%
				tidyr::gather(variable_col, correlation, -groupingVariables, -variable_) %>%
				dplyr::rename(variable_row = variable_) %>%
				dplyr::filter(!excludeDiag | variable_col != variable_row)		
	}

	return(groupedCorr)
}

#' Computes the "cross autocorrelation" which the the cross-sectional correlation of a variable with a lagged
#' version of itself.  For example, the correlation of a factor at time T with the same factor at time T - 1
#' at the security level.  This gives a measure of the turnover of the factor.
#' NOTE: It's that there be no duplicates within date/identifierVariables so that the lag is properly
#' calculated.
#' @param data The data.frame containing the variables to tbe correlated
#' @param dateVariable The name of the variable that gives the temporal ordering of the data 
#' @param identifierVariables The identifiers that are uses to match variable values across time
#' @param variables The variables to correlate
#' @param withVariables The variables to correlate with.  If NULL then defaults to variables argument.
#' @param lags The lag time lengths
Toolbox.calcCrossAutocorrelation <- function(
		data, 
		dateVariable, 
		identifierVariables, 
		variables, 
		withVariable = NULL, 
		lags = 1
) {
	
	# If withVariable is NULL then use "value_" which is each variable itself, making this into an "autocorrelation".
	# Othewise the lead/lag is performed against a fixed variable: withVariable
	
	if(is.null(withVariable)) {
		withVariable = "value_"
	}
	
	dataT <- data %>%
		dplyr::select_at(c(dateVariable, identifierVariables, variables)) %>%
		dplyr::group_by_at(c(dateVariable, identifierVariables)) %>%
		Toolbox.pivotNarrow("variable_id", "value_", variables)

	if(withVariable != "value_") {
		dataT <- dataT %>%
				dplyr::left_join(data %>%
								dplyr::ungroup() %>%
								dplyr::select_at(c(dateVariable, identifierVariables, withVariable)),
						by = c(dateVariable, identifierVariables)
				)
	}
	
	dataT <- dataT  %>%
			dplyr::group_by_at(c(identifierVariables, "variable_id")) %>%
			dplyr::arrange_at(dateVariable)
	
	autoCorrs <- Toolbox.createList(length(lags))
	names(autoCorrs) <- as.character(lags)
	
	# TODO: Rather than lag the value, it may be best to lead the dateVariable and then join on that in order to avoid 
	#			problems related to duplicate values for values for a given identifierVariables
	
	for(lagIndex in seq_along(lags)) {
		lag <- lags[lagIndex]
		
		Toolbox.logFine("lag: %s", lag)
		
		autoCorrs[[lagIndex]] <- dataT %>%
			dplyr::mutate(value_lag_ = Toolbox.smartLag(value_, lag)) %>%
			dplyr::group_by_at(c("variable_id", dateVariable)) %>%
			dplyr::summarize(
					pearson = cor(eval(as.name(withVariable)), value_lag_, use = "pairwise.complete.obs", method = "pearson"),
					spearman = cor(eval(as.name(withVariable)), value_lag_, use = "pairwise.complete.obs", method = "spearman"),
					gamma = cov(eval(as.name(withVariable)), value_lag_, use = "pairwise.complete") / var(value_lag_, na.rm = TRUE),
					n = sum(is.finite(eval(as.name(withVariable))) * is.finite(value_lag_), na.rm = TRUE)
			) %>%
			dplyr::ungroup()
	}
	
	return(
			dplyr::bind_rows(autoCorrs, .id = "lag") %>%
				dplyr::mutate(lag = as.numeric(lag))
	)		
}

#' Calculates the chi-squared test statistic for the two given variables 
Toolbox.calcChiSquareP <- function(x, y) {
	xF <- factor(x)
	yF <- factor(y)
	
	if(length(levels(xF)) > 1 && length(levels(yF)) > 1) { 
		return(chisq.test(xF, yF)$p.value)
	} else {
		return(NA)
	}
}

#' returns a n x 1 matrix of cosine similarities, 
#' which is a correlation-type metric, which often arises in text parsing
#' @param matrix1 n x k
#' @param vector2 1 x k
#' @seealso https://en.wikipedia.org/wiki/Cosine_similarity
Toolbox.calcCosineSimilarityVector <- function(matrix1, vector2=NULL) {
	matrix1 %*% t(vector2) / (sqrt(rowSums(matrix1 ^ 2)) * sqrt(rowSums(vector2 ^ 2)))   
}


#' returns a n1 x n2 matrix of cosine similarities,
#' if matrix2 is missing, returns a n x n matrix of cosine similarities,  
#' which is a correlation-type metric, which often arises in text parsing
#' @param matrix1 w_words x n1_firms
#' @param matrix2 w_words x n2_firms
#' @seealso https://en.wikipedia.org/wiki/Cosine_similarity
Toolbox.calcCosineSimilarity <- function(matrix1, matrix2=NULL, termLabels = FALSE) {
	
	if (is.null(matrix2)){
		#much much faster
		cosSim <- crossprod(matrix1)/sqrt(tcrossprod(colSums(matrix1^2)))
	} else {
		cosSim <- lapply(seq(1:nrow(matrix2)), function(x){
		Toolbox.calcCosineSimilarityVector(matrix1,matrix2[x,] %>% t())
	}) %>% 
		bind_cols() %>% 
		as.data.frame() }

	if (termLabels) {
		#headers
		cosSim <- cosSim %>% 
			setNames(names(matrix2))
	}

	return(cosSim)

}


#' helper function which allows rapid call and convenient formatting of self-referential cosine similarities
#' @param textDictionary n x k dataframe of all n terms ("text") and k column measures or documents
#' @param terms 1 x j list of terms which should be contained within textDictionary
Toolbox.calcCosineSimilarityElements <- function(textDictionary, terms){
  
	#matrix formatting to get the transpositioning to format correctly
	textDf <- as.data.frame(textDictionary)

	rownames(textDf) <- textDictionary$text

	textDf <- subset(textDf, select = -c(text) ) %>% as.matrix()

	cosSimlAll <- lapply(terms, function(x){

		termMatrix <- textDictionary %>% filter(text %in% x) %>% dplyr::select(-text) %>% as.matrix()

		cosSiml <- Toolbox.calcCosineSimilarity(textDf,termMatrix) %>% as.data.frame() %>% setNames(x)

		cosSiml$text <- rownames(cosSiml) 
		#matrix formatting
		rownames(cosSiml) <- NULL

		return(cosSiml)

	}) %>% 
	purrr::reduce(dplyr::full_join, by = c("text")) 

	#reorder columns
	cosSimlAll <- cosSimlAll[,c("text", terms)]

	return(cosSimlAll)
  
}



#' condense a correlation matrix by interactively collapsing highest correlated firms
#' this has to be run interativly, and so is rather slow to run.
#' @param corrMat n x n correlation matrix (formatted as as.matrix, not as.data.frame)
#' @param nGroups the number of groups we wish to end up with
#' @param cleanGroupIds pivot security-wise the groups constitutents.  This results in a larger df.
Toolbox.condenseCorrMatrix <- function(corrMat, nGroups, cleanGroupIds=T){

	#enforce cor matrix to be upper triangular.  also kills the main diagonal
	#cuts in half the number of rows & eliminates the double-time covariance shuffle
	corrMat[upper.tri(corrMat)] <- NA

	#kill diagonal
	diag(corrMat) <- NA
	
	#eliminate zeros
	corrMat[corrMat == 0 ] <- NA

	#transpose.  Prepend an underscore to avoid misbehaved numeric names
	corrMatDf <- data.frame(corrMat) %>% 
		setNames(paste0("_",colnames(corrMat))) %>%
		dplyr::mutate(variable1 = paste0("_",colnames(corrMat))) %>%
		tidyr::gather("variable2","corr",-variable1) %>%
		dplyr::filter(corr != 1, !is.na(corr), variable1 != variable2) %>%
		#add a fuzz value to avoid correlations with exactly the same value
		dplyr::mutate(corr = corr + (.Machine$double.eps * dplyr::row_number())) %>%   
		data.table::setDT()
	

	while(length(unique(c(unique(corrMatDf$variable1),unique(corrMatDf$variable2)))) > nGroups) {
		
		Toolbox.logDebug("nGroups = %s", length(unique(c(unique(corrMatDf$variable1),unique(corrMatDf$variable2)))))
		
		#grab ids associated with the largest correlation
		newIds <- c(as.matrix(corrMatDf[corr == max(corr)][,-3]))
		newIdsStack = paste0(sort(newIds),collapse="|")
		
		#slightly faster data.table version
		
		corrMatDfMod <- corrMatDf %>%
			# subset the corr matrix to rows/columns associated with max correlation
			.[variable1 %in% newIds | variable2 %in% newIds] %>%
			# eliminate the value itself
			.[!(variable1 %in% newIds & variable2 %in% newIds)] %>%
			# define var1 as the high-correlation securities
			.[, var1 := ifelse(variable1 %in% newIds,variable1,variable2)] %>%
			# define var2 as all others
			.[, var2 := ifelse(variable1 == var1,variable2,variable1)] %>%
			# average all correlations within the max groups
			.[, .(corr=mean(corr)), by=var2]  %>%
			# define new grouping variables.  Separated by |
			.[, var1 := newIdsStack] %>%
			data.table::setnames(c("variable1","corr","variable2")) %>%
			data.table::setcolorder(c("variable1","variable2","corr"))
		
		#create new cor matrix by stacking modified rows onto an anti-joined existing correlations
		corrMatDf <- data.table::rbindlist(list(corrMatDfMod,
												corrMatDf[!variable1 %in% newIds][!variable2 %in% newIds]))
	}
	
	#remove prepended "_"
	corrMatDf <- corrMatDf %>% 
		dplyr::as_tibble() %>%
		dplyr::mutate(variable1 = stringr::str_remove_all(variable1,"\\_"), variable2 = stringr::str_remove_all(variable2,"\\_"))
	
	if (cleanGroupIds) {
		
		#giant string split and transpose for all groups
		groupIds <- data.frame(stringr::str_split(unique(unique(corrMatDf$variable1),unique(corrMatDf$variable2)),"\\|",simplify=T),
								group = unique(unique(corrMatDf$variable1),unique(corrMatDf$variable2))) %>% 
			dplyr::mutate(group_id = row_number()) %>%
			tidyr::gather("variable","identifier",-group,-group_id) %>%
			dplyr::filter(!is.na(identifier), identifier !="") %>%
			dplyr::select(-"variable") 
		# option to count number of obs each group
		# %>% dplyr::group_by(group_id) %>% dplyr::mutate(n=n()) %>% dplyr::ungroup()
		
		corrMatDf <- corrMatDf %>% 
			dplyr::inner_join(groupIds %>% rename(variable1 = group, group_id1 = group_id, identifier1 = identifier)) %>%
			dplyr::inner_join(groupIds %>% rename(variable2 = group, group_id2 = group_id, identifier2 = identifier))
		
	}
	
	return(corrMatDf)
	
	# #dplyr version -- easier to read but slower
	# corrMatDf <- corrMatDf %>% 
	#	 #try to speed up a tad
	#	 dtplyr::lazy_dt() %>% 
	#	 #subset the corr matrix to rows/columns associated with max  correlation
	#	 filter(variable1 %in% newIds | variable2 %in% newIds) %>%
	#	 #eliminate the value itself
	#	 filter(!(variable1 %in% newIds & variable2 %in% newIds)) %>%
	#	 #define var1 as the high-correlation securities
	#	 mutate(var1 = ifelse(variable1 %in% newIds,variable1,variable2)) %>%
	#	 #define var2 as all others
	#	 mutate(var2 = ifelse(variable1==var1,variable2,variable1)) %>%
	#	 #average all correlatoins within the max groups
	#	 group_by(var2) %>%
	#	 summarise(corr=mean(corr)) %>%
	#	 #define new grouping variables.  Separated by |
	#	 mutate(var1 = newIdsStack) %>%
	#	 #rename and stack with now-max-variable-free original data
	#	 rename(variable1=var1, variable2=var2) %>%
	#	 as_tibble() %>%
	#	 bind_rows(corrMatDf %>% filter(!variable1 %in% newIds, !variable2 %in% newIds))
	
	
}


#' given a correlation matrix, assign peers based on correlations.
#' peers are based either on a floor or # of names (ranked by correlation)
#' if both inputs are >0, this is an OR condition, not an and condition
#' to get the _and_ condition, run this separately and inner_join
#' @param corrMat n x n correlation matrix (formatted as as.matrix, not as.data.frame)
#' @param minCorr assign peers for all security-wise correlations above this value
#' @param maxPeers assign top number of peers
Toolbox.peerifyCorrMatrix <- function(corrMat, minCorr=0, maxPeers=0) {
	
	#remove cor == 1 items
	corrMatNoOnes <- (corrMat == 1)
	corrMatNoOnes[corrMatNoOnes==TRUE] <- 1
	corrMatNoOnes <- -corrMatNoOnes + 1
	corrMatNoOnes <- corrMat * corrMatNoOnes

	corrMatNoOnesDf <- data.frame(corrMatNoOnes) %>% 
		#dance to get the colnames to admit numeric names
		setNames(paste0("_",colnames(corrMatNoOnes))) %>%
		dplyr::mutate(identifier = paste0("_",colnames(corrMatNoOnes))) %>%
		tidyr::gather("peer","corr",-identifier) %>% 
		dplyr::filter(identifier != peer)
	
	if (minCorr>0){
		
		corrPeersMinCorr <- corrMatNoOnesDf %>%
			dplyr::filter(corr > minCorr) %>%
			dplyr::distinct(identifier,peer, .keep_all = TRUE) %>%
			dplyr::group_by(identifier) %>%
			dplyr::mutate(n_peers = n()) %>%
			dplyr::ungroup() %>%
			dplyr::mutate(identifier = stringr::str_remove(identifier,"\\_"),
							peer = stringr::str_remove(peer,"\\_")) %>%
			dplyr::arrange(identifier, peer) 
		
	}
	
	if (maxPeers>0) {
		
		corrPeersmaxPeers <- corrMatNoOnesDf %>% 
			dplyr::group_by(identifier) %>%
			dplyr::slice_max(corr, n = maxPeers) %>%
			dplyr::mutate(n_peers = n()) %>%
			dplyr::ungroup() %>%
			dplyr::mutate(identifier = stringr::str_remove(identifier,"\\_"),
							peer = stringr::str_remove(peer,"\\_")) %>%
			dplyr::arrange(identifier, peer)
		
	}
	
	if (maxPeers>0 & minCorr >0) {
		
		corrPeersBoth <- corrPeersMinCorr %>% 
			dplyr::filter(n_peers >= maxPeers) %>%
			dplyr::bind_rows(corrPeersmaxPeers %>% 
								dplyr::filter(identifier %in% 
									(corrPeersMinCorr %>% 
										dplyr::filter(n_peers < maxPeers) %>%  
										dplyr::select(identifier) %>% 
										unique() %>% 
										as.matrix() %>% 
										c()
									)
								)
							)
		return(corrPeersBoth)
		
	}
	
	if (minCorr>0 & maxPeers==0){
		return(corrPeersMinCorr)   
	}
	
	if (minCorr==0 & maxPeers>0){
		return(corrPeersmaxPeers)   
	}
	
	
}


#' wrapper around the glasso::glasso package which estimates sparse precision matrices 
#' @param retMatrix t_x_n matrix of the returns the glasso is using
#' @param rho the regularization intensity parameter; referred to as "lambda" in lasso
#' @param nobs number of observations 
#' @param penalizeDiagonal whether to penalize the diagonal elements
#' @param mbApproximation this approximation enables the process to run much faster
#' @param zeroConstraints inherited transformation.  Default is inherited to be "averageOverrideZero"
#' @param returnInvCM whether to return the precision matrix (default) or the covariance itself
#' @param initialValueInvCm initial values of the precision matrix.  Can help with convergence
Toolbox.glassoFit <- function(
		retMatrix, 
		rho, 
		nobs = nrow(retMatrix), 
		penalizeDiagonal = FALSE, 
		mbApproximation = TRUE,
		zeroConstraints = NULL,
		returnInvCM = TRUE,
		initialValueInvCm = NULL ) {
  
	start <- "cold"
	initialValueCm <- initialValueInvCm

	if(!is.null(initialValueInvCm)){
		start <- "warm"
		initialValueCm <- solve(Toolbox.positiveDefinitify(initialValueCm))
	}

	glassoFit <- glasso::glasso(s = var(retMatrix),
					rho = rho,
					nobs = nobs,
					penalize.diagonal = penalizeDiagonal,
					approx = mbApproximation,
					zero = zeroConstraints,
					start = start,
					wi.init = initialValueInvCm,
					w.init = initialValueCm)

	if(returnInvCM){
		glassoFitCm <- glassoFit$wi
	} else {
		glassoFitCm <- glassoFit$w
	}

	#return inverse covariance matrix
	return(glassoFitCm)
	
}

#tune rho parameter in glasso via bisection.  Iterate until the median firm number of estimated elements is above a threshold 
#run over a subset of columns to speed up the process
#' @param retMatrix t_x_n matrix of the returns the glasso is using
#' @param medianThreshold the threshold at which the loop stops
#' @param crossValidateGroups number of groups to cross-validate over
#' @param convergenceStep iteration speed.   A value closer to 1 will cause this to run slower.  Bisection = .5
#' @param crossValidateSetSeed whether to set seeds in the cross validation.  Handy for tying out results, since glasso has a stochastic start
#' @param ... arguments to pass through to Toolbox.glassoFit
Toolbox.glassoTuneRho <- function(
		retMatrix, 
		medianThreshold = 20, 
		crossValidateGroups = 10,
		convergenceStep = .75,
		crossValidateSetSeed = TRUE,
		...) {
  
  #instantiate rho with upper bound of 1
	currentRho <- 1

	cvConverge <- FALSE

	while (!cvConverge){
	
	Toolbox.logFine("current Rho parameter is: %s", currentRho)

	cvList <- list(crossValidateGroups)

	for (i in 1:crossValidateGroups){
	
		if(crossValidateSetSeed) {
			set.seed(i)
		}

		#set the cols grabbed equal to cols / crossValidateGroups

		sampleCols <- sample(ncol(retMatrix), round(ncol(retMatrix) / crossValidateGroups))

		cvList[[i]] <- Toolbox.glassoFit(retMatrix[,sampleCols], rho = currentRho, ...) %>% 
			Toolbox.countNonzeroMatrix() %>% 
			median()
	  
	}

	Toolbox.logFine("current average CV median is: %s", mean(unlist(cvList)))

	if(mean(unlist(cvList)) > medianThreshold){
	  
		cvConverge = TRUE
	  
	}

	currentRho <- currentRho * convergenceStep

	}

	#return rho plus root k correction for size of matrix.
	#empirically derived
	return(currentRho * sqrt(crossValidateGroups))

}
