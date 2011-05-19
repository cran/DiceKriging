drop.response <- function(formula, data) {

	#formula <- as.formula(formula)
	tt <- terms(formula, data=data)
	if (length(attr(tt, "term.labels"))!=0) {
		tt0 <- delete.response(tt)
		formula <- reformulate(attr(tt0, "term.labels"))
	} else formula <- ~1
	return(formula)
}
