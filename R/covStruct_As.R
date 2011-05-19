## ******************
##   FUNCTION - as 
## ******************

setAs(from="covIso", to="covTensorProduct", 
	def=function(from, to) {
		to <- new("covTensorProduct")
		from.names <- slotNames(from)
		for (i in 1:length(from.names)) {
			slot(to, from.names[i]) <- slot(from, from.names[i])
		}
		to@range.n <- from@d
		to@param.n <- to@range.n 
		to@range.val <- rep(from@range.val, from@d)
		return(to)
}
)


extract.covIso <- function(from) {
		to <- new("covIso")
		from.names <- setdiff(slotNames(from), c("knots", "eta"))
		for (i in 1:length(from.names)) {
			slot(to, from.names[i]) <- slot(from, from.names[i])
		}
		to@range.val <- 1  #rep(1, from@d)
#		to@range.names <- "theta"
		return(to)
}


## fonction ci-dessous : uniquement pour validation (aucun sens sinon...)
as.covIso <- function(object){
	object.iso <- new("covIso")
	object.names <- setdiff(slotNames(object), c("range.n", "shape.n", "shape.names", "shape.val"))
	for (i in 1:length(object.names)) {
		slot(object.iso, object.names[i]) <- slot(object, object.names[i])
	}
	object.iso@param.n <- 1
	object.iso@range.names <- object@range.names[1]
	object.iso@range.val <- object@range.val[1]
	return(object.iso)
}

