<<<<<<< .mine
#---------------------
# TUKEY BIWEIGHT
#---------------------
# Robust average of vector x with missing values ignored in the calculation.
tukey <- function(x, c=5, epsilon=1e-4, na.rm=TRUE) {
	x <- x[!is.na(x)]
	return <- tukey.biweight(x, c=c, epsilon=epsilon)
}
=======
#---------------------
# TUKEY BIWEIGHT
#---------------------
# Robust average of vector x with missing values ignored in the calculation.
tukey <- function(x,c=5,epsilon=1e-4,na.rm=TRUE) {
	x <- x[!is.na(x)]
	return <- tukey.biweight(x,c=c,epsilon=epsilon)
}>>>>>>> .r101298
