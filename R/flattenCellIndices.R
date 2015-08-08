#------------------------------
# FLATTENING OF CELLS
#------------------------------

flattenCellIndices <- function(cells) {

	## Flatten cell data and clean the names
	cells <- unlist(cells)						
	names(cells) <- gsub("[.](groups|indices.*)", "", names(cells))
	unitNames <- gsub("[.].*", "", names(cells))				
	groupNames <- gsub(".*[.]", "", names(cells))
	groupNames <- paste(unitNames, groupNames, sep=".")

	return(data.frame(unitNames=unitNames, groupNames=groupNames, cell=cells))	

}

