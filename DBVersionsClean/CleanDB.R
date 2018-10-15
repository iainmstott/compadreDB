
################################################################################
### CLEAN UP THE COM(P)ADRE DATABASES
### This code calculates variables from the com(p)adre databases to help choose 
### a 'clean' database: models which aren't problematic, which don't violate 
### certain assumptions, which are of a better quality.

#_______________________________________________________________________________
# PACKAGES
# load the necessary packages. Install if needed:
# install.packages(devtools)
# devtools::install_github("iainmstott/popdemo/1.3-0/popdemo")
# devtools::install_github("github.com/jonesor/Rcompadre")
library(popdemo)
library(RCompadre)


#_______________________________________________________________________________
# DATA
# load the data
load("path/to/data/file.Rdata")
# Choose from below depending on which database you're working with:
# DB <- compadre
# DB <- comadre


#_______________________________________________________________________________
# NAs 
# find matrices with NAs and add to metadata
hasNA <- sapply(compadre$mat, function(x){ any(is.na(x$matA)) })
DB$metadata <- data.frame(DB$metadata, hasNA)


#_______________________________________________________________________________
# IRREDUCIBILITY
# find which matrcices are reducible (skipping NA matrices) and add to metadata
Irreducible <- mapply(
    function(na, M){ ifelse(na, FALSE, isIrreducible(M$matA)) },
    na = hasNA, M = DB$mat,
    SIMPLIFY = TRUE
)
DB$metadata <- data.frame(DB$metadata, Irreducible)


#_______________________________________________________________________________
# AVERAGED SURVIVAL
# work out whether there appears to be averaging of survival 

# Survival per stage (column sums of U matrix)
stageSurv <- lapply(lapply(DB$mat, function(M){ M$matU }), colSums)

# find the groups of consecutive average survival using run length encoding
# (rle)
survGroups_rle <- lapply(stageSurv, rle)
# lengths only without corresponding values
survGroupsSummary <- lapply(survGroups_rle, function(G){ 
                                        grpSum <- G$lengths 
                                        names(grpSum) <- NULL
                                        grpSum
                                    }
                                )

# find the different groups the stages belong to
survGroups <- lapply(survGroupsSummary, 
                     function(S){ 
                         unlist(lapply(seq_along(S), function(G){
                             rep(G, S[[G]])
                         }))
                     })

# find the maximum consecutive average in each matrix
maxConsecSurv <- sapply(survGroupsSummary, max)

# extract comadre matrix dimensions
comDim <- DB$metadata$MatrixDimension

# take the lists of consecutive survival and work out a "traffic light"
# system to categorise them. Find the length of the consecutive averages
# for each matrix, then find the maximum for each matrix.
# "No averaging" is considered to be up to 2 consecutive survival of the
# same value.
# G If there is no averaging, (in this case there are NAs in the consecutive
#   survival), then they're rated "G" (green).
# Y If there is averaging, but the number of stages with the same average
#   survival doesn't exceed half the matrix, then the matrix is rated "Y"
#   (yellow).
# R If more than half of the stages have the same average value, then the matrix
#   is rated "R" (red).
survRYG <- character(length(maxConsecSurv))
survRYG[maxConsecSurv %in% c(1,2)] <- "G"
survRYG[maxConsecSurv >= 3 & maxConsecSurv < comDim/2] <- "Y"
survRYG[maxConsecSurv >= comDim / 2] <- "R"

# NAs in the correct places
splitNA <- !(comadre$metadata$MatrixSplit %in% "Divided")
stageSurv[splitNA] <- NA
survGroups_rle[splitNA] <- NA
survGroupsSummary[splitNA] <- NA
survGroups[splitNA] <- NA
maxConsecSurv[splitNA] <- NA
survRYG[splitNA] <- NA

# add to metadata
DB$metadata <- data.frame(DB$metadata, AverageSurv = survRYG)


#_______________________________________________________________________________
# ENHANCED MATRIX COMPOSITE
# Find out whether 'Mean' and 'Pooled' matrices are composites of matrices that
# are in com(p)adre, or matrices that are not found in com(p)adre. From the 
# perspective of pseudoreplication, it's important to differentiate.

objects <- awesomecode(data, etc)
