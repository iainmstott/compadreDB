
################################################################################
### CLEAN UP THE COM(P)ADRE DATABASES
### This code calculates variables from the com(p)adre databases to help choose 
### a 'clean' database: models which aren't problematic, which don't violate 
### certain assumptions, which are of a better quality.

#_______________________________________________________________________________
# PACKAGES
# load the necessary packages. Install if needed:
### NOTE!!! Replace "VERSION" in popdemo path with desired popdemo version
# install.packages(devtools)
# devtools::install_github("iainmstott/popdemo/VERSION/popdemo")
# devtools::install_github("github.com/jonesor/Rcompadre")

# load packages
library(popdemo)
library(Rcompadre)
library(ape)


#_______________________________________________________________________________
# DATA --- EDIT AS NEEDED
# load the data
#load("path/to/data/file.Rdata")
load("DBVersions/COMPADRE_v.5.0.0.Rdata")
# Choose from below depending on which database you're working with:
DB <- compadre


#_______________________________________________________________________________
# NAs 
# find if matrices have NAs and add to metadata
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



# First we need to work out which U matrices we can get survival 
# estimates from.
# Some U matrices contain NAs: we can't calculate stuff from these. matUna 
# is a variable that is TRUE if a U matrix contains at least one NA, FALSE
# otherwise.
matUna <- sapply(DB$mat, function(M){ any(is.na(M$matU)) })
# Some U matrices don't contain any nonzero entries: we can't calculate 
# stuff from these. matUzero is a variable that is TRUE if all entries in 
# matU are zero, FALSE otherwise.
matUzero <- sapply(DB$mat, function(M){ all(M$matU %in% 0) })
# The rest of the U matrices don't contain NAs and do contain nonzero entries:
# we can calculate stuff from these. matUnonzero is a variable that is TRUE
# if matU contains at least one nonzero number and no NAs, FALSE otherwise.
matUnonzero <- sapply(DB$mat, function(M){ 
    !any(is.na(M$matU)) & any(M$matU > 0) 
})
# note that these variables are connected to the "MatrixSplit" variable, 
# as a U matrix that indivisible contains all NAs. However a matrix can 
# be 'Divided' and contain either NAs, or no nonzero entries:
which(matUna + (DB$metadata$MatrixSplit == "Divided") == 2)
which(matUzero + (DB$metadata$MatrixSplit == "Divided") == 2)
# a matrix can also be recorded as indivisible but have a nonzero U matrix:
which(matUnonzero + (DB$metadata$MatrixSplit == "Indivisible") == 2)

# Survival per stage (column sums of U matrix). stageSurv is a list containing 
# the column sums of matU (survival per stage)
stageSurv <- lapply(lapply(DB$mat, function(M){ M$matU }), colSums)
# find the groups of consecutive average survival using run length encoding
# (rle)
survGroups_rle <- lapply(stageSurv, 
                        function(sS){ if(!is.null(sS)) rle(sS) } )
# lengths only, without corresponding values
survGroupsSummary <- lapply(survGroups_rle, function(G){ 
                                        grpSum <- G$lengths 
                                        names(grpSum) <- NULL
                                        grpSum
                                    }
                                )
# find the different 'groups' of averages the stages belong to. survGroups 
# assigns groups to stageSurv. For example, a matrix with stageSurv
# of c(0.1, 0.12, 0.5, 0.5, 0.5, 0.9) would have fecGroups of
# c(1, 2, 3, 3, 3, 4)
survGroups <- lapply(survGroupsSummary, 
                     function(S){ 
                         unlist(lapply(seq_along(S), function(G){
                             rep(G, S[[G]])
                         }))
                     })

# identify stages with zero fecundity (rather than being "averages" of 
# survival, these have merely not been observed), and flag in the survGroups
# variable. For example, a matrix with stageSurv
# of c(0.1, 0.12, 0.5, 0.5, 0.5, 0, 0, 0.9) would have survGroupsZero of
# c(1, 2, 3, 3, 3, 0, 0, 5)
survGroupsZero <- mapply(function(sG, sS, mUz){ 
                             if(is.null(sG) | mUz) sGZ <- NULL
                             if(!(is.null(sG)) & !mUz){ 
                                 sGZ <- sG
                                 sGZ[sS %in% 0] <- 0
                             }
                             sGZ
                         },
                         sG = survGroups, sS = stageSurv, mUz = matUzero)

#collapse these so that zero survival groups aren't included. For example,
# a matrix with stageSurv of c(0.1, 0.12, 0.5, 0.5, 0.5, 0, 0, 0.9) 
# would have survGroupsNoZero of c(1, 2, 3, 3, 3, 5)
survGroupsNoZero <- lapply(survGroupsZero, function(sGZ) sGZ[!(sGZ %in% 0)] )
# run length encoding on new no-zero groups
survGroupsNoZero_rle <- lapply(survGroupsNoZero, 
                            function(sGNZ){ if(!is.null(sGNZ)) rle(sGNZ) } )
# lengths only, without corresponding values
survGroupsNoZeroSummary <- lapply(survGroupsNoZero_rle, function(G){ 
                                    grpSum <- G$lengths 
                                    names(grpSum) <- NULL
                                    grpSum
                                }
                            )

# find the maximum number of consecutive averaged stages in each matrix
# including zeroes:
maxConsecSurv <- sapply(survGroupsSummary, 
                              function(sGS){
                                      if(!is.null(sGS)) mCS <- max(sGS)
                                      if(is.null(sGS)) mCS <- NA
                                      mCS
                              })
# excluding zeroes:
maxConsecSurvNoZero <- sapply(survGroupsNoZeroSummary, 
                              function(sGNZS){
                                      if(!is.null(sGNZS)) mCSNZ <- max(sGNZS)
                                      if(is.null(sGNZS)) mCSNZ <- NA
                                      mCSNZ
                              })

# extract number of stages
nStages <- DB$metadata$MatrixDimension

# take the lists of consecutive survival and work out a "traffic light"
# system to categorise them. Find the length of the consecutive averages
# for each matrix, then find the maximum for each matrix.
# G If a all stages have different survival, the matrix is rated "G" (green).
# Y If there is apparent averaging (2 or more consecutive survival values the 
#   same), but the number of stages with the same value of survival doesn't 
#   exceed half the number of stages, then the matrix is rated "Y" (yellow).
# R If half or more of the fecund stages have the same average value, then 
#   the matrix is rated "R" (red), unless there is only one stage.
# including Zeroes:
survRYG <- character(length(maxConsecSurv))
survRYG[maxConsecSurv %in% 1] <- "G"
survRYG[maxConsecSurv >= 2 & maxConsecSurv <= nStages / 2] <- "Y"
survRYG[maxConsecSurv >= 2 & maxConsecSurv > nStages / 2] <- "R"
# excluding zeroes:
survRYGNoZero <- character(length(maxConsecSurvNoZero))
survRYGNoZero[maxConsecSurvNoZero %in% 1] <- "G"
survRYGNoZero[maxConsecSurvNoZero >= 2 & maxConsecSurvNoZero <= nStages / 2] <- "Y"
survRYGNoZero[maxConsecSurvNoZero >= 2 & maxConsecSurvNoZero > nStages / 2] <- "R"

# make very sure NAs in the correct places
stageSurv[!matUnonzero] <- NA
survGroups_rle[!matUnonzero] <- NA
survGroupsSummary[!matUnonzero] <- NA
survGroups[!matUnonzero] <- NA
survGroupsZero[!matUnonzero] <- NA
survGroupsNoZero[!matUnonzero] <- NA
survGroupsNoZero_rle[!matUnonzero] <- NA
survGroupsNoZeroSummary[!matUnonzero] <- NA
maxConsecSurv[!matUnonzero] <- NA
maxConsecSurvNoZero[!matUnonzero] <- NA
survRYG[!matUnonzero] <- NA
survRYGNoZero[!matUnonzero] <- NA

# add to metadata
DB$metadata$UcolSum <- stageSurv
DB$metadata$AvSurvGrp <- survGroupsZero
DB$metadata$AvSurvMax <- maxConsecSurv
DB$metadata$AvSurvMaxNoZero <- maxConsecSurvNoZero
DB$metadata$AvSurvRYG <- survRYG
DB$metadata$AvSurvRYGNoZero <- survRYGNoZero



#_______________________________________________________________________________
# AVERAGED FECUNDITY
# work out whether there appears to be averaging of fecundity 

# First we need to work out which F matrices we can get fecundity 
# estimates from.
# Some F matrices contain NAs: we can't calculate stuff from these. matFna 
# is a variable that is TRUE if an F matrix contains at least one NA, FALSE
# otherwise.
matFna <- sapply(DB$mat, function(M){ any(is.na(M$matF)) })
# Some F matrices don't contain any nonzero entries: we can't calculate 
# stuff from these. matFzero is a variable that is TRUE if all entries in 
# matF are zero, FALSE otherwise.
matFzero <- sapply(DB$mat, function(M){ all(M$matF %in% 0) })
# The rest of the F matrices don't contain NAs and do contain nonzero entries:
# we can calculate stuff from these. matFnonzero is a variable that is TRUE
# if matF contains at least one nonzero number and no NAs, FALSE otherwise.
matFnonzero <- sapply(DB$mat, function(M){ 
    !any(is.na(M$matF)) & any(M$matF > 0) 
})
# note that these variables are connected to the "MatrixSplit" variable, 
# as a matrix that indivisible contains all NAs. However a matrix can 
# be 'Divided' and contain either NAs, or no nonzero entries:
which(matFna + (DB$metadata$MatrixSplit == "Divided") == 2)
which(matFzero + (DB$metadata$MatrixSplit == "Divided") == 2)
# a matrix can also be recorded as indivisible but have a nonzero F matrix:
which(matFnonzero + (DB$metadata$MatrixSplit == "Indivisible") == 2)

# Fecundity per stage (column sums of F matrix). stageFec is a list containing 
# the column sums of matF (the total sexual reproduction per stage)
stageFec <- lapply(lapply(DB$mat, function(M){ M$matF }), colSums)
# First fecund stage. firstFec is a variable with the number of the first stage
# where stageFec is greater than zero.
firstFec <- mapply( function(sF, Fna, Fz, Fnz){ 
                        if(Fna | Fz) return(NA)
                        if(Fnz){
                            return(min(which(!(sF %in% 0))))
                        }
                    }, 
            sF = stageFec, Fna = matFna, Fz = matFzero, Fnz = matFnonzero)
# extract only fecund stages (those to RHS of first fecund stage). stageFecMature
# is a list containing stageFec for only those stages after and including 
# firstFec.
stageFecMature <- mapply(function(sF, fF){ 
                             if(!is.na(fF)) sF[fF:length(sF)]
                         }, 
                         sF = stageFec, fF = firstFec)
# find the groups of consecutive average fecundity using run length encoding
# (rle)
fecGroups_rle <- lapply(stageFecMature, 
                        function(sF){ if(!is.null(sF)) rle(sF) } )
# lengths only, without corresponding values
fecGroupsSummary <- lapply(fecGroups_rle, function(G){ 
                                        grpSum <- G$lengths 
                                        names(grpSum) <- NULL
                                        grpSum
                                    }
                                )
# find the different 'groups' of averages the stages belong to. fecGroups 
# assigns groups to stageFecMature. For example, a matrix with stageFecMature
# of c(0.1, 0.12, 0.5, 0.5, 0.5, 0, 0, 0.9) would have fecGroups of
# c(1, 2, 3, 3, 3, 4, 4, 5)
fecGroups <- lapply(fecGroupsSummary, 
                     function(S){ 
                         unlist(lapply(seq_along(S), function(G){
                             rep(G, S[[G]])
                         }))
                     })

# identify stages with zero fecundity (rather than being "averages" of 
# fecundity, these have merely not been observed), and flag in the fecGroups
# variable. For example, a matrix with stageFecMature
# of c(0.1, 0.12, 0.5, 0.5, 0.5, 0, 0, 0.9) would have fecGroupsZero of
# c(1, 2, 3, 3, 3, 0, 0, 5)
fecGroupsZero <- mapply(function(fG, sFM, mFz){ 
                            if(is.null(fG) | mFz) fGZ <- NULL
                            if(!(is.null(fG)) & !mFz){ 
                                fGZ <- fG
                                fGZ[sFM %in% 0] <- 0
                            }
                            fGZ
                        },
                        fG = fecGroups, sFM = stageFecMature, mFz = matFzero)

#collapse these so that zero fecundity groups aren't included. For example,
# a matrix with stageFecMature of c(0.1, 0.12, 0.5, 0.5, 0.5, 0, 0, 0.9) 
# would have fecGroupsNoZero of c(1, 2, 3, 3, 3, 5)
fecGroupsNoZero <- lapply(fecGroupsZero, function(fGZ) fGZ[!(fGZ %in% 0)] )
# run length encoding on new no-zero groups
fecGroupsNoZero_rle <- lapply(fecGroupsNoZero, 
                            function(fGNZ){ if(!is.null(fGNZ)) rle(fGNZ) } )
# lengths only, without corresponding values
fecGroupsNoZeroSummary <- lapply(fecGroupsNoZero_rle, function(G){ 
                                    grpSum <- G$lengths 
                                    names(grpSum) <- NULL
                                    grpSum
                                }
                            )

# find the maximum number of consecutive averaged stages in each matrix
# including zeroes:
maxConsecFec <- sapply(fecGroupsSummary, 
                              function(fGS){
                                      if(!is.null(fGS)) mCF <- max(fGS)
                                      if(is.null(fGS)) mCF <- NA
                                      mCF
                              })
# excluding zeroes:
maxConsecFecNoZero <- sapply(fecGroupsNoZeroSummary, 
                             function(fGNZS){
                                 if(!is.null(fGNZS)) mCFNZ <- max(fGNZS)
                                 if(is.null(fGNZS)) mCFNZ <- NA
                                 mCFNZ
                            })

# extract number of fertile stages
nFstages <- DB$metadata$MatrixDimension - firstFec + 1 
# extract number of reproducing stages
# nRstages <- sapply(stageFec, function(sF) length(which(sF != 0)))

# take the lists of consecutive fecundity and work out a "traffic light"
# system to categorise them. Find the length of the consecutive averages
# for each matrix, then find the maximum for each matrix.
# G If a all stages have different fertility, the matrix is rated "G" (green).
# Y If there is apparent averaging (2 or more consecutive fertility values the 
#   same), but the number of stages with the same value of fertility doesn't 
#   exceed half the number of fecund stages, then the matrix is rated "Y" (yellow).
# R If half or more of the fecund stages have the same average value, then 
#   the matrix is rated "R" (red), unless there is only one fecund stage.
# (NOTE: these could be changed to use nRstages instead of nFstages)
# including Zeroes:
fecRYG <- character(length(maxConsecFec))
fecRYG[maxConsecFec %in% 1] <- "G"
fecRYG[maxConsecFec >= 2 & maxConsecFec <= nFstages / 2] <- "Y"
fecRYG[maxConsecFec >= 2 & maxConsecFec > nFstages / 2] <- "R"
# excluding zeroes:
fecRYGNoZero <- character(length(maxConsecFecNoZero))
fecRYGNoZero[maxConsecFecNoZero %in% 1] <- "G"
fecRYGNoZero[maxConsecFecNoZero >= 2 & maxConsecFecNoZero <= nFstages / 2] <- "Y"
fecRYGNoZero[maxConsecFecNoZero >= 2 & maxConsecFecNoZero > nFstages / 2] <- "R"

# make very sure NAs in the correct places
stageFec[!matFnonzero] <- NA
firstFec[!matFnonzero] <- NA
stageFecMature[!matFnonzero] <- NA
fecGroups_rle[!matFnonzero] <- NA
fecGroupsSummary[!matFnonzero] <- NA
fecGroups[!matFnonzero] <- NA
fecGroupsZero[!matFnonzero] <- NA
fecGroupsNoZero[!matFnonzero] <- NA
fecGroupsNoZero_rle[!matFnonzero] <- NA
fecGroupsNoZeroSummary[!matFnonzero] <- NA
maxConsecFec[!matFnonzero] <- NA
maxConsecFecNoZero[!matFnonzero] <- NA
nFstages[!matFnonzero] <- NA
fecRYG[!matFnonzero] <- NA
fecRYGNoZero[!matFnonzero] <- NA

# add to metadata
DB$metadata$FcolSum <- stageFec
DB$metadata$FirstF <- firstFec
DB$metadata$FnStages <- nFstages
# DB$metadata$RnStages <- nRstages
DB$metadata$AvFecGrp <- fecGroupsZero
DB$metadata$AvFecMaxNoZero <- maxConsecFecNoZero
DB$metadata$AvFecRYG <- fecRYG
DB$metadata$AvFecRYGNoZero <- fecRYGNoZero


#_______________________________________________________________________________
# NEW DATABASE

compadre <- DB

# remove surplus objects
rm(list = ls()[!ls() %in% "compadre"])

# load the phylogeny
phy <- read.tree("Phylogeny/COMPADRE-COMADRE_Phylo_June_16_2019.tre")

# save the .Rdata object
save.image("DBVersionsClean/compadre_v.5.0.0_Clean.Rdata")
