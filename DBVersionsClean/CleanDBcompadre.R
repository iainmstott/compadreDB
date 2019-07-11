
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
cdb_DB <- as_cdb(compadre)


#_______________________________________________________________________________
# NAs, zeroes, sums, singularity, ergodicity
# run cdb_flag to check whether there are common problems with the matrices
DB <- cdb_flag(cdb_DB)



#_______________________________________________________________________________
# AVERAGED SURVIVAL
# work out whether there appears to be averaging of survival 


# We want to work with U matrices that don't contain NAs and do contain nonzero entries:
# we can calculate stuff from these. DB$check_nonzero_U is a variable that is TRUE
# if matU contains at least one nonzero number and no NAs, FALSE otherwise.
DB$check_nonzero_U <- !DB$check_NA_U & !DB$check_zero_U

# note that these variables are connected to the "MatrixSplit" variable, 
# as a U matrix that indivisible contains all NAs. However a matrix can 
# be 'Divided' and contain either NAs, or no nonzero entries:
which(DB$check_NA_U + (DB$MatrixSplit == "Divided") == 2)
which(DB$check_zero_U + (DB$MatrixSplit == "Divided") == 2)
# a matrix can also be recorded as indivisible but have a nonzero U matrix:
which(DB$check_nonzero_U + (DB$MatrixSplit == "Indivisible") == 2)

# Survival per stage (column sums of U matrix). stageSurv is a list containing 
# the column sums of matU (survival per stage)
stageSurv <- lapply(matU(DB), colSums)
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
                         sG = survGroups, sS = stageSurv,
                         mUz = DB$check_zero_U)

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
nStages <- DB$MatrixDimension

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
survRYGNoZero[maxConsecSurvNoZero >= 2 &
              maxConsecSurvNoZero <= nStages / 2] <- "Y"
survRYGNoZero[maxConsecSurvNoZero >= 2 &
              maxConsecSurvNoZero > nStages / 2] <- "R"

# make very sure NAs in the correct places
stageSurv[!DB$check_nonzero_U] <- NA
survGroups_rle[!DB$check_nonzero_U] <- NA
survGroupsSummary[!DB$check_nonzero_U] <- NA
survGroups[!DB$check_nonzero_U] <- NA
survGroupsZero[!DB$check_nonzero_U] <- NA
survGroupsNoZero[!DB$check_nonzero_U] <- NA
survGroupsNoZero_rle[!DB$check_nonzero_U] <- NA
survGroupsNoZeroSummary[!DB$check_nonzero_U] <- NA
maxConsecSurv[!DB$check_nonzero_U] <- NA
maxConsecSurvNoZero[!DB$check_nonzero_U] <- NA
survRYG[!DB$check_nonzero_U] <- NA
survRYGNoZero[!DB$check_nonzero_U] <- NA

# add to metadata
DB$UcolSum <- stageSurv
DB$AvSurvGrp <- survGroupsZero
DB$AvSurvMax <- maxConsecSurv
DB$AvSurvMaxNoZero <- maxConsecSurvNoZero
DB$AvSurvRYG <- survRYG
DB$AvSurvRYGNoZero <- survRYGNoZero


#_______________________________________________________________________________
# AVERAGED FECUNDITY
# work out whether there appears to be averaging of fecundity 

# We want to work with F matrices that don't contain NAs and do contain nonzero entries:
# we can calculate stuff from these. DB$check_nonzero_F is a variable that is TRUE
# if matU contains at least one nonzero number and no NAs, FALSE otherwise.
DB$check_zero_F <- sapply(matF(DB), function(M){ all(M %in% 0) })

DB$check_nonzero_F <- !DB$check_NA_F & !DB$check_zero_F

# note that these variables are connected to the "MatrixSplit" variable, 
# as an F matrix that indivisible contains all NAs. However a matrix can 
# be 'Divided' and contain either NAs, or no nonzero entries:
which(DB$check_NA_F + (DB$MatrixSplit == "Divided") == 2)
which(DB$check_zero_F + (DB$MatrixSplit == "Divided") == 2)
# a matrix can also be recorded as indivisible but have a nonzero F matrix:
which(DB$check_nonzero_F + (DB$MatrixSplit == "Indivisible") == 2)

# Fecundity per stage (column sums of F matrix). stageFec is a list containing 
# the column sums of matF (the total sexual reproduction per stage)
stageFec <- lapply(matF(DB), colSums)
# First fecund stage. firstFec is a variable with the number of the first stage
# where stageFec is greater than zero.
firstFec <- mapply( function(sF, Fna, Fz, Fnz){ 
                        if(Fna | Fz) return(NA)
                        if(Fnz){
                            return(min(which(!(sF %in% 0))))
                        }
                    }, 
                    sF = stageFec, Fna = DB$check_NA_F,
                    Fz = DB$check_zero_F, Fnz = DB$check_nonzero_F)
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
                        fG = fecGroups, sFM = stageFecMature,
                        mFz = DB$check_zero_F)

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
nFstages <- DB$MatrixDimension - firstFec + 1 
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
stageFec[!DB$check_nonzero_F] <- NA
firstFec[!DB$check_nonzero_F] <- NA
stageFecMature[!DB$check_nonzero_F] <- NA
fecGroups_rle[!DB$check_nonzero_F] <- NA
fecGroupsSummary[!DB$check_nonzero_F] <- NA
fecGroups[!DB$check_nonzero_F] <- NA
fecGroupsZero[!DB$check_nonzero_F] <- NA
fecGroupsNoZero[!DB$check_nonzero_F] <- NA
fecGroupsNoZero_rle[!DB$check_nonzero_F] <- NA
fecGroupsNoZeroSummary[!DB$check_nonzero_F] <- NA
maxConsecFec[!DB$check_nonzero_F] <- NA
maxConsecFecNoZero[!DB$check_nonzero_F] <- NA
nFstages[!DB$check_nonzero_F] <- NA
fecRYG[!DB$check_nonzero_F] <- NA
fecRYGNoZero[!DB$check_nonzero_F] <- NA

# add to metadata
DB$FcolSum <- stageFec
DB$FirstF <- firstFec
DB$FnStages <- nFstages
# DB$RnStages <- nRstages
DB$AvFecGrp <- fecGroupsZero
DB$AvFecMaxNoZero <- maxConsecFecNoZero
DB$AvFecRYG <- fecRYG
DB$AvFecRYGNoZero <- fecRYGNoZero


#_______________________________________________________________________________
# NEW DATABASE

compadre <- DB

# remove surplus objects
rm(list = ls()[!ls() %in% "compadre"])

# load the phylogeny
phy <- read.tree("Phylogeny/COMPADRE-COMADRE_Phylo_June_16_2019.tre")

# save the .Rdata object
save.image("DBVersionsClean/compadre_v.5.0.0_Clean.Rdata")
