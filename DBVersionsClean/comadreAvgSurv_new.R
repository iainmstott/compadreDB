################################################################################
# DATA SELECTION STUFF FOR COMADRE MATRICES

#_______________________________________________________________________________
# PACKAGES
library(popdemo)
library(Rcompadre)


#_______________________________________________________________________________
# DATA

# load comadre version 2.0.1
load("C:/Dropbox/Work/Data/comadre/COMADRE_v2.0.1.Rdata")

# Examples to test the code on:
# 461 (Coragyps atratus, the Black Vulture). Has only juvenile averaging
# 646 (Tympanuchus cupido, the Pinnated grouse). Has only adult averaging
# 990 (Capra ibex, the Alpine Ibex). Has all stages averaged
# 1754 (Caretta caretta, the loggerhead turtle). Has many consecutive averages
# 1927 (HIV). Has no averaging
# test <- logical(dim(metadata(comadreS4))[2])
# test[c(461, 646, 990, 1754, 1927)] <- TRUE
# comadreS4 <- subsetDB(comadreS4, test)


#_______________________________________________________________________________
# FINDING AVERAGES OF SURVIVAL ACROSS STAGES

# find out numerically which may be averaged
# Survival (sum per stage)
stageSurv <- lapply(lapply(comadre$mat, function(M){ M$matU }), colSums)

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
comDim <- comadreS4$metadata$MatrixDimension

# take the lists of consecutive survival and work out a "traffic light"
# system to categorise them. Find the length of the consecutive averages
# for each matrix, then find the maximum for each matrix.
# If there is no averaging, (in this case there are NAs in the consecutive
# survival), then they're rated "G" (green).
# If there is averaging, but the number of stages with the same average
# survival doesn't exceed half the matrix, then the matrix is rated "Y"
# (yellow).
# If more than half of the stages have the same average value, then the matrix
# is rated "R" (red).
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

# proportion very naughty averaging
length(which(survRYG %in% "R")) / dim(comadreS4$metadata)[1]
# proportion a little naughty averaging
length(which(survRYG %in% "Y")) / dim(comadreS4$metadata)[1]
# proportion of very good sciencing
length(which(survRYG %in% "G")) / dim(comadreS4$metadata)[1]


#_______________________________________________________________________________
# OTHER PROBLEMS WITH MATRICES

# find matrices with NAs and add to metadata
na <- sapply(comadre$mat, function(x){ any(is.na(x$matA)) })
comadre$metadata <- data.frame(comadre$metadata, na)

# find which matrcices are reducible (skipping NA matrices) and add to metadata
irred <- mapply(
    function(na, com){ ifelse(na, FALSE, isIrreducible(com$matA)) },
    na = na, com = comadre$mat,
    SIMPLIFY = TRUE
)

# find which matrices are greater than 2x2
bigenough <- sapply( comadre$mat,
    function(com){ dim(com$matA)[1] > 2}
)

# We want to use only Individual matrices (not mean or pooled etc.),
# to avoid replication of data. But there are also NAs to avoid of course.
indM <- comadre$metadata$MatrixComposite %in% "Individual"
#...BUT we also want to use matrices that only have one entry per species
#and these could be recorded as many different things (Mean, Pooled, Individual)
numM <- table(comadre$metadata$SpeciesAccepted) #number of matrices per species
lonelySpp <- names(numM)[numM %in% 1] #names of species with only one matrix
single <- comadre$metadata$SpeciesAccepted %in% lonelySpp #corresponding rows of metadata
#checks:
#1 are they really only single species?
length(which(single))==length(lonelySpp)
#2 what are the matrix types for the singles?
unique(comadre$metadata$MatrixComposite[single])
#no seasonals, good


#_______________________________________________________________________________
# SELECTING MATRICES TO USE


#-------------------------------------------------------------------------------
# 1: ONLY GREEN MATRICES
# find which matrices we can use: average survival green, don't have NAs, 
# are irreducible, are larger than 2x2, are individual or single, are
# mammal or bird or reptile, have been divided into U-F-C matrices, have a 
# yearly time step, haven't had a treatment applied, and are from wild 
# populations.
usefulG <- (comadre$metadata$avgSRYG %in% "G" & 
            !na & 
            irred &
            bigenough &
            (indM %in% TRUE | single %in% TRUE) &
            (comadre$metadata$Class %in% "Mammalia" | 
             comadre$metadata$Class %in% "Aves" |
             comadre$metadata$Class %in% "Reptilia") &              #2 NAs!
            comadre$metadata$MatrixSplit %in% "Divided" & 
            comadre$metadata$AnnualPeriodicity %in% "1" &           # LOTS OF NAs!
            comadre$metadata$MatrixTreatment %in% "Unmanipulated" & # LOTS OF NAs!
            comadre$metadata$MatrixCaptivity %in% "W"               # LOTS OF NAs!
)
# How many are useful?
length(which(usefulG))

# subset the database to just the useful matrices
comadreG <- subsetDB(comadre, usefulG)
# How many species are useful?
( usefulSpG <- unique(comadreG$metadata$SpeciesAccepted) )
length(usefulSpG)
# How many families are useful?
( usefulFamG <- unique(comadreG$metadata$Family) )
length(usefulFamG)

# subset the database to just mammals
comadreG_M <- subsetDB(comadreG, comadreG$metadata$Class %in% "Mammalia")
# How many matrices are useful?
dim(comadreG_M$metadata)[1]
# How many species are useful?
( usefulSpG_M <- unique(comadreG_M$metadata$SpeciesAccepted) )
length(usefulSpG_M)
# How many families are useful?
( usefulFamG_M <- unique(comadreG_M$metadata$Family) )
length(usefulFamG_M)

# subset the database to just birds
comadreG_A <- subsetDB(comadreG, comadreG$metadata$Class %in% "Aves")
# How many matrices are useful?
dim(comadreG_A$metadata)[1]
# How many species are useful?
( usefulSpG_A <- unique(comadreG_A$metadata$SpeciesAccepted) )
length(usefulSpG_A)
# How many species are useful?
( usefulFamG_A <- unique(comadreG_A$metadata$Family) )
length(usefulFamG_A)

# subset the database to just reptiles
comadreG_R <- subsetDB(comadreG, comadreG$metadata$Class %in% "Reptilia")
# How many matrices are useful?
dim(comadreG_R$metadata)[1]
# How many species are useful?
( usefulSpG_R <- unique(comadreG_R$metadata$SpeciesAccepted) )
length(usefulSpG_R)
# How many families are useful?
( usefulFamG_R <- unique(comadreG_R$metadata$Family) )
length(usefulFamG_R)


#-------------------------------------------------------------------------------
# 2: GREEN AND YELLOW MATRICES
# find which matrices we can use: average survival green or yellow, don't have 
# NAs, are irreducible, are larger than 2x2, are individual or single, are
# mammal or bird or reptile, have been divided into U-F-C matrices, have a 
# yearly time step, haven't had a treatment applied, and are from wild 
# populations.
usefulGY <- (comadre$metadata$avgSRYG %in% c("G", "Y") & 
             !na & 
             irred &
             bigenough &
             (indM %in% TRUE | single %in% TRUE) &
             (comadre$metadata$Class %in% "Mammalia" | 
              comadre$metadata$Class %in% "Aves" |
              comadre$metadata$Class %in% "Reptilia") &              #2 NAs!
             comadre$metadata$MatrixSplit %in% "Divided" & 
             comadre$metadata$AnnualPeriodicity %in% "1" &           # LOTS OF NAs!
             comadre$metadata$MatrixTreatment %in% "Unmanipulated" & # LOTS OF NAs!
             comadre$metadata$MatrixCaptivity %in% "W"               # LOTS OF NAs!
)
# How many are useful?
length(which(usefulGY))

# subset the database to just the useful matrices
comadreGY <- subsetDB(comadre, usefulGY)
# How many species are useful?
( usefulSpGY <- unique(comadreGY$metadata$SpeciesAccepted) )
length(usefulSpGY)
# How many families are useful?
( usefulFamGY <- unique(comadreGY$metadata$Family) )
length(usefulFamGY)

# subset the database to just mammals
comadreGY_M <- subsetDB(comadreGY, comadreGY$metadata$Class %in% "Mammalia")
# How many matrices are useful?
dim(comadreGY_M$metadata)[1]
# How many species are useful?
( usefulSpGY_M <- unique(comadreGY_M$metadata$SpeciesAccepted) )
length(usefulSpGY_M)
# How many families are useful?
( usefulFamGY_M <- unique(comadreGY_M$metadata$Family) )
length(usefulFamGY_M)

# subset the database to just birds
comadreGY_A <- subsetDB(comadreGY, comadreGY$metadata$Class %in% "Aves")
# How many matrices are useful?
dim(comadreGY_A$metadata)[1]
# How many species are useful?
( usefulSpGY_A <- unique(comadreGY_A$metadata$SpeciesAccepted) )
length(usefulSpGY_A)
# How many species are useful?
( usefulFamGY_A <- unique(comadreGY_A$metadata$Family) )
length(usefulFamGY_A)

# subset the database to just reptiles
comadreGY_R <- subsetDB(comadreGY, comadreGY$metadata$Class %in% "Reptilia")
# How many matrices are useful?
dim(comadreGY_R$metadata)[1]
# How many species are useful?
( usefulSpGY_R <- unique(comadreGY_R$metadata$SpeciesAccepted) )
length(usefulSpGY_R)
# How many families are useful?
( usefulFamGY_R <- unique(comadreGY_R$metadata$Family) )
length(usefulFamGY_R)


#-------------------------------------------------------------------------------
# 3: GREEN, YELLOW AND RED MATRICES
# find which matrices we can use: including all matrices with average survival, 
# don't have NAs, are irreducible, are larger than 2x2, are individual or 
# single, are mammal or bird or reptile, have been divided into U-F-C matrices, 
# have a yearly time step, haven't had a treatment applied, and are from wild 
# populations.
usefulGYR <- (!na & 
              irred &
              bigenough &
              (indM %in% TRUE | single %in% TRUE) &
              (comadre$metadata$Class %in% "Mammalia" | 
               comadre$metadata$Class %in% "Aves" |
               comadre$metadata$Class %in% "Reptilia") &              #2 NAs!
              comadre$metadata$MatrixSplit %in% "Divided" & 
              comadre$metadata$AnnualPeriodicity %in% "1" &           # LOTS OF NAs!
              comadre$metadata$MatrixTreatment %in% "Unmanipulated" & # LOTS OF NAs!
              comadre$metadata$MatrixCaptivity %in% "W"               # LOTS OF NAs!
)
# How many are useful?
length(which(usefulGYR))

# subset the database to just the useful matrices
comadreGYR <- subsetDB(comadre, usefulGYR)
# How many species are useful?
( usefulSpGYR <- unique(comadreGYR$metadata$SpeciesAccepted) )
length(usefulSpGYR)
# How many families are useful?
( usefulFamGYR <- unique(comadreGYR$metadata$Family) )
length(usefulFamGYR)

# subset the database to just mammals
comadreGYR_M <- subsetDB(comadreGYR, comadreGYR$metadata$Class %in% "Mammalia")
# How many matrices are useful?
dim(comadreGYR_M$metadata)[1]
# How many species are useful?
( usefulSpGYR_M <- unique(comadreGYR_M$metadata$SpeciesAccepted) )
length(usefulSpGYR_M)
# How many families are useful?
( usefulFamGYR_M <- unique(comadreGYR_M$metadata$Family) )
length(usefulFamGYR_M)

# subset the database to just birds
comadreGYR_A <- subsetDB(comadreGYR, comadreGYR$metadata$Class %in% "Aves")
# How many matrices are useful?
dim(comadreGYR_A$metadata)[1]
# How many species are useful?
( usefulSpGYR_A <- unique(comadreGYR_A$metadata$SpeciesAccepted) )
length(usefulSpGYR_A)
# How many species are useful?
( usefulFamGYR_A <- unique(comadreGYR_A$metadata$Family) )
length(usefulFamGYR_A)

# subset the database to just reptiles
comadreGYR_R <- subsetDB(comadreGYR, comadreGYR$metadata$Class %in% "Reptilia")
# How many matrices are useful?
dim(comadreGYR_R$metadata)[1]
# How many species are useful?
( usefulSpGYR_R <- unique(comadreGYR_R$metadata$SpeciesAccepted) )
length(usefulSpGYR_R)
# How many families are useful?
( usefulFamGYR_R <- unique(comadreGYR_R$metadata$Family) )
length(usefulFamGYR_R)





################################################################################
# OLD AVERAGE SURVIVAL CODE (WITHOUT USING RLE)

# load packages
library(popdemo)
library(Rcompadre)

# load comadre version 2.0.1
load("C:/Dropbox/Work/Data/comadre/COMADRE_v2.0.1.Rdata")

# Examples to test the code on:
# 461 (Coragyps atratus, the Black Vulture). Has only juvenile averaging
# 646 (Tympanuchus cupido, the Pinnated grouse). Has only adult averaging
# 990 (Capra ibex, the Alpine Ibex). Has all stages averaged
# 1754 (Caretta caretta, the loggerhead turtle). Has many consecutive averages
# 1927 (HIV). Has no averaging
 test <- logical(dim(comadre$metadata)[2])
 test[c(461, 646, 990, 1754, 1927)] <- TRUE
 comadre <- subsetDB(comadre, test)

# find out numerically which may be averaged
# Survival (sum per stage)
stageSurv <- lapply( comadre$mat,
    function(com){ colSums(com$matU) }
)
# Differences between survival values of each stage 
# (length s-1 where s is number of stages)
survDiff <- lapply(stageSurv, diff)
# Logical: are differences equal to zero? 
# (i.e. survival of both stages is the same)
# s-1 length
survZeroDiff <- lapply(survDiff, function(sDiff){ !(as.logical(sDiff)) })
# Index of which survival values are the same
# e.g. 1 2 3 would mean differences between 1 & 2, 2 & 3, and 3 & 4 are the same
survWhichZeroDiff <- lapply(survZeroDiff, which)
# NAs in the correct places
splitNA <- !(comadre$metadata$MatrixSplit %in% "Divided")
stageSurv[splitNA] <- NA
survDiff[splitNA] <- NA
survZeroDiff[splitNA] <- NA
survWhichZeroDiff[splitNA] <- NA

# locate 'breaks' in consecutive survival (i.e. stages where survival stops
# being the same). This helps identify different consecutive averages in the 
# same model.
survBreaks <- lapply(survWhichZeroDiff, function(nzDiff){
        c(0, which(diff(nzDiff) != 1), length(nzDiff))}
)

# Now turn this into lists of lists of stages which have consecutive average 
# survival.
consecSurv <- mapply(function(sBreaks, wzDiff){
        sapply(seq(length(sBreaks) - 1), function(i){
            stages1 <- wzDiff[(sBreaks[i] + 1):sBreaks[i + 1]]
            stages2 <- c(stages1, stages1[length(stages1)] + 1)
            list(stages2)
        })
    },
    sBreaks = survBreaks, wzDiff = survWhichZeroDiff,
    SIMPLIFY = TRUE
)

# take the lists of consecutive survival and work out a "traffic light"
# system to categorise them. Find the length of the consecutive averages
# for each matrix, then find the maximum for each matrix.
# If there is no averaging, (in this case there are NAs in the consecutive
# survival), then they're rated "G" (green).
# If there is averaging, but the number of stages with the same average
# survival doesn't exceed half the matrix, then the matrix is rated "Y"
# (yellow).
# If more than half of the stages have the same average value, then the matrix
# is rated "R" (red).
survRYG <- mapply( function(comDim, conS){
        maxConsecSurv <- max(sapply(conS, length))
        if(any(is.na(unlist(conS)))) maxConsecSurv <- 1
        if(maxConsecSurv %in% c(1,2)) RYG <- "G"
        if(maxConsecSurv >= 3 & maxConsecSurv < comDim/2) RYG <- "Y"
        if(maxConsecSurv >= comDim/2) RYG <- "R"
        RYG
    },
    comDim = comadre$metadata$MatrixDimension, conS = consecSurv,
    SIMPLIFY = TRUE
)
# note... might want to look at the survival of juveniles and adults 
# separately.

# proportion very naughty averaging
length(which(survRYG %in% "R") / dim(comadre$metadata)[1]
# proportion a little naughty averaging
length(which(survRYG %in% "Y") / dim(comadre$metadata)[1]
# proportion of very good sciencing
length(which(survRYG %in% "G") / dim(comadre$metadata)[1]


