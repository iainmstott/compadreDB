library(dplyr)
library(stringr)
library(Rcompadre)




### COMPADRE ###################################################################

# load development, current and legacy versions
load("DBVersionsX/COMPADRE_v.X.X.X.Rdata")
cpxxx <- as_cdb(compadre)
rm(compadre)
load("DBVersions/COMPADRE_v.5.0.0.Rdata")
cp500 <- as_cdb(compadre)
rm(compadre)
load("DBVersions/COMPADRE_v.4.0.1.Rdata")
cp401 <- as_cdb(compadre)
rm(compadre)

( cpxxx_vers <- table(cpxxx$VersionRelease, useNA = "always") )
names(cpxxx_vers)[is.na(names(cpxxx_vers))] <- "NA"
sum(cpxxx_vers) - sum(cpxxx_vers[c("X", "NA")]) # 8628
nrow(cp500) # 9121
# website: 8616



# match up the SpeciesAuthor from 2.0.1 to 3.0.0
# create a numeric (2.0.1 dim) which will be the indices used to match
# create a vector of TRUE (2.0.1 dim) which will be used to flag problem matrices
# for each SpeciesAuthor... (obtained using unique(SpeciesAuthor))
## check that all the matching numbers are consecutive (diff == 1)
### if they aren't consecutive, skip indices & flag all
### if they are, create consecutive index from smallest to largest number
###   returned, put in the right place and ignore flag

# dim for both versions
dim401 <- dim(cp401)[1]
dim500 <- dim(cp500)[1]

# SpeciesAuthor for both versions
SA401 <- cp401$SpeciesAuthor
SA500 <- cp500$SpeciesAuthor
# unique SA for 2.0.1
unique_SA401 <- unique(SA401)
num_SA401 <- length(unique_SA401)

# list of matrices for both versions
matA401 <- matA(cp401)
matA500 <- matA(cp500)

### set up for loop

## lists
# list of integers to place each SA in 401
SA401_m_SA401 <- list()
# list of integers to match SA from 401 to 500
SA500_m_SA401 <- list()
# list of characters to hold each matrix string (401) by SA
matstring401 <- list()
# list of characters to hold each matrix string (500) by SA
matstring500 <- list()

## species data (length of unique SA values)
# dataframe of species info for SA in 401: number of matrices for that SA in 
# 401 and 500, whether matrices for that SA in 401 and 500 are consecutive rows,
# whether number of matrices in for that SA 401 and 500 are equal to one another
cp401_SAdata <- data.frame(unique_SA401,
                         num401 = numeric(num_SA401), num500 = numeric(num_SA401),
                         numequal = logical(num_SA401),
                         consec401 = logical(num_SA401), consec500 = logical(num_SA401),
                         check1 = logical(num_SA401)
)

## matrix data (length of num matrices in 2.0.1 and 3.0.0 respectively)
# dataframe of matrix infor for matrices in 401: 
# splitpos indicate whether matrices for that SA are consecutive in both 
# 401 and 500, 
# numequal indicates whether the number of matrices for that SA 
# is equal in 401 and 500, 
# match_SA500 is the index matching the that matrix to its supposed 
# corresponding row for that matrix in 500 (only for any matrix where check1 
# is false),
# check indicates whether the matrix needs checking in the 2nd round 
# (numequal is false, splitpos is true, &/or 
# character string version of matrix doesn't match in 401 and 500).
SA_matdata401 <- data.frame(SA401, SA_in500 = logical(dim401),
                            numequal = logical(dim401),
                            splitpos = logical(dim401),
                            match_SA500 = numeric(dim401),
                            check = logical(dim401)
)
# dataframe of matrix info for matrices in 500: 
# same as above, except match_SA401 reverse-matches the matrix to its supposed 
# corresponding row in 401.
SA_matdata500 <- data.frame(SA500, SA_in401 = logical(dim500),
                            numequal = logical(dim500),
                            splitpos = logical(dim500),
                            match_SA401 = numeric(dim500),
                            check = logical(dim500)
)

### loop to fill in all the above empty information 
### (sp = one SpeciesAuthor from 401)
for (sp in seq_along(unique_SA401)) {

    ## lists
    # positions of sp in 401
    SA401_m_SA401[[sp]] <- which(SA401 %in% unique_SA401[sp])
    names(SA401_m_SA401)[sp] <- unique_SA401[sp]
    # positions of sp in 500
    SA500_m_SA401[[sp]] <- which(SA500 %in% unique_SA401[sp])
    names(SA500_m_SA401)[sp] <- unique_SA401[sp]
    matstring401[[sp]] <- sapply(matA401[SA401_m_SA401[[sp]]], mat_to_string)
    names(matstring401)[sp] <- unique_SA401[sp]
    matstring500[[sp]] <- sapply(matA500[SA500_m_SA401[[sp]]], mat_to_string)
    names(matstring500)[sp] <- unique_SA401[sp]

    ## Species data frame
    # how many matrices for sp in each of the two versions?
    cp401_SAdata[sp, c("num401", "num500")] <-
        c(length(SA401_m_SA401[[sp]]), length(SA500_m_SA401[[sp]]))
    # are the number of matrices for sp in each data base equal?
    cp401_SAdata$numequal[sp] <- cp401_SAdata[sp, "num401"] == cp401_SAdata[sp, "num500"]
    # are all of the index values consecutive?
    # (If yes hopefully the matrices can be matched over in chunks)
    # Note: species with only 1 matrix are automatically consecutive
    #       species with 0 matrix are NA
    if (cp401_SAdata[sp, "num401"] > 1) {
        cp401_SAdata[sp, "consec401"] <- all(diff(SA401_m_SA401[[sp]]) == 1)
    }
    if (cp401_SAdata[sp, "num500"] > 1) {
        cp401_SAdata[sp, "consec500"] <- all(diff(SA500_m_SA401[[sp]]) == 1)
    }
    if (cp401_SAdata[sp, "num401"] == 1) cp401_SAdata[sp, "consec401"] <- TRUE
    if (cp401_SAdata[sp, "num500"] == 1) cp401_SAdata[sp, "consec500"] <- TRUE
    if (cp401_SAdata[sp, "num401"] == 0) cp401_SAdata[sp, "consec401"] <- NA
    if (cp401_SAdata[sp, "num500"] == 0) cp401_SAdata[sp, "consec500"] <- NA
    # First check: are there equal numbers of matrices for sp in both 401 and
    # 500, with all matrices on consecutive rows in both data bases?
    cp401_SAdata[sp, "check1"] <- !all(c(cp401_SAdata$numequal[sp],
                                       cp401_SAdata$consec401[sp],
                                       cp401_SAdata$consec500[sp]))

    ## matrix data frames
    # is the species in the other data base?
    SA_matdata401[SA401_m_SA401[[sp]], "SA_in500"] <- cp401_SAdata[sp, "num500"] > 0
    SA_matdata500[SA500_m_SA401[[sp]], "SA_in401"] <- cp401_SAdata[sp, "num401"] > 0
    # are the number of matrices for each sp equal in 401 and 500?
    SA_matdata401[SA401_m_SA401[[sp]], "numequal"] <- cp401_SAdata$numequal[sp]
    SA_matdata401[SA500_m_SA401[[sp]], "numequal"] <- cp401_SAdata$numequal[sp]
    # are the matrices on consecutive rows in 401 / 500?
    SA_matdata401[SA401_m_SA401[[sp]], "splitpos"] <-
        any(cp401_SAdata$consec401[sp], cp401_SAdata$consec500[sp])
    SA_matdata500[SA500_m_SA401[[sp]], "splitpos"] <-
        any(cp401_SAdata$consec401[sp], cp401_SAdata$consec500[sp])
    # 1. match the rows for sp in 500 to the rows in 401 and vice versa
    #    (NA if check1 is TRUE as there is no easy match)
    # 2. check: TRUE if check1 is TRUE, also TRUE if check1 is FALSE but 
    #    the matrices, represented as strings, are not the same.
    if (cp401_SAdata$check1[sp]) {
        SA_matdata401[SA401_m_SA401[[sp]], "match_SA500"] <- 0
        SA_matdata500[SA500_m_SA401[[sp]], "match_SA401"] <- 0
        SA_matdata401[SA401_m_SA401[[sp]], "check"] <- TRUE
        SA_matdata500[SA500_m_SA401[[sp]], "check"] <- TRUE
    }
    if (!cp401_SAdata$check1[sp]) {
        SA_matdata401[SA401_m_SA401[[sp]], "match_SA500"] <- SA500_m_SA401[[sp]]
        SA_matdata500[SA500_m_SA401[[sp]], "match_SA401"] <- SA401_m_SA401[[sp]]
        SA_matdata401[SA401_m_SA401[[sp]], "check"] <- matstring401[[sp]] != matstring500[[sp]]
        SA_matdata500[SA500_m_SA401[[sp]], "check"] <- matstring401[[sp]] != matstring500[[sp]]
    }
}
# fix variables in SA_matdata500
SA_matdata500[!SA_matdata500$SA_in401, 
              c("numequal", "splitpos", "match_SA401", "check")] <- NA

## This gives three data frames:
# 1. SA_sppdata: each row is a unique SpeciesAuthor from 2.0.1.
# variables are: 
# - unique_SA401: the SpeciesAuthor from 2.0.1
# - num401: the number of matrices for that SA in 2.0.1
# - num500: the number of matrices for that SA in 3.0.0
# - numequal: are the 2 above variables equal?
# - consec401: are the matrices for that SA on consecutive rows in 2.0.1?
# - consec500: are the matrices for that SA on consecutive rows in 3.0.0?
# - check1: do the matrices for that SA need checking to match them up?
#   (this is TRUE if any of the above 3 variables are FALSE)
# 2. SA_matdata401: each row is unique matrix from 2.0.1.
# variables are:
# - SA401: the SpeciesAuthor from 2.0.1
# - SA_in500: is the SpeciesAuthor in 3.0.0?
# - numequal: are the numbers of matrices for the SA equal in 2.0.1 and 3.0.0?
# - splitpos: are the matrices on consecutive rows in 401 / 500?
# - match_SA500: corresponding rows for matrix in 3.0.0 (if numequal is FALSE
# - and splitpos is FALSE)
# - check: do the matrices for that SA need checking to match them up?
#   (this is TRUE if numequal or splitpos are FALSE, or if
#   the string-converted matrices don't match in the 2 DBs).
# 3. SA_matdata500: same as SA_matdata401 but with SpeciesAuthor from 3.0.0
# (SA500), SA found in 401? (SA_in401), and corresponding rows for matrix in
# 2.0.1 (if numequal is FALSE and splitpos is FALSE) (match_SA401). For any
# SA not found in 401, all variables apart from SA500 and SA_in401 are NA.


# Any problems?
all(cp401_SAdata$consec401) #FALSE
all(cp401_SAdata$consec500) #FALSE
# ...yes

# which species have different numbers in different DBs?
cp_equalF <- data.frame(which = which(!cp401_SAdata$numequal),
                        SA = unique_SA401[!cp401_SAdata$numequal])
cp_equalF
# 2   3   6   9  13  28  31  32  45  50  51  52  62  64  66  67  69  72
# 82  90  92  94  97 102 107 108 110 111 112 122 129 130 134 136 137 145
# 165 168 170 172 173 176 180 183 204 212 213 214 217 219 228 231 235 243
# 245 252 253 254 256 261 272 275 287 291 293 302 307 317 323 336 337 338
# 342 343 346 347 348 351 357 364 373 386 395 428 436 439 440 462 463 472
# 478 480 483 484 489 491 492 511 518 541 548 549 550 557 559 560 564 565
# 570 573 574 579 584 588 599 606 607 612 615 638 643 668 670 672 673 687
# 713 714 716 727 741 742 743 753 773# are any of the SA not in 3.0.0 at all?

which(cp401_SAdata$num500 %in% 0)
# 32 111 176 228 291 357 463 548 615 713 714 727

# Which of the SA matrices aren't consecutive in 401
consecF401 <- data.frame(which = which(!(cp401_SAdata$consec401 %in% TRUE)),
                         consec401 = cp401_SAdata$consec401[!(cp401_SAdata$consec401 %in% TRUE)],
                         SA = unique_SA401[!(cp401_SAdata$consec401 %in% TRUE)])
consecF401
# 213

# Which aren't consecutive in 500
consecF500 <- data.frame(which = which(!(cp401_SAdata$consec500 %in% TRUE)),
                         consec500 = cp401_SAdata$consec500[!(cp401_SAdata$consec500 %in% TRUE)],
                         SA = unique_SA401[!(cp401_SAdata$consec500 %in% TRUE)])
consecF500
# 3   32  50  51  52  94  102 110 111 165 168 176 204 228 245 254 261 275 291
# 323 357 364 436 463 478 480 484 548 559 573 574 579 584 612 615 713 714 727
# 741 743

# Which species need all their matrices checking?
cp401_SAcheck <- data.frame(which = which(cp401_SAdata$check1 %in% TRUE),
                      check1 = cp401_SAdata$check1[cp401_SAdata$check1 %in% TRUE],
                      SA = unique_SA401[cp401_SAdata$check1 %in% TRUE])
cp401_SAcheck
# 2   3   6   9   13  28  31  32  45  50  51  52  62  64  66  67  69  72
# 82  90  92  94  97  102 107 108 110 111 112 122 129 130 134 136 137 145
# 165 168 170 172 173 176 180 183 204 212 213 214 217 219 228 231 235 243
# 245 252 253 254 256 261 272 275 287 291 293 302 307 317 323 336 337 338
# 342 343 346 347 348 351 357 364 373 386 395 428 436 439 440 462 463 472
# 478 480 483 484 489 491 492 511 518 541 548 549 550 557 559 560 564 565
# 570 573 574 579 584 588 599 606 607 612 615 638 643 668 670 672 673 687
# 713 714 716 727 741 742 743 753 773

# how many matrices in total need checking?
length(which(SA_matdata500$check %in% TRUE))
# 1897






### COMADRE ####################################################################

# load comadre development, current and legacy versions
load("DBVersionsX/COMADRE_v.X.X.X.Rdata")
cmxxx <- as_cdb(comadre)
rm(comadre)
load("DBVersions/COMADRE_v.3.0.0.Rdata")
cm300 <- as_cdb(comadre)
rm(comadre)
load("DBVersions/COMADRE_v.2.0.1.Rdata")
cm201 <- as_cdb(comadre)
rm(comadre)

( cmxxx_vers <- table(cmxxx$VersionRelease, useNA = "always") )
names(cmxxx_vers)[is.na(names(cmxxx_vers))] <- "NA"
sum(cmxxx_vers) - sum(cmxxx_vers[c("X", "NA")]) # 2281
nrow(cm300) # 2277
# website: 2277

# match up the SpeciesAuthor from 2.0.1 to 3.0.0
# create a numeric (2.0.1 dim) which will be the indices used to match
# create a vector of TRUE (2.0.1 dim) which will be used to flag problem matrices
# for each SpeciesAuthor... (obtained using unique(SpeciesAuthor))
## check that all the matching numbers are consecutive (diff == 1)
### if they aren't consecutive, skip indices & flag all
### if they are, create consecutive index from smallest to largest number
###   returned, put in the right place and ignore flag

# dim for both versions
dim201 <- dim(cm201)[1]
dim300 <- dim(cm300)[1]

# SpeciesAuthor for both versions
SA201 <- cm201$SpeciesAuthor
SA300 <- cm300$SpeciesAuthor
# unique SA for 2.0.1
unique_SA201 <- unique(SA201)
num_SA201 <- length(unique_SA201)

# list of matrices for both versions
matA201 <- matA(cm201)
matA300 <- matA(cm300)

### set up for loop

## lists
# list of integers to place each SA in 201
SA201_m_SA201 <- list()
# list of integers to match SA from 201 to 300
SA300_m_SA201 <- list()
# list of characters to hold each matrix string (201) by SA
matstring201 <- list()
# list of characters to hold each matrix string (300) by SA
matstring300 <- list()

## species data (length of unique SA values)
# dataframe of species info for SA in 201: number of matrices for that SA in 
# 201 and 300, whether matrices for that SA in 201 and 300 are consecutive rows,
# whether number of matrices in for that SA 201 and 300 are equal to one another
cm201_SAdata <- data.frame(unique_SA201,
                         num201 = numeric(num_SA201), num300 = numeric(num_SA201),
                         numequal = logical(num_SA201),
                         consec201 = logical(num_SA201), consec300 = logical(num_SA201),
                         check1 = logical(num_SA201)
)

## matrix data (length of num matrices in 2.0.1 and 3.0.0 respectively)
# dataframe of matrix infor for matrices in 201: 
# splitpos indicate whether matrices for that SA are consecutive in both 
# 201 and 300, 
# numequal indicates whether the number of matrices for that SA 
# is equal in 201 and 300, 
# match_SA300 is the index matching the that matrix to its supposed 
# corresponding row for that matrix in 300 (only for any matrix where check1 
# is false),
# check indicates whether the matrix needs checking in the 2nd round 
# (numequal is false, splitpos is true, &/or 
# character string version of matrix doesn't match in 201 and 300).
SA_matdata201 <- data.frame(SA201, SA_in300 = logical(dim201),
                            numequal = logical(dim201),
                            splitpos = logical(dim201),
                            match_SA300 = numeric(dim201),
                            check = logical(dim201)
)
# dataframe of matrix info for matrices in 300: 
# same as above, except match_SA201 reverse-matches the matrix to its supposed 
# corresponding row in 201.
SA_matdata300 <- data.frame(SA300, SA_in201 = logical(dim300),
                            numequal = logical(dim300),
                            splitpos = logical(dim300),
                            match_SA201 = numeric(dim300),
                            check = logical(dim300)
)

### loop to fill in all the above empty information 
### (sp = one SpeciesAuthor from 201)
for (sp in seq_along(unique_SA201)) {

    ## lists
    # positions of sp in 201
    SA201_m_SA201[[sp]] <- which(SA201 %in% unique_SA201[sp])
    names(SA201_m_SA201)[sp] <- unique_SA201[sp]
    # positions of sp in 300
    SA300_m_SA201[[sp]] <- which(SA300 %in% unique_SA201[sp])
    names(SA300_m_SA201)[sp] <- unique_SA201[sp]
    matstring201[[sp]] <- sapply(matA201[SA201_m_SA201[[sp]]], mat_to_string)
    names(matstring201)[sp] <- unique_SA201[sp]
    matstring300[[sp]] <- sapply(matA300[SA300_m_SA201[[sp]]], mat_to_string)
    names(matstring300)[sp] <- unique_SA201[sp]

    ## Species data frame
    # how many matrices for sp in each of the two versions?
    cm201_SAdata[sp, c("num201", "num300")] <-
        c(length(SA201_m_SA201[[sp]]), length(SA300_m_SA201[[sp]]))
    # are the number of matrices for sp in each data base equal?
    cm201_SAdata$numequal[sp] <- cm201_SAdata[sp, "num201"] == cm201_SAdata[sp, "num300"]
    # are all of the index values consecutive?
    # (If yes hopefully the matrices can be matched over in chunks)
    # Note: species with only 1 matrix are automatically consecutive
    #       species with 0 matrix are NA
    if (cm201_SAdata[sp, "num201"] > 1) {
        cm201_SAdata[sp, "consec201"] <- all(diff(SA201_m_SA201[[sp]]) == 1)
    }
    if (cm201_SAdata[sp, "num300"] > 1) {
        cm201_SAdata[sp, "consec300"] <- all(diff(SA300_m_SA201[[sp]]) == 1)
    }
    if (cm201_SAdata[sp, "num201"] == 1) cm201_SAdata[sp, "consec201"] <- TRUE
    if (cm201_SAdata[sp, "num300"] == 1) cm201_SAdata[sp, "consec300"] <- TRUE
    if (cm201_SAdata[sp, "num201"] == 0) cm201_SAdata[sp, "consec201"] <- NA
    if (cm201_SAdata[sp, "num300"] == 0) cm201_SAdata[sp, "consec300"] <- NA
    # First check: are there equal numbers of matrices for sp in both 201 and
    # 300, with all matrices on consecutive rows in both data bases?
    cm201_SAdata[sp, "check1"] <- !all(c(cm201_SAdata$numequal[sp],
                                       cm201_SAdata$consec201[sp],
                                       cm201_SAdata$consec300[sp]))

    ## matrix data frames
    # is the species in the other data base?
    SA_matdata201[SA201_m_SA201[[sp]], "SA_in300"] <- cm201_SAdata[sp, "num300"] > 0
    SA_matdata300[SA300_m_SA201[[sp]], "SA_in201"] <- cm201_SAdata[sp, "num201"] > 0
    # are the number of matrices for each sp equal in 201 and 300?
    SA_matdata201[SA201_m_SA201[[sp]], "numequal"] <- cm201_SAdata$numequal[sp]
    SA_matdata201[SA300_m_SA201[[sp]], "numequal"] <- cm201_SAdata$numequal[sp]
    # are the matrices on consecutive rows in 201 / 300?
    SA_matdata201[SA201_m_SA201[[sp]], "splitpos"] <-
        any(cm201_SAdata$consec201[sp], cm201_SAdata$consec300[sp])
    SA_matdata300[SA300_m_SA201[[sp]], "splitpos"] <-
        any(cm201_SAdata$consec201[sp], cm201_SAdata$consec300[sp])
    # 1. match the rows for sp in 300 to the rows in 201 and vice versa
    #    (NA if check1 is TRUE as there is no easy match)
    # 2. check: TRUE if check1 is TRUE, also TRUE if check1 is FALSE but 
    #    the matrices, represented as strings, are not the same.
    if (cm201_SAdata$check1[sp]) {
        SA_matdata201[SA201_m_SA201[[sp]], "match_SA300"] <- 0
        SA_matdata300[SA300_m_SA201[[sp]], "match_SA201"] <- 0
        SA_matdata201[SA201_m_SA201[[sp]], "check"] <- TRUE
        SA_matdata300[SA300_m_SA201[[sp]], "check"] <- TRUE
    }
    if (!cm201_SAdata$check1[sp]) {
        SA_matdata201[SA201_m_SA201[[sp]], "match_SA300"] <- SA300_m_SA201[[sp]]
        SA_matdata300[SA300_m_SA201[[sp]], "match_SA201"] <- SA201_m_SA201[[sp]]
        SA_matdata201[SA201_m_SA201[[sp]], "check"] <- matstring201[[sp]] != matstring300[[sp]]
        SA_matdata300[SA300_m_SA201[[sp]], "check"] <- matstring201[[sp]] != matstring300[[sp]]
    }
}
# fix variables in SA_matdata300
SA_matdata300[!SA_matdata300$SA_in201, 
              c("numequal", "splitpos", "match_SA201", "check")] <- NA

## This gives three data frames:
# 1. cm201_SAdata: each row is a unique SpeciesAuthor from 2.0.1.
# variables are: 
# - unique_SA201: the SpeciesAuthor from 2.0.1
# - num201: the number of matrices for that SA in 2.0.1
# - num300: the number of matrices for that SA in 3.0.0
# - numequal: are the 2 above variables equal?
# - consec201: are the matrices for that SA on consecutive rows in 2.0.1?
# - consec300: are the matrices for that SA on consecutive rows in 3.0.0?
# - check1: do the matrices for that SA need checking to match them up?
#   (this is TRUE if any of the above 3 variables are FALSE)
# 2. SA_matdata201: each row is unique matrix from 2.0.1.
# variables are:
# - SA201: the SpeciesAuthor from 2.0.1
# - SA_in300: is the SpeciesAuthor in 3.0.0?
# - numequal: are the numbers of matrices for the SA equal in 2.0.1 and 3.0.0?
# - splitpos: are the matrices on consecutive rows in 201 / 300?
# - match_SA300: corresponding rows for matrix in 3.0.0 (if numequal is FALSE
# - and splitpos is FALSE)
# - check: do the matrices for that SA need checking to match them up?
#   (this is TRUE if numequal or splitpos are FALSE, or if
#   the string-converted matrices don't match in the 2 DBs).
# 3. SA_matdata300: same as SA_matdata201 but with SpeciesAuthor from 3.0.0
# (SA300), SA found in 201? (SA_in201), and corresponding rows for matrix in
# 2.0.1 (if numequal is FALSE and splitpos is FALSE) (match_SA201). For any
# SA not found in 201, all variables apart from SA300 and SA_in201 are NA.




# Any problems?
all(cm201_SAdata$consec201) #FALSE
all(cm201_SAdata$consec300) #FALSE
# ...yes

# which species have different numbers in different DBs?
cm_equalF <- data.frame(which = which(!cm201_SAdata$numequal),
                        SA = unique_SA201[!cm201_SAdata$numequal])
cm_equalF
# 32  33 126 133 289 298 306 315 322 363 393 404 440 444 472

# are any of the SA not in 3.0.0 at all?
which(cm201_SAdata$num300 %in% 0)
# 33

# Which of the SA matrices aren't consecutive in 201
consecF201 <- data.frame(which = which(!(cm201_SAdata$consec201 %in% TRUE)),
                         consec201 = cm201_SAdata$consec201[!(cm201_SAdata$consec201 %in% TRUE)],
                         SA = unique_SA201[!(cm201_SAdata$consec201 %in% TRUE)])
consecF201
# 186 208

# Which aren't consecutive in 300
consecF300 <- data.frame(which = which(!(cm201_SAdata$consec300 %in% TRUE)),
                         consec300 = cm201_SAdata$consec300[!(cm201_SAdata$consec300 %in% TRUE)],
                         SA = unique_SA201[!(cm201_SAdata$consec300 %in% TRUE)])
consecF300
# 33 186 208 298 306 315 322 393 404

# Which species need all their matrices checking?
cm201_SAcheck <- data.frame(which = which(cm201_SAdata$check1 %in% TRUE),
                      check1 = cm201_SAdata$check1[cm201_SAdata$check1 %in% TRUE],
                      SA = unique_SA201[cm201_SAdata$check1 %in% TRUE])
cm201_SAcheck
# 32  33 126 133 186 208 289 298 306 315 322 363 393 404 440 444 472

# how many matrices in total need checking?
length(which(SA_matdata300$check %in% TRUE))
# 164


