library(dplyr)
library(stringr)
library(Rcompadre)
library(hashids)

load("DBVersions/COMPADRE_v.5.0.0.Rdata")
compadre_v.5.0.0 <- as_cdb(compadre)
rm(compadre)
load("DBVersions/COMPADRE_v.3.2.1.Rdata")
compadre_v.3.2.1 <- as_cdb(compadre)
rm(compadre)

# strip out rows where SpeciesAccepted is NA
cp500 <- subset(compadre_v.5.0.0, !is.na(compadre_v.5.0.0$SpeciesAccepted))
cp321 <- subset(compadre_v.3.2.1, !is.na(compadre_v.3.2.1$SpeciesAccepted))

# short (3&3) species names
cp500_short <- sapply(cp500$SpeciesAccepted, 
                      function(sp){
                          spbin <- str_split(sp, " ")[[1]][1:2]
                          sp3 <- str_sub(spbin, 1, 3)
                          sp3[is.na(sp3)] <- "XXX"
                          sp3x <- str_replace_all(sp3, "[^A-Za-z]", "X")
                          out <- paste(sp3x, collapse = "")
                          out
                      })

# get matrices and convert into strings of numbers ready for hashid
cp500_matA <- matA(cp500)
cp500_matA12 <- lapply(cp500_matA, round, digits = 12) # 12 is an arbitrary number
cp500_matAstr <- sapply(cp500_matA12, mat_to_string)
cp500_integer <- str_replace_all(cp500_matAstr, "[^1-9]", "")
### *** the integers aren't unique, even when combined with the species...

# create smaller integer strings so hashid doesn't break
integer_length <- str_length(cp500_integer)
short_length <- 11 # 11 (arbitrary) is the length of integer used to get hashid
length_by <- integer_length %/% short_length
cp500_integer_short <- mapply(function(int, il, lb){ 
                                  ith <- seq(1, il, lb)
                                  paste0(mapply(str_sub, int, ith, ith), collapse="")
                               }, 
                               int = cp500_integer, 
                               il = integer_length,
                               lb = length_by)

# create hashids
hid_set <- hashid_settings(salt = "CompadreDB", min_length = 12)
cp500_id <- paste(cp500_short, hashid(cp500_integer_short))

### *** the above code hasn't been tested, but should create a short(ish) ID for
### each matrix. These won't be unique, given the note above: problem is that 
### there are plenty of cases where the species and the matrix are the same, so
### creating the same ID for different rows of the data.

### *** need to repeat for comadre

load("DBVersions/COMADRE_v.3.0.0.Rdata")
comadre_v.3.0.0 <- as_cdb(comadre)
rm(comadre)
load("DBVersions/COMADRE_v.2.0.1.Rdata")
comadre_v.2.0.1 <- as_cdb(comadre)
rm(comadre)

# shorter names
cm300 <- comadre_v.3.0.0
cm201 <- comadre_v.2.0.1
