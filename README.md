com(p)adreDB
==========

This fork is mainly for my own use, but contains some useful code and pointers
to certain 'upgrades' available for the databases. The data themselves aren't 
available online in this remote (they're not 'official' versions of the database),
but can be obtained by emailing me iainmstott[at]gmail.com.

***

This repository contains code for the COMPADRE Plant Matrix Database and COMADRE Animal Matrix Database.
It will evolve into an R package but will initially contain scripts for interacting with the COMPADRE and COMADRE databases which are released as structured R list objects.

The data can be downloaded from www.compadre-db.org, or www.comadre-db.org.
Follow news at our blog https://compadredb.wordpress.com, and our Twitter accounts @compadreDB @comadreDB.


The structure of the compadre and comadre data objects
-----------------------------------------
    compadre/comadre +-- metadata {dataframe with ca 60 columns and one row one for each set of matrices (mat A, matU, matF, matC)
             |
             +-- matrixClass {list with one entry for each set of matrices. Each entry is a data frame with 3 columns: MatrixClassOrganized, MatrixClassAuthor, MatrixClassNumber.}
             |
             |-- mat {list with one entry for each row of metadata}
             |     |
             |     +-- matA {matrix}
             |     +-- matU {matrix}
             |     +-- matF {matrix}
             |     \-- matC {matrix}
             |
              \-- version {a vector with version information}

