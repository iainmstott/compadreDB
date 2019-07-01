Versions of the COMPADRE Plant Matrix Database and COMADRE Animal Matrix Database 'cleaned up'
==============================================================================================

This fork is mainly for my own use, but contains some useful code and pointers
to certain 'upgrades' available for the databases. The data themselves aren't 
available online in this remote (they're not 'official' versions of the database),
but can be obtained by emailing me iainmstott[at]gmail.com.

***

This folder contains code to 'clean up' the com(p)adre databases. It should 
run on any database version (and both current and previous versions are 
supplied in this folder).

The code adds several variables to the metadata which can be used to exclude matrices 
that suffer from a number of common problems. The problems cover issues that 
can interfere with coding, violations of model assumptions, and potential issues with 
how the original data have been used to build the model. All the problems addressed 
here can be calculated from the raw matrices, and yield the following variables:
* **hasNA**: does the matrix contain NAs? (TRUE / FALSE)
* **Irreducible**: is the matrix irreducible? (TRUE / FALSE) - in most
                   cases, it's desirable that the matrix is irreducible 
                   (Irreducible == TRUE)

* **AvSurvRYG**: is there averaging of survival in the matrix? (_R_ / _Y_ / _G_) 
it's common to average vital rates across ages or stages where individuals are 
hard to distinguish, individuals are hard to identify by age, or otherwise 
detailed information isn't available.
    - _R_: Red, where survival (column sums of the U matrix) has the same 
       average value consecutively across > half of the (st)ages.
    - _Y_: Yellow, where survival (column sums of the U matrix) has the same 
       average value consecutively across >= 2 (st)ages (i.e. there is at 
       least some averaging), but identical consecutive values are never more than
       half the number of (st)ages (i.e. a 3x3 matrix with 2 identical consecutive 
       survival values is red, whereas a 4x4 matrix with 2 identical consecutive survival 
       values is yellow).
    - _G_: Green, where there is no averaging (all survival values are different).
* **AvFecRYG**: is there averaging of fecundity in the matrix? (_R_ / _Y_ / _G_) 
it's common to average vital rates across ages or stages where individuals are 
hard to distinguish, individuals are hard to identify by age, or otherwise 
detailed information isn't available.
    - _R_: Red, where fecundity (column sums of the F matrix) has the same 
       average value consecutively across > half of reproductive (st)ages.
    - _Y_: Yellow, where fecundity (column sums of the F matrix) has the same 
       average value consecutively across >= 2 of reproductive (st)ages (i.e. there is at 
       least some averaging), but identical consecutive values are never more than
       half the number of (st)ages (i.e. a matrix with 2 identical consecutive 
       fecundity values of a total of 3 reproductive stages is red, whereas a matrix with 2 identical consecutive fecundity values of a total of 4 reproductive stages is yellow).
    - _G_: Green, where there is no averaging (all fecundity values are different).

The following variables are used to determine the **AvSurvRYG** and **AvFecRYG**
variables, and may be useful for coming up for bespoke criteria for deciding 
whether to include matrices with averaged vital rates:
* **UcolSum**: list column. Column sums of U matrices, i.e. survival of each stage.
* **AvSurvGrp**: groupings of identical consecutive survival values. E.g. if the stage 
  survival of an organism is `c(0.5, 0.5, 0.7, 0.8, 0.8, 0.8)` then the groups would 
  be recorded as `c(1, 1, 2, 3, 3, 3)`.
* **AvSurvMax**: the maximum consecutive number of stages with identical survival
  values. In the above example this would be 3.
* **FcolSum**: list column. Column sums of F matrices, i.e. fecundity of each stage.
* **FirstF**: the number of the first reproductive stage.
* **FnStages**: the number of reproductive stages.
* **AvFecGrp**: groupings of identical consecutive fecundity values. E.g. if the stage 
  fecundity of an organism is `c(0.2, 0.5, 0.5, 1.2, 1.2)` then the groups would 
  be recorded as `c(1, 2, 2, 3, 3)`.
* **AvFecMax**: the maximum consecutive number of stages with identical fecundity
  values. In the above example this would be 2.

Perhaps these extra variables will be useful for working with the data. The code 
provided should work with any past database versions, and any subsets of the 
database.  

Possible variabes to be included in the future:
* **inPhylo**: is the matrix in the phylogeny? (TRUE / FALSE)

(to be completed)
