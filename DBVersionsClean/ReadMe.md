Versions of the COMPADRE Plant Matrix Database and COMADRE Animal Matrix Database 'cleaned up'
==============================================================================================

This folder contains code to 'clean up' the com(p)adre databases. Not all data 
is perfect, and sometimes models have certain problems or defy certain 
assumptions. The problems addressed here can all be calculated from the raw 
matrices, and yield the following variables:
* **hasNA**: does the matrix contain NAs? (TRUE / FALSE)
* **Irreducible**: is the matrix irreducible? (TRUE / FALSE) - in most
                   cases, it's desirable that the matrix is irreducible 
                   (Irreducible == TRUE)
* **AverageSurv**: is there averaging of survival in the matrix? (_R_ / _Y_ / _G_) 
it's common to average vital rates across ages or stages where individuals are hard to distinguish or 
detailed information isn't available.
    - _R_: Red, where survival (column sums of the U matrix) has the same 
       average value consecutively across >= half of the (st)ages.
    - _Y_: Yellow, where survival (column sums of the U matrix) has the same 
       average value consecutively across >= 3 (st)ages (i.e. there is at 
       least some averaging), but identical consecutive values never more than
       half the number of (st)ages.
    - _G_: Green, where there is no averaging.
* **MatrixCompositePlus**: extension to the MatrixComposite variable. The 
                           original variable doesn't distinguish between 
                           mean/pooled matrices which are calculated from 
                           other matrices in the database, and those which 
                           are calculated from matrices which _aren't_ in the 
                           database (and therefore in fact fulfil the criteria 
                           of being an individual matrix as outlined in the 
                           user documents). MatrixCompositePlus contains the 
                           variables:
    - _Individual_: Same as for MatrixComposite, representing a single study x 
                species x population x treatment x period
    - _IndividualMean_: A mean matrix which is calculated from other matrices              which _aren't_ in the com(p)adre database, and therefore is in             fact the sole matrix for its study x species x population x 
                treatment x period.
    - _MultipleMean_: A mean matrix which is calulated from other matrices                  which _are_ in the com(p)adre database.
    - _IndividualPooled_: A pooled matrix which is calculated from data used                 to build other matrices which _aren't_ in the com(p)adre 
                  database, and therefore is in fact the sole matric for its 
                  study x species x population x treatment x period.
    - _MultiplePooled_: A pooled matrix which is calculated from data used to build other matrices which _are_ in the com(p)adre database.
    - _Seasonal_: Same for MatrixComposite, a matrix that does not descrive a 
              full annual transition.
* **inPhylo**: is the matrix in the phylogeny? (TRUE / FALSE)

These are variables that have been useful for me in working with the data, 
and perhaps they'll be useful for you too.






~~SpeciesAuthor~~  
SpeciesAccepted  
~~CommonName~~  
~~Genus~~  
~~Family~~  
~~Order~~  
~~Class~~  
~~Phylum~~  
~~Kingdom~~  
~~OrganismType~~  
~~DicotMonoc~~  
~~AngioGymno~~  
Authors  
Journal  
YearPublication  
~~DOI.ISBN~~  
~~AdditionalSource~~  
~~StudyDuration~~  
~~StudyStart~~  
~~StudyEnd~~  
~~AnnualPeriodicity~~  
NumberPopulations  
~~MatrixCriteriaSize~~  
~~MatrixCriteriaOntogeny~~  
~~MatrixCriteriaAge~~  
MatrixPopulation  
~~Lat~~  
~~Lon~~  
~~Altitude~~  
~~Country~~  
~~Continent~~  
~~Ecoregion~~  
StudiedSex  
~~MatrixComposite~~  
MatrixTreatment  
MatrixCaptivity  
MatrixStartYear  
MatrixStartSeason  
MatrixStartMonth  
MatrixEndYear  
MatrixEndSeason  
MatrixEndMonth  
~~MatrixSplit~~  
~~MatrixFec~~  
~~Observation~~  
~~MatrixDimension~~  
~~SurvivalIssue~~  
