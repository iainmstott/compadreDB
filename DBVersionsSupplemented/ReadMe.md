Versions of the COMPADRE Plant Matrix Database and COMADRE Animal Matrix Database supplemented with extra data
==============================================================================================================

This fork is mainly for my own use, but contains some useful code and pointers
to certain 'upgrades' available for the databases. The data themselves aren't 
available online in this remote (they're not 'official' versions of the database),
but can be obtained by emailing me iainmstott[_at_]gmail[_dot_]com.

***

This folder includes code to help 'clean up' the databases, by adding variables
that indicate problematic models. 

The below info is available for both compadre and comadre up to versions 
3.2.1 and 2.0.1 respectively. Average matrices calculated by the com(p)adre team
are removed as of versions 5.0.0 and 3.0.0 respectively.

Extra data includes:
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


The below info is available only for compadre version 3.2.1.

**MetadataSupplement.csv**
* **OrganismType** (for compadre version 3.2.1 and lower): the life history categories defined compadre 4.0.0 and comadre 2.0.0 and higher. For information on this variable please see the [compadre user guide](https://github.com/jonesor/compadreDB/blob/master/COMPADRE-UserGuide/COMPADRE-UserGuide.pdf)

* **SuccessionType**: the life history categories defined in Silvertown +_et al_.
                     (1993) J. Ecol, 81, 465-476. These are:
  - _Algae_: Algae
  - _Monocarp_: Monocarpic (semelparous) plants
  - _Open_: iteroparous herbs of open habitats
  - _Closed_: iteroparous herbs of closed (e.g. forest) habitats
  - _Shrub_: Shrubs
  - _Tree_: Trees

* **Endangered**: whether the species is endangered or not, based on [IUCN red list category](http://www.iucnredlist.org/about/introduction) This may include:
  - _TRUE_ (EX, )
  - Extinct in the wild (_EW_)
  - Critically Endangered (_CR_)
  - Endangered (_EN_)
  - Vulnerable (_VU_)
  - Near Threatened (_NT_)
  - Least Concern (_LC_)
  - Data Deficient (_DD_)
  - Not Evaluated (_NE_)

**SeedsSupplement.csv**
* **SeedsProblem**: does the matrix suffer from the 'seeds problem'? 
                    (TRUE / FALSE) - in most cases, it's desirable that the
                    matrix does not have the seeds problem.



At the moment this data is only collected for plants (compadre) and available for version 3.2.1.