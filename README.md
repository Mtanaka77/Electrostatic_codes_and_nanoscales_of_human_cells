## Charge inversion and nanoscale phenomena for living cells ##

Simple models of giant charge inversion or nanoscale phenomena are created 
for living cells, and molecular dynamics simulations are execcuted (Refs. 1-4).
"DNA in nanopores" is simulated for DNA transportation where counterion condensation 
and coion depletion are the key in the nanoscale pores of human cells (Ref. 5).

The molecular dynamics simulations are made in several settings, and they are
discussed in the articles. Gaint charge inversion in this case is about macroions 
surrounded by the electrolyte of multivalent counterions and monovalent coions.
The simulation code @chginv3.f03 with parminv3.h and CIMV3_config.START3 is used.
Description of the code is given, and main subroutines are moldyn, realteil and 
p3m_perform, 

The run is utilized by MPI and FFTW3 packages, and PDF's are shown with figures.
In the articles (Ref. 1), the figures are illustrated partly in color pictures.  


## References ##

1. M. Tanaka and A.Yu Grosberg, J.Chem.Phys., vol.115, 567-574 (2001).
2. M. Tanaka and A.Yu Grosberg, Euro.Phys.J., E7, 371-379 (2002).
3. M. Tanaka, Phys.Reviews., E68, 061501 (2003).
4. M. Tanaka, J. Physics: Condensed Matter, vol.16, S2127-2134 (2004).
5. Y.Rabin and M.Tanaka, Phys.Rev.Lett., vol.94, 148103 (2005).


