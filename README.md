## Charge inversion and nanoscale phenomena for living cells ##

Giant charge inversion and/or nanoscale phenomena for living cells are studied 
using fortran 2003, and molecular dynamics simulations are execcuted (Refs. 1-5).
"DNA in nanopores" in Ref. 5 is simulated for DNA transportation, where 
counterion condensation and coion depletion are the key of the nanoscale pores 
of human cells.

The molecular dynamics simulation codes are created in several settings, and they are
discussed in the articles. Gaint charge inversion is about the macroions 
surrounded by the electrolyte of multivalent counterions and monovalent coions.
The simulation code is @chginv5.f03 with the paramer file parm_inv15.h and 
the configure file CIMV15_config.START3.
The equation of motion has the Lengevin thermostat on top of Coulomb and Lennard-Jones
forces, as are the first 40 lines of the title, author, equation of motion, and 
structure of the code. Other important remarks are explained thereafter.
Main subroutines are RUN_MD, moldyn, realteil, p3m_perform, and Gopen graphic packages.

The simulation run is utilized in terms of MPI and FFTW3 packages, and PDF's 
of Ref. 1 are shown with figures. We can see the occurrence of charge inversion 
by looking at sharp peaks in charge distribution functions around the macroions. 
From the figures, there are two requirements for the giant inversion, namely,
(1) the multivalent counterions and monovalent coions exist around macroions, and 
(2) the Coulomb energy is bigger than thermal energy of such macroions.

In the articles (Refs. 1-5), the figures are best illustrated, partly 
in color pictures. The PDF file cimv15_773a.pdf shows a stationary macroion,
for which the first peak of counterions exists by cimv15_78.pdf close to the 
macroion's surface.

## References ##
1. M. Tanaka and A.Yu Grosberg, J.Chem.Phys., vol.115, 567-574 (2001).
2. M. Tanaka and A.Yu Grosberg, Euro.Phys.J., E7, 371-379 (2002).
3. M. Tanaka, Phys.Reviews., E68, 061501 (2003).
4. M. Tanaka, J. Physics: Condensed Matter, vol.16, S2127-2134 (2004).
5. Y.Rabin and M.Tanaka, Phys.Rev.Lett., vol.94, 148103 (2005).


