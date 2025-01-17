## Charge inversion and nanoscale phenomena for living cells ##

Giant charge inversion and/or nanoscale phenomena for living cells are studied 
using the codes of fortran 2003 and MPI, and molecular dynamics simulations are execcuted (Refs. 1-5).
Especially, "DNA in nanopores" is simulated for DNA transportation through nanepores, where 
counterion condensation and coion depletion are the keys of the nanoscale pores of human cells.

The molecular dynamics simulation codes are created in several settings, and they are
discussed in the articles. Gaint charge inversion is about the macroions 
surrounded by the electrolyte of multivalent counterions and monovalent coions.
The simulation code is @chginv3.f03 with the paramer file parm_inv13.h and 
the configure file CIMV13_config.START3.
The equation of motion in Eq.(1) of Ref.1 has the Lengevin thermostat on top of 
Coulomb and Lennard-Jones forces of the righthand side.
The first 40 lines of the code are the title, author, equation of motion, and 
structure of the code. Other important remarks are explained thereafter.
Main subroutines are RUN_MD, moldyn, realteil, p3m_perform, and Gopen graphic packages.

The simulation code utilizes the MPI and FFTW3 packages, and PDF's of Ref. 1 are shown 
with figures. We can see the occurrence of giant charge inversion by looking at 
sharp peaks in charge distribution functions around the macroions. 
From the figures, there are two requirements for the giant inversion, namely,
(1) the multivalent counterions and monovalent coions exist around macroions, and 
(2) the Coulomb energy is larger than thermal energy of macroions and counterions.

In the listed articles, the figures are best illustrated, partly 
in color pictures. We use counterions with Zp=3 and the radius 2.0 Ang and coions 
with Zn=-1 and 5.2 Ang. The PDF file cimv13_773a.pdf in the Langevin thermostat 
shows a rather stationary macroion. The first peak of counterions at the top-left 
panels exists very close to the macroion's surface at R= Rmac to (Rmac+5 Ang) of 
cimv13_78.pdf - only the 20 Angstrom regime from the macroion is plotted. 

(* CGS system: t_unit= 0.01d-12 s, a_unit= 1.00d-08 cm, w_unit= 1.6605d-24 g, 
e_unit= 4.803d-10 esu, and \epsilon=78 at 27 deg Celsius.)

Everyone is welcome to copy and/or rewrite this simulation code under GPL-3
license, by keeping initial 40 lines of this author of the code intact. 

## References ##
1. M. Tanaka and A.Yu Grosberg, J.Chem.Phys., vol.115, 567-574 (2001).
2. M. Tanaka and A.Yu Grosberg, Euro.Phys.J., E7, 371-379 (2002).
3. M. Tanaka, Phys.Reviews., E68, 061501 (2003).
4. M. Tanaka, J. Physics: Condensed Matter, vol.16, S2127-2134 (2004).
5. Y.Rabin and M.Tanaka, Phys.Rev.Lett., vol.94, 148103 (2005).


