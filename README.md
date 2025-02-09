## Charge inversion and nanoscale phenomena for living cells ##

Giant charge inversion and nanoscale phenomena for living cells are studied 
using the Fortran 2003 code, and molecular dynamics simulations are execcuted (Refs. 1-5).
Especially, "DNA in nanopores" is simulated for DNA transportation through nanepores, where 
counterion condensation and coion depletion are the keys of the nanoscale human cells.

### Simulation of charge inversion phenomena ###

Charge inversion is about the macroions surrounded by the electrolyte of multivalent counterions 
and monovalent coions.
The molecular dynamics simulation code to study charge inversion is created in different settings, 
and is discussed in the articles.
The equation of motion in Eq.(1) of Ref.1 has the Lengevin thermostat on top of 
Coulomb and Lennard-Jones forces of the righthand side.
The short-range forces Eq.(1) and the long-range forces in p3m_perform routines are used 
for the best accuracy of molecular dynamics simulations (Ref. 6).

The simulation code is @chginv3.f03 with the paramer file parm_inv13.h and 
the configure file CIMV13_config.START3.
The first 55 lines of the code are the title, author, equation of motion, and 
the structure of the code. Other important remarks are explained thereafter.
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

### Short-range and long-range interactions ###

The Coulomb force by the positive charge $ Z_{p}*e_unit $ and negative charge $ Z_{n}*e_unit $
is written $ -\nabla(Z_{p}*Z_{n}*e_unit^2/ r) $ in an infinite system (Ref. 6).
However, the Coulomb forces in the periodic boundary system are calculated by 
the short-range and long-range forces for the maximum accuracy (Ref. 7).
The key of the short-range force is, 
$ prefactr * Z_{p}*Z_{n} *e_unit^2 (erfc/r +2 *alpha/sqrt(pi))*exp(-ar**2)/r $.
Such forces were formulated in the C language, which were rewritten for large speedup 
in the Fortran language (Ref. 8).

We note that the particle-in-cell method uses space meshes for averaged density in nuclear fusion, 
which is a different area from the molecular dynamics simulation.

(* CGS unit system: w_unit= 1.673d-24 g, e_unit= 4.803d-10 esu, t_unit= 1.00d-15 s, 
a_unit= 1.00d-08 cm, and \epsilon=78 at 27 deg Celsius.)

### Coulomb and P3M simulation code ###

The short-range forces and the long-range forces in p3m_perform routines are included 
to study the charge inversion phenomena at high accuracy (Ref. 1).
Everyone is welcome to copy and rewrite this simulation code under 
GNU General Public License v3.0, by keeping top 55 lines of the
author's code intact. 

## References ##
1. M. Tanaka and A.Yu Grosberg, J.Chem.Phys., vol.115, 567-574 (2001).
2. M. Tanaka and A.Yu Grosberg, Euro.Phys.J., E7, 371-379 (2002).
3. M. Tanaka, Phys.Reviews., E68, 061501 (2003).
4. M. Tanaka, J. Physics: Condensed Matter, vol.16, S2127-2134 (2004).
5. Y. Rabin and M. Tanaka, Phys.Rev.Lett., vol.94, 148103 (2005).
6. M. Tanaka and M. Murakami, Comp.Phys.Commun., 241, pp. 56-63 (2019).
7. M. Deserno and C. Holm, J.Chem.Phys. 109, 7694–7701 (1998).
8. M. Tanaka, Collaboration at Institute of Polymerforchung, University of Mainz, Germany (1999).

