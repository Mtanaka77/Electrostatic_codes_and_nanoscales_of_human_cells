## Molecular dynamics simulation for electrostatic living cells ##

The charge inversion and ion transport phonemena through nanopores are studied for living cells. 
The electrostatic molecular dynamics simulations are execcuted (Refs. 1-5).
We first talk about the electostatic code to understand and properly execute molecular dynamics simulation. 
Then back to the thema, the charge inversion and "DNA in nanopores" are simulated for DNA transport through nanopores.
The former deals with the periodic boundary system, while the latter treats the 3D non-periodic boundary problem with short-range and long-range electrostatic interactions. 
The counterion condensation and coion depletion are the keys of the nanoscale human cells.


### (1) What is electrostatic molecular dynamics simulation ###

First of all, we have an easy electrostatic simulation code of three dimensions that first does 
Fortran 2003 compilation, and then does parallel execution. It uses MPI v.3 and the Ghostviewer script. 
We use the Linux operating system of 6 cores of 3 GMz, typically in our desktop workstation. 

The program is divided with, (i) the parallelization and parameters setups, (ii) the initialization of
ions and electrons by /init/, (iii) the main loop of simulation run /moldyn/, where the important 
subroutine is the calculation of ions and electrons dynamics, /cl_forces/,
and (iv) closing of the code /gclose/. The mass of ion is 1836 for that of electron 1. 
The open boundary system is used for which the periodic boundary system /reflect/ is skipped.
The run takes about one hour, and graphic outputs are generated in time as expl07.77a.ps.
It should be converted to expl07.77a.pdf to view on the screen or draw on letter paper.

We can see in the figure that ions are spattered but they stay roughly at constant radii,
where electrons tend to have the similar trend because Coulomb forces make the positive and negative 
species togather in the sphere.

Code: @md3-para7.f03 with the parameter file paramE7.h, ca. 1600 lines.


### (2) Simulation of charge inversion phenomena ###

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


### (3) DNA transport through nanepores ###

The transport of DNA with counterions and coions is studied where a narrow nanopore along the z-direction seperates wide downside and upside regions (Ref. 5).  The cylinder of the pore is assumed 1.5 nm wide and 4.0 nm high, embedded in 4.0 nm in x,y directions and 12.0 nm in the z direction. The short-range Coulomb and Lennard-Jones forces are treated, i.e., 
m dv/dt = (q'q/Gamma *r^2) (grad r/r) - fgm *(2 r(i)-r(i+1)-r(i-1)) +48 *(epsil_LJ/kT) grad[(sigma/r_ij)^12 -(sigma/r_ij)^6].
The long-range forces with the meshes of (i,j,k) coordinates are solved by the Poisson equation, i.e.,
div(epsilon grad [pot(i,j,k)]) = - 4*pi *rho(i,j,k). 

The simulation code is named @nanoporAPF.f03 (ca 9,900 lines with graphics), and the paramter file paramAPF.h and the configuration file PORV41_config.start3. The used subroutines are: RUN_MD, moldyn, sht-forces, LJ-forces, sprmul, reflect_endpl, init, poissn, emcof3, cresmd, and graphics. It has N_x=N_y=80 and N_z=120 meshes, 14,000 particles, and a test run takes 15 minutes/6 cores (3.0 GHz) for t=800 with dt=0.01 (x 10^-14 s).  

The file porv41.777.pdf shows four plots of potential plot, particles of DNA and ions, those of all particles 
(every 100 of water), and the velocity plot. One can see the DNA moving toward the positive z direction with time
to the cell.  



### Coulomb and P3M or bound (non periodic) simulation code ###

The short-range forces and the long-range forces in p3m_perform routines, or non-periodic code are included 
to study the charge inversion or the DNA tranport phenomena at high accuracy (Refs. 1, 5).
Everyone is welcome to copy and rewrite this simulation codes under 
GNU General Public License v3.0, by keeping top (2) 55 lines or (3) 100 lines of the
author's codes intact. 


## References ##
1. M. Tanaka and A.Yu Grosberg, J.Chem.Phys., vol.115, 567-574 (2001).
2. M. Tanaka and A.Yu Grosberg, Euro.Phys.J., E7, 371-379 (2002).
3. M. Tanaka, Phys.Reviews., E68, 061501 (2003).
4. M. Tanaka, J. Physics: Condensed Matter, vol.16, S2127-2134 (2004).
5. Y. Rabin and M. Tanaka, Phys.Rev.Lett., vol.94, 148103 (2005).
6. M. Tanaka and M. Murakami, Comp.Phys.Commun., 241, pp. 56-63 (2019).
7. M. Deserno and C. Holm, J.Chem.Phys. 109, 7694–7701 (1998).
8. M. Tanaka, Collaboration at Institute of Polymerforchung, University of Mainz, Germany (1999).

