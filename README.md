## Molecular electrostatic dynamics for living human cells ##

To know the behavior of the long-range electrics (Poisson equation) in the living cells, 
the five-atom water model TIP5P is used (the Sec.3 below). You can skip Sec.1 and Sec.2 
in your choice.

The charge inversion and DNA transport phenomena are studied for living cells, and the 
molecular electrostatic dynamics is simulated (Refs. 1-5, 10).
We first mention the electrostatic code to execute molecular dynamics, 
and the charge inversion and finally the DNA transport are simulated.
The former deals with the periodic boundaries, while the latter treats the non-periodic 
boundaries with short-range and long-range electrostatic interactions in living cells. 

### (1) What is electrostatic molecular dynamics simulation ? ###

First of all, we have an easy electrostatic simulation code of three dimensions 
that first does Fortran 2003 compilation, and then does parallel execution. 
It uses MPI v.3 and the Ghostviewer script. 
We use the Linux operating system of 6 cores of 3 GMz, typically in our desktop workstation. 

The program is divided with, (i) the parallelization and parameters setups, (ii) the initialization of
ions and electrons by /init/, (iii) the main loop of simulation run /moldyn/, where the important 
subroutine is the calculation of ions and electrons dynamics, /cl_forces/,
and (iv) closing of the code /gclose/. The mass of ion is 1836 for that of electron is unity 1. 
The open boundary system is used for which the periodic boundary system /reflect/ is skipped.
The run takes about one hour, and graphic outputs are generated in time as expl07.77a.ps.
It should be converted to expl07.77a.pdf to view on the screen or draw on letter paper.

We can see in the figure that ions are spattered but they stay roughly at a constant radius.
Whereas electrons tend to initially expand faster but have the similar radius because 
Coulomb forces make the positive and negative species together within the sphere.

Code: @md3-para7.f03 with the parameter file paramE7.h, ca. 1600 lines.


### (2.A) Simulation of charge inversion phenomena ###

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

### (2.B) Short-range and long-range interactions ###

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


### (3) DNA molecule transport through nanopores ###

The transport of DNA with counterions and coions is studied, where a narrow nanopore in the z-direction seperates downside and upside regions (Ref. 5, Ref. 10 below). The cylinder of the pore is assumed 1.5 nm wide and 5.6 nm high, which is embedded in 7.4 nm in the x,y directions and 15.5 nm in the z direction. The short-range Coulomb and Lennard-Jones forces are treated, i.e., 
m_i dv_i/dt = Sum_i (Gamma q'q_i/r_ij^2) (grad r/r) - fgm *(2 r(i)-r(i+1)-r(i-1)) 
+48 *(epsil_LJ/kT) grad[(sigma/r_ij)^12 -(sigma/r_ij)^6] +q_i *E
with the coupling constant Bjerrum= Gamma * a_LJ =7 at T= 300 K.

The long-range potential forces with the mesh (i,j,k) coordinates are solved by the Poisson equation, i.e.,
div(eps(i,j,k) [grad pot(i,j,k)]) = - 4*pi* Gamma*rho(i,j,k). The operator positions of div, grad and eps(i,j,k) are important. There are large potentials of positive V_top and negative V_bot drops at the end plates, and they are small otherwise.

The simulation code is named @nanoAPGa.f03 (ca 9,500 lines with graphics), and the parameter file paramAPGa.h and the configuration file PORW21_config.start3. The used subroutines are: RUN_MD, moldyn, sht-forces, sprmul, reflecest_endpl, init, poissn, escof3, cresmd, and graphics. There are many input items to run the simulation; they are the nanopore sizes of R_pore and H_pore, respectively, the number of DNA, counterions and coions, the Cartesian meshes of the Poisson equation, the time step of dt, the potential values of top and bottom plates, and the Gamma value at 300 K, etc. It has N_x= N_y= 80 and N_z= 100 meshes for L_x= L_y= 74 Ang and L_z= 155 Ang, that is approximately 14,000 particles (Ref. 9, 10), and 95,000 particles for the five-atom TIP5P model of water (Ref. 10). 

A run of 15 minutes by the 6 cores/3.0 GHz machine is executed for t= 500 with the time step dt=0.01 (x 10^-14 s), while for 10 hours in the TIP5P simulation code. 
This file porw21.773a.pdf for the small potential gap V_top-V_bot= 1.5 kT with the dielectric constant eps=3 in the pore region shows four plots of potentials, particles of DNA and ions, those of all particles (every 5 water molecules), and the velocity distributions. One can see that the DNA chain moves toward the positive z direction into the cell volume. Moreover, the low dielectric constant eps(\r) in the pore makes the DNA blob more concentrated because counterions easily find negatively-charged DNA, which accelerates it upward to the cell region.

The water molecules are important by using the five-water water model @nanoWatPa.f03 because they are accelerated to the z-direction by electric potentials, shown here in porw31.773a.pdf. 

### (4) Coulomb-P3M (periodic) or Coulomb-Poisson (non periodic) simulations ###

The Coulomb forces are initiallized for water and hydrate (Ref. 9, "vitroid").
In this directory, the short-range and long-range Coulomb forces are used for a large system
in periodic p3m_perform routines with high accuracy in Sec.(2), and the non-periodic boundary 
conditions are utilized in Sec.(3). 

[1] Periodic system, the prototype

* @md3-para7.f03, with the parameter file paramE7.h (ca. 1,600 lines)

[2] Periodic system, the charge inversion

* @chginv3.f03, with the paramer file parm_inv13.h and 
the configure file CIMV13_config.START3 (ca. 5,100 lines)

[3] The non-periodic and periodic/bounded systems, the transport of DNA molecule

* @nanoAPGa.f03 (ca. 9,500 lines with graphics), with 
the parameter file paramAPGa.h and the configuration file PORV21_config.start3.

* The five-atom TIP5P water model using the bounded Poisson equation is utilized, @nanoWatPa.f03 
(11,000 lines), with paramWatP.h and PORW31_config.start3. In every time step it calculates
the Coulomb forces of HH-O-LL, L being the dummy charge sites, on top of the time-consuming Poisson equation.
It runs about 10 hours for t= 500 by 6 cores/3.0 GHz. 


Anyone may read, copy and rewrite these simulation codes under 
GNU General Public License v3.0, by keeping intact the first top 55 lines in Sec.[2] 
and 100 lines in Sec.[3] of the codes.


## References ##
1. M. Tanaka and A.Yu Grosberg, J.Chem.Phys., vol.115, 567-574 (2001).
2. M. Tanaka and A.Yu Grosberg, Euro.Phys.J., E7, 371-379 (2002).
3. M. Tanaka, Phys.Reviews., E68, 061501 (2003).
4. M. Tanaka, J. Physics: Condensed Matter, vol.16, S2127-2134 (2004).
5. Y. Rabin and M. Tanaka, Phys.Rev.Lett., vol.94, 148103 (2005).
6. M. Tanaka and M. Murakami, Comp.Phys.Commun., 241, pp. 56-63 (2019).
7. M. Deserno and C. Holm, J.Chem.Phys. 109, 7694–7701 (1998).
8. M. Tanaka, Collaboration at Institute of Polymerforchung, University of Mainz, Germany (1999).
9. M. Matsumoto, Water of GenIce, https://github.com/vitroid/.
10. The simulation codes of this directory are updated in 2025.


