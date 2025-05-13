 The Fortram module "hoecs.f90" can be used to calculate via 
 a DFT approach the nonlinear elastic constants of a solid 
 material. The method implemented in this module allow to 
 calculate elastic constants up to 8th order of a solid 
 material of arbitrary symmetry.  

 Compilation of this Fortran module requires a standard 
 fortran compiler and the BLAS and LAPACK libraries (see 
 Mafefile). The executable resulting from the compilation 
 works like a standard linux command, accepting several 
 arguments. Type "./hoecs.x -h" to see the list of possible 
 arguments. To become familiar with the operations this 
 module can perform, untar the compressed archive "examples.tgz", 
 and read the README files in the directories "./Si", "./Mg", 
 and "./SiO2". These directories contain example applications 
 to silicon, magnesium, and quartz, respectively. Input files, 
 scripts, and results are also available to replicate the 
 tasks described in the README files.

 Send emails to "angelo.bongiorno@csi.cuny.edu" for inquires, 
 to report bugs, or to request clarifications.
