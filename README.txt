#This program has been written by Rémi Khatib.
#It calculates dipole moments and polarizability
#from atomic positions (.xyz file)
#It is based on on the DID model and served for
#the article:
#Khatib et al., J. Phys. Chem. C, 120, 18665, 2016

#March 2020:
#
#Code modified by Hossam Elgabarty in a number of ways:
#
#1. Treat a system with arbitrary topology, one simply reads
#   a psf topology file which defines "molecules" (psf residues), additionally, several input blocks,
#   one per molecule kind, are also needed to define the gas-phase dipole moment,
#   polarizability, hyperpolarizability, together with the molecular coordinate frame
#   that corresponds to these tensors.
#
#2. OMP parallelization (In progress)
#
# ---------   TO BE DONE  ------------------
#2. Implement Torii's extended dipole model:
#   Chem. Phys. Lett., 353, 431, 2002 (See the beta-mu term in eq. 5)
#
#3. Periodic boundary conditions for dipolar interactions:
#   - A. Aguado and P. A. Madden, J. Chem. Phys. 119, 7471, 2003
#   - T. Laino and J. Hutter, J. Chem. Phys. 129, 074102, 2008
#
#4. Solving the SCF on a grid by a collocation of Legendre functions.
#
#
#To handle psf topologies, the code now requires linking against the loos library,
#see http://grossfieldlab.github.io/loos/


#=========
#Make file
#=========
#Compilation
> make

#Key words in the makefile to read the Yuki's input and produce
#Yuki's output. If commented, you will have the readable version
YUKI='-D YUKI'

#Warnings vs production compilation
#In the makefile, comment or uncomment the flag
DEBUG='warning'

#Remove the object files
> make clean

#Remove even more files
> make very_clean

#Make a tag for emacs (etag)
> make tag


#=====
#Tests
#=====
#The input files used for the article
#Khatib et al., J. Phys. Chem. C, 120, 18665, 2016
#are listed and commented in the __tests directory
#
#There are 2 sets of input. Indeed, there is one with
#readable output (__tests/readable) and another one with
#a kind of output used by Yuki Nagata (__tests/yuki).
#There is linear relation ship between the 2 outputs
#of the dipole / polarizability (see the code for more details).
#The readable/yuki input differs by 1 line only but
#if the makefile is not the good one, you may have nothing.
#Moreover, one extra file is produced with the readble version in
#order to print the position of the center of mass
#Please, refer to the makefile section (YUKI='-D YUKI')
#
#For each set of input, there are 4 kind of input files One for each model.
#Please, refer to the article to have a precise
#description of the models
#
#In the __test directory can be found 2 extra files
# -1to3_points_pol.ods (for libreoffice calc) where I
#  describe how o go from a single point description
#  to a triple point description
# -pos.xyz which is here just for the tests, since it
#  10 steps of an SPC/E dynamic
#
#
#Tests for the new features by Hossam Elgabarty:
#...
#...
#...
#...
