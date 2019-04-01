##############################################
#                                            #
#        Makefile for compiling SYFMC        #
#                                            #
##############################################

compiler = ifort
f90s = General.f90 MCBasis.f90 Nonbonded.f90 ConfigurationalEnergy.f90 MCAnalyzer.f90 Trial.f90 Evolve.f90 Main.f90
exe = SYFMC.exe
# For sse2: Alex2 seems not to support -static option, so cannot use -fast
flags = -m64 -msse2 -ipo -O3 -no-prec-div -fp-model fast=2
# For avx2:
#flags = -fast -march=core-avx2

$(exe): $(f90s)
    $(compiler) $(flags) $^ -o $(exe)

clean:
    rm $(exe)
    rm *.mod