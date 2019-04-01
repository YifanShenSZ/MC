##############################################
#                                            #
#        Makefile for compiling SYFMC        #
#                                            #
##############################################

compiler = ifort
f90s = General.f90 MCBasis.f90 Nonbonded.f90 ConfigurationalEnergy.f90 MCAnalyzer.f90 Trial.f90 Evolve.f90 Main.f90
exe = SYFMC.exe
flags = -fast -march=core-avx2

$(exe): $(f90s)
    $(compiler) $(flags) $^ -o $(exe)

clean:
    rm $(exe)
    rm *.mod
