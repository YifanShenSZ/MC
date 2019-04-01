!Monte Carlo simulation package for canonical ensemble
!This version considers no vibrations (flag: probably will add some day)
!Available ensembles: NVT, NPT, mVT, RVT, RPT
!Axis order: cartesian x-y-z, cylindrical z-x-y. Originate at the centre of the box.
!If 1 dimensional period, periodic boundary for z
!If 2 dimensional period, periodic boundary for x and y
!IO unit: temperature in K, pressure in bar, energy in kJ/mol, length in A, charge in e, density and concentration in mol/L
!Computation unit: atomic unit
!
!Tips for input:
!The origin of the input molecule xyz file should be the centre of mass
!For A+B=C reaction, B has to be a small molecule for detailed balance issue
!
!What to output: (these are just my personal interests, feel free to add more)
!NVT1D, density + pressure tensor
!NVT2D, density + pressure tensor
!NVT3D, 
!NPT  , concentration + chemical potential
!mVT1D, concentration + density + pressure tensor
!mVT2D, concentration + density + pressure tensor
!mVT3D, concentration
!RVT1D, concentration + density + pressure tensor
!RVT2D, concentration + density + pressure tensor
!RVT3D, concentration
!RPT3D, concentration + chemical potential
program main
    use General
    use MCBasis
    use ConfigurationalEnergy
    use MCAnalyzer
    use Evolve
    implicit none
    logical::AdvancedInput,ContinueAnalyzation
    integer::i,j,k,ii,jj
    !Density refers to inhomogeneous concentration
    real*8,allocatable,dimension(:)::Concentration,WidomChemicalPotential,Virial2D,TotalDensity
    real*8,allocatable,dimension(:,:)::Virial1D,MoleculeNumber
!---------- Initialize ----------
    write(*,'(A77)')'Monte Carlo simulation package for canonical ensemble, copy right Yifan Shen'
    call ReadInput(AdvancedInput)
    if(AdvancedInput) then
        open(unit=99,file='AdvancedInput',status='old')
            !reltol is fixed to be 1d-2 for NewStructure, 1d-4 for Optimize, 1d-6 for Analyze
            namelist /parameters/ abstol,MaxIteration,WallElement,WallBondLength,LayerDistance,&
                                      WallAtomDensity,epsilonwall,sigmawall,&!MCBasis
                                  bucketnumbers,&!MCAnalyzer
                                  Evolve_FollowFreq,Evolve_AdjustStep,Evolve_adjust!Evolve
            read(99,nml=parameters)
        close(99)
    end if
    call InitializeMCBasis(ContinueAnalyzation)
    call InitializeConfigurationalEnergy()
    if(jobtype=='Analyze') then
        call InitializeMCAnalyzer(ContinueAnalyzation)
    end if
    !For NewStructure, energy related MCBasis variables are uninitialized:
    !mtp, mtpold, intertp, intertpold, AverageEnergy, AverageGroundEnergy
    if(jobtype=='NewStructure') then
        select case(PeriodicDims)
            case(1)
                AverageEnergy=EnergyTotalAverageUpdated1D()
                !Write the wall structure
                !If you want to view the molecules and the wall together, manually add the wall atoms to the bottom of Structure.xyz
                open(unit=99,file='WallStructure.xyz',status='replace')
                    write(99,*)size(Wall,2)
                    write(99,*)
                    do i=1,size(Wall,2)
                        write(99,'(A2,F20.15,F20.15,F20.15)')WallElement,wall(1,i)/AInAU,wall(2,i)/AInAU,wall(3,i)/AInAU
                    end do
                close(99)
            case(2)
                AverageEnergy=EnergyTotalAverageUpdated2D()
                !Write the wall structure
                !If you want to view the molecules and the wall together, manually add the wall atoms to the bottom of Structure.xyz
                !2D wall structure is not needed to account for interactions
                call GenerateWall2D()
                wall=wall/AInAU
                open(unit=99,file='WallStructure.xyz',status='replace')
                    write(99,*)size(Wall,2)
                    write(99,*)
                    do i=1,size(Wall,2)
                        write(99,'(A2,F20.15,F20.15,F20.15)')WallElement,wall(1,i),wall(2,i),wall(3,i)
                    end do
                close(99)
                deallocate(wall)
            case(3)
                AverageEnergy=EnergyTotalAverageUpdated3D()
        end select
        AverageGroundEnergy=AverageEnergy
        mtpold=mtp
        intertpold=intertp
    end if
!------------- End --------------

!------------ Evolve ------------
    call showtime()
    write(*,*)'Initial average energy =',AverageEnergy/kJmolInAU,'kJ/mol'
    write(*,'(A12)')'Evolving...'
    !Run job
    select case(jobtype)
        case('NewStructure')
            select case(PeriodicDims)
                case(1)
                    call NVT1D()
                    !We traded accuracy for speed when initializing new strucutre
                    !For future use consistency we recompute the final energy at default accuracy level
                    reltol=1d-4
                    call EwaldSummationDumpFactor1D()
                    AverageEnergy=EnergyTotalAverageUpdated1D()
                case(2)
                    call NVT2D()
                    !We traded accuracy for speed when initializing new strucutre
                    !For future use consistency we recompute the final energy at default accuracy level
                    reltol=1d-4
                    call EwaldSummationDumpFactor2D()
                    AverageEnergy=EnergyTotalAverageUpdated2D()
                case(3)
                    call NVT3D()
                    !We traded accuracy for speed when initializing new strucutre
                    !For future use consistency we recompute the final energy at default accuracy level
                    reltol=1d-4
                    call EwaldSummationDumpFactor3D()
                    AverageEnergy=EnergyTotalAverageUpdated3D()
            end select
            AverageGroundEnergy=AverageEnergy
        case('Optimize')
            select case(ensemble)
                case('NVT')
                    select case(PeriodicDims)
                        case(1)
                            call NVT1D()
                        case(2)
                            call NVT2D()
                        case(3)
                            call NVT3D()
                    end select
                case('NPT')
                    if(PeriodicDims<3) then
                        stop 'It is impossible for a non-bulk system to hold a constant pressure'
                    else
                        call NPT()
                    end if
                case('mVT')
                    select case(PeriodicDims)
                        case(1)
                            call mVT1D()
                        case(2)
                            call mVT2D()
                        case(3)
                            call mVT3D()
                    end select
                case('mPT')
                    !From Gibbs phase rule: f = 2 + c - p, where f = # degree of freedom, c = # components, p = # phases
                    !f intensive properties have determined all intensive properties of the system,
                    !you just need 1 additional extensive degree of freedom to determine the system
                    !Although extensive properties computed from the simulation system is meaningless,
                    !they still must be determined in order to avoid theoretical issue
                    stop 'Specifying all intensive properties overdetermines the system'
                case('RVT')
                    select case(PeriodicDims)
                        case(1)
                            call RVT1D()
                        case(2)
                            call RVT2D()
                        case(3)
                            call RVT3D()
                    end select
                case('RPT')
                    if(PeriodicDims<3) then
                        stop 'It is impossible for a non-bulk system to hold a constant pressure'
                    else
                        call RPT()
                    end if
                case default
                    stop 'Unsupported ensemble'
            end select
        case('Analyze')
            select case(ensemble)
                case('NVT')
                    select case(PeriodicDims)
                        case(1)
                            call NVT1DA(Virial1D,MoleculeNumber)
                        case(2)
                            call NVT2DA(Virial2D,MoleculeNumber)
                        case(3)
                            stop 'Bulk NVT ensemble is trivial'
                    end select
                case('NPT')
                    if(PeriodicDims<3) then
                        stop 'It is impossible for a non-bulk system to hold a constant pressure'
                    else
                        call NPTA(Concentration,WidomChemicalPotential)
                    end if
                case('mVT')
                    select case(PeriodicDims)
                        case(1)
                            call mVT1DA(Concentration,Virial1D,MoleculeNumber)
                        case(2)
                            call mVT2DA(Concentration,Virial2D,MoleculeNumber)
                        case(3)
                            call mVT3DA(Concentration)
                    end select
                case('mPT')
                    !See explanation in Optimize section
                    stop 'Specifying all intensive properties overdetermines the system'
                case('RVT')
                    select case(PeriodicDims)
                        case(1)
                            call RVT1DA(Concentration,Virial1D,MoleculeNumber)
                        case(2)
                            call RVT2DA(Concentration,Virial2D,MoleculeNumber)
                        case(3)
                            call RVT3DA(Concentration)
                    end select
                case('RPT')
                    if(PeriodicDims<3) then
                        stop 'It is impossible for a non-bulk system to hold a constant pressure'
                    else
                        call RPTA(Concentration,WidomChemicalPotential)
                    end if
                case default
                    stop 'Unsupported ensemble'
            end select
        case default
            stop 'Program abort: unsupported job type'
    end select
!------------- End --------------

!------------ Output ------------
    write(*,*)'Writing final output...'
    call WriteCheckPoint()
    !Output analyzation
    if(jobtype=='Analyze') then
        allocate(TotalDensity(bucketnumbers))!Work space for density + pressure tensor
        select case(ensemble)
            case('NVT')
                select case(PeriodicDims)
                    case(1)
                        call MoleculeNumber2Density1D(MoleculeNumber,TotalDensity)
                        call Virial2Pressure1D(Virial1D,TotalDensity)
                    case(2)
                        call MoleculeNumber2Density2D(MoleculeNumber,TotalDensity)
                        call Virial2Pressure2D(Virial2D,TotalDensity)
                end select
            case('NPT')
                Concentration=Concentration*(AInAU*AInAU*AInAU*1d27)/NAvogadro!Convert to human unit
                open(unit=99,file='Concentration.txt',status='replace')
                    write(99,*)Concentration
                close(99)
                !To save cpu time evolve subroutine returns Widom insertion parameter B rather than chemical potential
                forall(i=1:MoleculeKinds)
                    WidomChemicalPotential(i)=-temperature*Log(WidomChemicalPotential(i)/Concentration(i)&
                        *(MoleculeInput(i).mass*temperature/pim2)**1.5d0)/kJmolInAU!Convert to human unit
                end forall
                open(unit=99,file='ChemicalPotential.txt')
                    write(99,*)WidomChemicalPotential
                close(99)
            !V type ensembles hold volume constant, so to save cpu time,
            !evolve subroutines return ensemble average molecule number rather than concentration
            case('mVT')
                select case(PeriodicDims)
                    case(1)
                        Concentration=Concentration/volume*(AInAU*AInAU*AInAU*1d27)/NAvogadro
                        open(unit=99,file='Concentration.txt',status='replace')
                            write(99,*)Concentration
                        close(99)
                        call MoleculeNumber2Density1D(MoleculeNumber,TotalDensity)
                        call Virial2Pressure1D(Virial1D,TotalDensity)
                    case(2)
                        Concentration=Concentration/volume*(AInAU*AInAU*AInAU*1d27)/NAvogadro
                        open(unit=99,file='Concentration.txt',status='replace')
                            write(99,*)Concentration
                        close(99)
                        call MoleculeNumber2Density2D(MoleculeNumber,TotalDensity)
                        call Virial2Pressure2D(Virial2D,TotalDensity)
                    case(3)
                        Concentration=Concentration/volume*(AInAU*AInAU*AInAU*1d27)/NAvogadro
                        open(unit=99,file='Concentration.txt',status='replace')
                            write(99,*)Concentration
                        close(99)
                end select
            !V type ensembles hold volume constant, so to save cpu time,
            !evolve subroutines return ensemble average molecule number rather than concentration
            case('RVT')
                select case(PeriodicDims)
                    case(1)
                        Concentration=Concentration/volume*(AInAU*AInAU*AInAU*1d27)/NAvogadro
                        open(unit=99,file='Concentration.txt',status='replace')
                            write(99,*)Concentration
                        close(99)
                        call MoleculeNumber2Density1D(MoleculeNumber,TotalDensity)
                        call Virial2Pressure1D(Virial1D,TotalDensity)
                    case(2)
                        Concentration=Concentration/volume*(AInAU*AInAU*AInAU*1d27)/NAvogadro
                        open(unit=99,file='Concentration.txt',status='replace')
                            write(99,*)Concentration
                        close(99)
                        call MoleculeNumber2Density2D(MoleculeNumber,TotalDensity)
                        call Virial2Pressure2D(Virial2D,TotalDensity)
                    case(3)
                        Concentration=Concentration/volume*(AInAU*AInAU*AInAU*1d27)/NAvogadro
                        open(unit=99,file='Concentration.txt',status='replace')
                            write(99,*)Concentration
                        close(99)
                end select
            case('RPT')
                Concentration=Concentration*(AInAU*AInAU*AInAU*1d27)/NAvogadro!Convert to human unit
                open(unit=99,file='Concentration.txt',status='replace')
                    write(99,*)Concentration
                close(99)
                !To save cpu time evolve subroutine returns Widom insertion parameter B rather than chemical potential
                forall(i=1:MoleculeKinds)
                    WidomChemicalPotential(i)=-temperature*Log(WidomChemicalPotential(i)/Concentration(i)&
                        *(MoleculeInput(i).mass*temperature/pim2)**1.5d0)/kJmolInAU!Convert to human unit
                end forall
                open(unit=99,file='ChemicalPotential.txt')
                    write(99,*)WidomChemicalPotential
                close(99)
        end select
    end if
    write(*,*)'Mission complete'
!------------- End --------------
end program main