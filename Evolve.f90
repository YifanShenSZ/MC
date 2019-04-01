!Translational trial move has anround 50% acceptance ratio
!Rotational trial move is uniform
!Volume change probability is 1%
!Particle insertion probability is 30% during optimization while 1% during analyzation
!Reactive trial probability is 30% during optimization while 1% during analyzation
module Evolve
    use General
    use MCBasis
    use ConfigurationalEnergy
    use MCAnalyzer
    use Trial
    implicit none

!Evolve module only parameter
    !Every how many steps write check point, try adjust translational trial step,
    !adjust only when such trial has been made for enough times 
    integer::Evolve_FollowFreq=1000000,Evolve_AdjustStep=1000,Evolve_adjust=100

contains
!----------- Optimization evolution -----------
    !NVT
    !For 1D
    subroutine NVT1D()
        logical::accept
        integer::i,j,k
        real*8::direction
        !Trial for TotalSteps times
        do k=1,TotalSteps
            !Trial move: relaxation
            call Relax1D(i,j,direction,accept)
            !Follow evolving progress
            if(mod(k,Evolve_FollowFreq)==0) then
                call showtime()
                write(*,*)'step =',k
                AverageEnergy=EnergyTotalAverageCurrentOutput1D()
                if(AverageEnergy<AverageGroundEnergy) then
                    AverageGroundEnergy=AverageEnergy
                end if
                write(*,*)'Save check point files'
                call WriteCheckPoint()
            end if
            !Adjust translational trial step
            if(mod(k,Evolve_AdjustStep)==0) then
                do i=1,MoleculeKinds
                    do j=1,3
                        if(TrialCount(j,i)>Evolve_adjust) then
                            if(TrialSuccessCount(j,i)<TrialCount(j,i)/2) then
                                TrialStep(j,i)=TrialStep(j,i)*0.95d0
                            else
                                TrialStep(j,i)=TrialStep(j,i)*1.05d0
                                if(TrialStep(j,i)>MaxTranStep(j)) then
                                    TrialStep(j,i)=MaxTranStep(j)
                                end if
                            end if
                            TrialSuccessCount(j,i)=0
                            TrialCount(j,i)=0
                        end if
                    end do
                end do
            end if
        end do
    end subroutine NVT1D
    !For 2D
    subroutine NVT2D()
        logical::accept
        integer::i,j,k
        real*8::direction
        !Trial for TotalSteps times
        do k=1,TotalSteps
            !Trial move: relaxation
            call Relax2D(i,j,direction,accept)
            !Follow evolving progress
            if(mod(k,Evolve_FollowFreq)==0) then
                call showtime()
                write(*,*)'step =',k
                AverageEnergy=EnergyTotalAverageCurrentOutput2D()
                if(AverageEnergy<AverageGroundEnergy) then
                    AverageGroundEnergy=AverageEnergy
                end if
                write(*,*)'Save check point files'
                call WriteCheckPoint()
            end if
            !Adjust translational trial step
            if(mod(k,Evolve_AdjustStep)==0) then
                do i=1,MoleculeKinds
                    do j=1,3
                        if(TrialCount(j,i)>Evolve_adjust) then
                            if(TrialSuccessCount(j,i)<TrialCount(j,i)/2) then
                                TrialStep(j,i)=TrialStep(j,i)*0.95d0
                            else
                                TrialStep(j,i)=TrialStep(j,i)*1.05d0
                                if(TrialStep(j,i)>MaxTranStep(j)) then
                                    TrialStep(j,i)=MaxTranStep(j)
                                end if
                            end if
                            TrialSuccessCount(j,i)=0
                            TrialCount(j,i)=0
                        end if
                    end do
                end do
            end if
        end do
    end subroutine NVT2D
    !For 3D
    subroutine NVT3D()
        logical::accept
        integer::i,j,k
        real*8::direction
        !Trial for TotalSteps times
        do k=1,TotalSteps
            !Trial move: relaxation
            call Relax3D(i,j,direction,accept)
            !Follow evolving progress
            if(mod(k,Evolve_FollowFreq)==0) then
                call showtime()
                write(*,*)'step =',k
                AverageEnergy=EnergyTotalAverageCurrentOutput3D()
                if(AverageEnergy<AverageGroundEnergy) then
                    AverageGroundEnergy=AverageEnergy
                end if
                write(*,*)'Save check point files'
                call WriteCheckPoint()
            end if
            !Adjust translational trial step
            if(mod(k,Evolve_AdjustStep)==0) then
                do i=1,MoleculeKinds
                    do j=1,3
                        if(TrialCount(j,i)>Evolve_adjust) then
                            if(TrialSuccessCount(j,i)<TrialCount(j,i)/2) then
                                TrialStep(j,i)=TrialStep(j,i)*0.95d0
                            else
                                TrialStep(j,i)=TrialStep(j,i)*1.05d0
                                if(TrialStep(j,i)>MaxTranStep(j)) then
                                    TrialStep(j,i)=MaxTranStep(j)
                                end if
                            end if
                            TrialSuccessCount(j,i)=0
                            TrialCount(j,i)=0
                        end if
                    end do
                end do
            end if
        end do
    end subroutine NVT3D
    
    !NPT
    subroutine NPT()
        logical::accept
        integer::i,j,k
        real*8::direction
        !Trial for TotalSteps times
        do k=1,TotalSteps
            !Trial move: volume change, relaxation
            call random_number(direction)
            if(direction<0.01d0) then
                call VolumeChange(accept)
            else
                call Relax3D(i,j,direction,accept)
            end if
            !Follow evolving progress
            if(mod(k,Evolve_FollowFreq)==0) then
                call showtime()
                write(*,*)'step =',k
                AverageEnergy=EnergyTotalAverageCurrentOutput3D()
                if(AverageEnergy<AverageGroundEnergy) then
                    AverageGroundEnergy=AverageEnergy
                end if
                write(*,*)'Save check point files'
                call WriteCheckPoint()
            end if
            !Adjust translational trial step
            if(mod(k,Evolve_AdjustStep)==0) then
                do i=1,MoleculeKinds
                    do j=1,3
                        if(TrialCount(j,i)>Evolve_adjust) then
                            if(TrialSuccessCount(j,i)<TrialCount(j,i)/2) then
                                TrialStep(j,i)=TrialStep(j,i)*0.95d0
                            else
                                TrialStep(j,i)=TrialStep(j,i)*1.05d0
                                if(TrialStep(j,i)>MaxTranStep(j)) then
                                    TrialStep(j,i)=MaxTranStep(j)
                                end if
                            end if
                            TrialSuccessCount(j,i)=0
                            TrialCount(j,i)=0
                        end if
                    end do
                end do
            end if
        end do
    end subroutine NPT
    
    !mVT
    !For 1D
    subroutine mVT1D()
        logical::accept
        integer::i,j,k
        real*8::direction
        !Trial for TotalSteps times
        do k=1,TotalSteps
            !Trial move: insertion or deletion, relaxation
            call random_number(direction)
            if(direction<0.3d0) then
                call Insertion1D(i,j,direction,accept)
            else
                call Relax1D(i,j,direction,accept)
            end if
            !Follow evolving progress
            if(mod(k,Evolve_FollowFreq)==0) then
                call showtime()
                write(*,*)'step =',k
                AverageEnergy=EnergyTotalAverageCurrentOutput1D()
                if(AverageEnergy<AverageGroundEnergy) then
                    AverageGroundEnergy=AverageEnergy
                end if
                write(*,*)'Save check point files'
                call WriteCheckPoint()
            end if
            !Adjust translational trial step
            if(mod(k,Evolve_AdjustStep)==0) then
                do i=1,MoleculeKinds
                    do j=1,3
                        if(TrialCount(j,i)>Evolve_adjust) then
                            if(TrialSuccessCount(j,i)<TrialCount(j,i)/2) then
                                TrialStep(j,i)=TrialStep(j,i)*0.95d0
                            else
                                TrialStep(j,i)=TrialStep(j,i)*1.05d0
                                if(TrialStep(j,i)>MaxTranStep(j)) then
                                    TrialStep(j,i)=MaxTranStep(j)
                                end if
                            end if
                            TrialSuccessCount(j,i)=0
                            TrialCount(j,i)=0
                        end if
                    end do
                end do
            end if
        end do
    end subroutine mVT1D
    !For 2D
    subroutine mVT2D()
        logical::accept
        integer::i,j,k
        real*8::direction
        !Trial for TotalSteps times
        do k=1,TotalSteps
            !Trial move: insertion or deletion, relaxation
            call random_number(direction)
            if(direction<0.3d0) then
                call Insertion2D(i,j,direction,accept)
            else
                call Relax2D(i,j,direction,accept)
            end if
            !Follow evolving progress
            if(mod(k,Evolve_FollowFreq)==0) then
                call showtime()
                write(*,*)'step =',k
                AverageEnergy=EnergyTotalAverageCurrentOutput2D()
                if(AverageEnergy<AverageGroundEnergy) then
                    AverageGroundEnergy=AverageEnergy
                end if
                write(*,*)'Save check point files'
                call WriteCheckPoint()
            end if
            !Adjust translational trial step
            if(mod(k,Evolve_AdjustStep)==0) then
                do i=1,MoleculeKinds
                    do j=1,3
                        if(TrialCount(j,i)>Evolve_adjust) then
                            if(TrialSuccessCount(j,i)<TrialCount(j,i)/2) then
                                TrialStep(j,i)=TrialStep(j,i)*0.95d0
                            else
                                TrialStep(j,i)=TrialStep(j,i)*1.05d0
                                if(TrialStep(j,i)>MaxTranStep(j)) then
                                    TrialStep(j,i)=MaxTranStep(j)
                                end if
                            end if
                            TrialSuccessCount(j,i)=0
                            TrialCount(j,i)=0
                        end if
                    end do
                end do
            end if
        end do
    end subroutine mVT2D
    !For 3D
    subroutine mVT3D()
        logical::accept
        integer::i,j,k
        real*8::direction
        !Trial for TotalSteps times
        do k=1,TotalSteps
            !Trial move: insertion or deletion, relaxation
            call random_number(direction)
            if(direction<0.3d0) then
                call Insertion3D(i,j,direction,accept)
            else
                call Relax3D(i,j,direction,accept)
            end if
            !Follow evolving progress
            if(mod(k,Evolve_FollowFreq)==0) then
                call showtime()
                write(*,*)'step =',k
                AverageEnergy=EnergyTotalAverageCurrentOutput3D()
                if(AverageEnergy<AverageGroundEnergy) then
                    AverageGroundEnergy=AverageEnergy
                end if
                write(*,*)'Save check point files'
                call WriteCheckPoint()
            end if
            !Adjust translational trial step
            if(mod(k,Evolve_AdjustStep)==0) then
                do i=1,MoleculeKinds
                    do j=1,3
                        if(TrialCount(j,i)>Evolve_adjust) then
                            if(TrialSuccessCount(j,i)<TrialCount(j,i)/2) then
                                TrialStep(j,i)=TrialStep(j,i)*0.95d0
                            else
                                TrialStep(j,i)=TrialStep(j,i)*1.05d0
                                if(TrialStep(j,i)>MaxTranStep(j)) then
                                    TrialStep(j,i)=MaxTranStep(j)
                                end if
                            end if
                            TrialSuccessCount(j,i)=0
                            TrialCount(j,i)=0
                        end if
                    end do
                end do
            end if
        end do
    end subroutine mVT3D
    
    !RVT
    !For 1D
    subroutine RVT1D()
        logical::accept
        integer::i,j,k,chosenp
        integer,dimension(2)::r,chosenr
        real*8::direction
        !Trial for TotalSteps times
        do k=1,TotalSteps
            !Trial move: reaction, relaxation
            call random_number(direction)
            if(direction<0.3d0) then
                !Select a reaction
                call random_number(direction)
                i=ceiling(direction*ReactionKinds)
                if(ReactionInput(i).NReactant==1) then
                    call ReactAB1D(i,r(1),r(2),chosenr(1),chosenr(2),direction,accept)
                else
                    call ReactABC1D(i,r,j,chosenr,chosenp,direction,accept)
                end if
            else
                call Relax1D(i,j,direction,accept)
            end if
            !Follow evolving progress
            if(mod(k,Evolve_FollowFreq)==0) then
                call showtime()
                write(*,*)'step =',k
                AverageEnergy=EnergyTotalAverageCurrentOutput1D()
                if(AverageEnergy<AverageGroundEnergy) then
                    AverageGroundEnergy=AverageEnergy
                end if
                write(*,*)'Save check point files'
                call WriteCheckPoint()
            end if
            !Adjust translational trial step
            if(mod(k,Evolve_AdjustStep)==0) then
                do i=1,MoleculeKinds
                    do j=1,3
                        if(TrialCount(j,i)>Evolve_adjust) then
                            if(TrialSuccessCount(j,i)<TrialCount(j,i)/2) then
                                TrialStep(j,i)=TrialStep(j,i)*0.95d0
                            else
                                TrialStep(j,i)=TrialStep(j,i)*1.05d0
                                if(TrialStep(j,i)>MaxTranStep(j)) then
                                    TrialStep(j,i)=MaxTranStep(j)
                                end if
                            end if
                            TrialSuccessCount(j,i)=0
                            TrialCount(j,i)=0
                        end if
                    end do
                end do
            end if
        end do
    end subroutine RVT1D
    !For 2D
    subroutine RVT2D()
        logical::accept
        integer::i,j,k,chosenp
        integer,dimension(2)::r,chosenr
        real*8::direction
        !Trial for TotalSteps times
        do k=1,TotalSteps
            !Trial move: reaction, relaxation
            call random_number(direction)
            if(direction<0.3d0) then
                !Select a reaction
                call random_number(direction)
                i=ceiling(direction*ReactionKinds)
                if(ReactionInput(i).NReactant==1) then
                    call ReactAB2D(i,r(1),r(2),chosenr(1),chosenr(2),direction,accept)
                else
                    call ReactABC2D(i,r,j,chosenr,chosenp,direction,accept)
                end if
            else
                call Relax2D(i,j,direction,accept)
            end if
            !Follow evolving progress
            if(mod(k,Evolve_FollowFreq)==0) then
                call showtime()
                write(*,*)'step =',k
                AverageEnergy=EnergyTotalAverageCurrentOutput2D()
                if(AverageEnergy<AverageGroundEnergy) then
                    AverageGroundEnergy=AverageEnergy
                end if
                write(*,*)'Save check point files'
                call WriteCheckPoint()
            end if
            !Adjust translational trial step
            if(mod(k,Evolve_AdjustStep)==0) then
                do i=1,MoleculeKinds
                    do j=1,3
                        if(TrialCount(j,i)>Evolve_adjust) then
                            if(TrialSuccessCount(j,i)<TrialCount(j,i)/2) then
                                TrialStep(j,i)=TrialStep(j,i)*0.95d0
                            else
                                TrialStep(j,i)=TrialStep(j,i)*1.05d0
                                if(TrialStep(j,i)>MaxTranStep(j)) then
                                    TrialStep(j,i)=MaxTranStep(j)
                                end if
                            end if
                            TrialSuccessCount(j,i)=0
                            TrialCount(j,i)=0
                        end if
                    end do
                end do
            end if
        end do
    end subroutine RVT2D
    !For 3D
    subroutine RVT3D()
        logical::accept
        integer::i,j,k,chosenp
        integer,dimension(2)::r,chosenr
        real*8::direction
        !Trial for TotalSteps times
        do k=1,TotalSteps
            !Trial move: reaction, relaxation
            call random_number(direction)
            if(direction<0.3d0) then
                !Select a reaction
                call random_number(direction)
                i=ceiling(direction*ReactionKinds)
                if(ReactionInput(i).NReactant==1) then
                    call ReactAB3D(i,r(1),r(2),chosenr(1),chosenr(2),direction,accept)
                else
                    call ReactABC3D(i,r,j,chosenr,chosenp,direction,accept)
                end if
            else
                call Relax3D(i,j,direction,accept)
            end if
            !Follow evolving progress
            if(mod(k,Evolve_FollowFreq)==0) then
                call showtime()
                write(*,*)'step =',k
                AverageEnergy=EnergyTotalAverageCurrentOutput3D()
                if(AverageEnergy<AverageGroundEnergy) then
                    AverageGroundEnergy=AverageEnergy
                end if
                write(*,*)'Save check point files'
                call WriteCheckPoint()
            end if
            !Adjust translational trial step
            if(mod(k,Evolve_AdjustStep)==0) then
                do i=1,MoleculeKinds
                    do j=1,3
                        if(TrialCount(j,i)>Evolve_adjust) then
                            if(TrialSuccessCount(j,i)<TrialCount(j,i)/2) then
                                TrialStep(j,i)=TrialStep(j,i)*0.95d0
                            else
                                TrialStep(j,i)=TrialStep(j,i)*1.05d0
                                if(TrialStep(j,i)>MaxTranStep(j)) then
                                    TrialStep(j,i)=MaxTranStep(j)
                                end if
                            end if
                            TrialSuccessCount(j,i)=0
                            TrialCount(j,i)=0
                        end if
                    end do
                end do
            end if
        end do
    end subroutine RVT3D

    !RPT
    subroutine RPT()
        logical::accept
        integer::i,j,k,chosenp
        integer,dimension(2)::r,chosenr
        real*8::direction
        !Trial for TotalSteps times
        do k=1,TotalSteps
            !Trial move: volume change, reaction, relaxation
            call random_number(direction)
            if(direction<0.01d0) then
                call VolumeChange(accept)
            else if(direction<0.31d0) then
                !Select a reaction
                call random_number(direction)
                i=ceiling(direction*ReactionKinds)
                if(ReactionInput(i).NReactant==1) then
                    call ReactAB3D(i,r(1),r(2),chosenr(1),chosenr(2),direction,accept)
                else
                    call ReactABC3D(i,r,j,chosenr,chosenp,direction,accept)
                end if
            else
                call Relax3D(i,j,direction,accept)
            end if
            !Follow evolving progress
            if(mod(k,Evolve_FollowFreq)==0) then
                call showtime()
                write(*,*)'step =',k
                AverageEnergy=EnergyTotalAverageCurrentOutput3D()
                if(AverageEnergy<AverageGroundEnergy) then
                    AverageGroundEnergy=AverageEnergy
                end if
                write(*,*)'Save check point files'
                call WriteCheckPoint()
            end if
            !Adjust translational trial step
            if(mod(k,Evolve_AdjustStep)==0) then
                do i=1,MoleculeKinds
                    do j=1,3
                        if(TrialCount(j,i)>Evolve_adjust) then
                            if(TrialSuccessCount(j,i)<TrialCount(j,i)/2) then
                                TrialStep(j,i)=TrialStep(j,i)*0.95d0
                            else
                                TrialStep(j,i)=TrialStep(j,i)*1.05d0
                                if(TrialStep(j,i)>MaxTranStep(j)) then
                                    TrialStep(j,i)=MaxTranStep(j)
                                end if
                            end if
                            TrialSuccessCount(j,i)=0
                            TrialCount(j,i)=0
                        end if
                    end do
                end do
            end if
        end do
    end subroutine RPT
!-------------------- End ---------------------
    
!------------ Analyzation evolution -----------
    !NVT
    !For 1D
    subroutine NVT1DA(Virial,MoleculeNumber)
        real*8,allocatable,dimension(:,:),intent(out)::Virial,MoleculeNumber
        logical::accept
        integer::i,j,k
        real*8::direction
        real*8,allocatable,dimension(:)::TotalDensity
        real*8,allocatable,dimension(:,:)::VirialTemp,MoleculeNumberTemp
        !Allocate memory to store your result, and set the counters to zero
            allocate(Virial(2,bucketnumbers))
            Virial=0d0
            allocate(VirialTemp(2,bucketnumbers))
            VirialTemp=0d0
            allocate(MoleculeNumber(bucketnumbers,MoleculeKinds))
            MoleculeNumber=0d0
            allocate(MoleculeNumberTemp(bucketnumbers,MoleculeKinds))
            MoleculeNumberTemp=0d0
            allocate(TotalDensity(bucketnumbers))
        !Trial for TotalSteps times
        do k=1,TotalSteps
            !Trial move: relaxation
            call Relax1D(i,j,direction,accept)
            if(accept) then
                call UpdateChosenVirial1D(i,j)
            end if
            !Collect data
            VirialTemp=VirialTemp+Virial1DInEachBin()
            MoleculeNumberTemp=MoleculeNumberTemp+MoleculeNumberInEachBin()
            !Follow evolving progress and add the analyzation
            if(mod(k,Evolve_FollowFreq)==0) then
                call showtime()
                write(*,*)'step =',k
                AverageEnergy=EnergyTotalAverageCurrentOutput1D()
                if(AverageEnergy<AverageGroundEnergy) then
                    AverageGroundEnergy=AverageEnergy
                end if
                write(*,*)'Save check point files'
                call WriteCheckPoint()
                call WriteVirialCheckPoint()
                !Add the analyzation
                Virial=Virial+VirialTemp/TotalSteps
                MoleculeNumber=MoleculeNumber+MoleculeNumberTemp/TotalSteps
                !Output current ensemble average
                VirialTemp=Virial*TotalSteps/k
                MoleculeNumberTemp=MoleculeNumber*TotalSteps/k
                call MoleculeNumber2Density1D(MoleculeNumberTemp,TotalDensity)
                call Virial2Pressure1D(VirialTemp,TotalDensity)
                !Reset the counters
                VirialTemp=0d0
                MoleculeNumberTemp=0d0
            end if
            !Adjust translational trial step
            if(mod(k,Evolve_AdjustStep)==0) then
                do i=1,MoleculeKinds
                    do j=1,3
                        if(TrialCount(j,i)>Evolve_adjust) then
                            if(TrialSuccessCount(j,i)<TrialCount(j,i)/2) then
                                TrialStep(j,i)=TrialStep(j,i)*0.95d0
                            else
                                TrialStep(j,i)=TrialStep(j,i)*1.05d0
                                if(TrialStep(j,i)>MaxTranStep(j)) then
                                    TrialStep(j,i)=MaxTranStep(j)
                                end if
                            end if
                            TrialSuccessCount(j,i)=0
                            TrialCount(j,i)=0
                        end if
                    end do
                end do
            end if
        end do
    end subroutine NVT1DA
    !For 2D
    subroutine NVT2DA(Virial,MoleculeNumber)
        real*8,allocatable,dimension(:),intent(out)::Virial
        real*8,allocatable,dimension(:,:),intent(out)::MoleculeNumber
        logical::accept
        integer::i,j,k
        real*8::direction
        real*8,allocatable,dimension(:)::VirialTemp,TotalDensity
        real*8,allocatable,dimension(:,:)::MoleculeNumberTemp
        !Allocate memory to store your result, and set the counters to zero
            allocate(Virial(bucketnumbers))
            Virial=0d0
            allocate(VirialTemp(bucketnumbers))
            VirialTemp=0d0
            allocate(MoleculeNumber(bucketnumbers,MoleculeKinds))
            MoleculeNumber=0d0
            allocate(MoleculeNumberTemp(bucketnumbers,MoleculeKinds))
            MoleculeNumberTemp=0d0
            allocate(TotalDensity(bucketnumbers))
        !Trial for TotalSteps times
        do k=1,TotalSteps
            !Trial move: relaxation
            call Relax2D(i,j,direction,accept)
            if(accept) then
                call UpdateChosenVirial2D(i,j)
            end if
            !Collect data
            VirialTemp=VirialTemp+Virial2DInEachBin()
            MoleculeNumberTemp=MoleculeNumberTemp+MoleculeNumberInEachBin()
            !Follow evolving progress and add the analyzation
            if(mod(k,Evolve_FollowFreq)==0) then
                call showtime()
                write(*,*)'step =',k
                AverageEnergy=EnergyTotalAverageCurrentOutput2D()
                if(AverageEnergy<AverageGroundEnergy) then
                    AverageGroundEnergy=AverageEnergy
                end if
                write(*,*)'Save check point files'
                call WriteCheckPoint()
                call WriteVirialCheckPoint()
                !Add the analyzation
                Virial=Virial+VirialTemp/TotalSteps
                MoleculeNumber=MoleculeNumber+MoleculeNumberTemp/TotalSteps
                !Output current ensemble average
                VirialTemp=Virial*TotalSteps/k
                MoleculeNumberTemp=MoleculeNumber*TotalSteps/k
                call MoleculeNumber2Density2D(MoleculeNumberTemp,TotalDensity)
                call Virial2Pressure2D(VirialTemp,TotalDensity)
                !Reset the counters
                VirialTemp=0d0
                MoleculeNumberTemp=0d0
            end if
            !Adjust translational trial step
            if(mod(k,Evolve_AdjustStep)==0) then
                do i=1,MoleculeKinds
                    do j=1,3
                        if(TrialCount(j,i)>Evolve_adjust) then
                            if(TrialSuccessCount(j,i)<TrialCount(j,i)/2) then
                                TrialStep(j,i)=TrialStep(j,i)*0.95d0
                            else
                                TrialStep(j,i)=TrialStep(j,i)*1.05d0
                                if(TrialStep(j,i)>MaxTranStep(j)) then
                                    TrialStep(j,i)=MaxTranStep(j)
                                end if
                            end if
                            TrialSuccessCount(j,i)=0
                            TrialCount(j,i)=0
                        end if
                    end do
                end do
            end if
        end do
    end subroutine NVT2DA
    
    !NPT
    subroutine NPTA(Concentration,Widom)
        real*8,allocatable,dimension(:),intent(out)::Concentration,Widom
        logical::accept
        integer::i,j,k
        real*8::direction,EnthalpyOfVaporization,EnthalpyOfVaporizationTemp
        real*8,allocatable,dimension(:)::ConcentrationTemp,WidomTemp
        !Allocate memory to store your result, and set the counters to zero
            allocate(Concentration(MoleculeKinds))
            Concentration=0d0
            allocate(ConcentrationTemp(MoleculeKinds))
            ConcentrationTemp=0d0
            allocate(Widom(MoleculeKinds))
            Widom=0d0
            allocate(WidomTemp(MoleculeKinds))
            WidomTemp=0d0
        !Trial for TotalSteps times
        do k=1,TotalSteps
            !Trial move: insertion or deletion, relaxation
            call random_number(direction)
            if(direction<0.01d0) then
                call VolumeChange(accept)
            else
                call Relax3D(i,j,direction,accept)
            end if
            !Collect data
            ConcentrationTemp=ConcentrationTemp+NExist/volume
            WidomTemp=WidomTemp+Widom3D()
            !Follow evolving progress and add the analyzation
            if(mod(k,Evolve_FollowFreq)==0) then
                call showtime()
                write(*,*)'step =',k
                AverageEnergy=EnergyTotalAverageCurrentOutput3D()
                if(AverageEnergy<AverageGroundEnergy) then
                    AverageGroundEnergy=AverageEnergy
                end if
                write(*,*)'Save check point files'
                call WriteCheckPoint()
                !Add the analyzation
                Concentration=Concentration+ConcentrationTemp/TotalSteps
                Widom=Widom+WidomTemp/TotalSteps
                !Output current ensemble average
                ConcentrationTemp=Concentration*TotalSteps/k*(AInAU*AInAU*AInAU*1d27)/NAvogadro
                open(unit=99,file='Concentration.txt',status='replace')
                    write(99,*)ConcentrationTemp
                close(99)
                forall(i=1:MoleculeKinds)
                    WidomTemp(i)=-temperature*Log(Widom(i)*TotalSteps/k/ConcentrationTemp(i)&
                        *(MoleculeInput(i).mass*temperature/pim2)**1.5d0)/kJmolInAU
                end forall
                open(unit=99,file='ChemicalPotential.txt')
                    write(99,*)WidomTemp
                close(99)
                !Reset the counter
                ConcentrationTemp=0d0
                WidomTemp=0d0
            end if
            !Adjust translational trial step
            if(mod(k,Evolve_AdjustStep)==0) then
                do i=1,MoleculeKinds
                    do j=1,3
                        if(TrialCount(j,i)>Evolve_adjust) then
                            if(TrialSuccessCount(j,i)<TrialCount(j,i)/2) then
                                TrialStep(j,i)=TrialStep(j,i)*0.95d0
                            else
                                TrialStep(j,i)=TrialStep(j,i)*1.05d0
                                if(TrialStep(j,i)>MaxTranStep(j)) then
                                    TrialStep(j,i)=MaxTranStep(j)
                                end if
                            end if
                            TrialSuccessCount(j,i)=0
                            TrialCount(j,i)=0
                        end if
                    end do
                end do
            end if
        end do
    end subroutine NPTA

    !mVT
    !For 1D
    subroutine mVT1DA(MoleculeCount,Virial,MoleculeNumber)
        real*8,allocatable,dimension(:),intent(out)::MoleculeCount
        real*8,allocatable,dimension(:,:),intent(out)::Virial,MoleculeNumber
        logical::accept
        integer::i,j,k,tp,chosen,LastOne
        real*8::direction
        real*8,allocatable,dimension(:)::MoleculeCountTemp,TotalDensity
        real*8,allocatable,dimension(:,:)::VirialTemp,MoleculeNumberTemp
        !Allocate memory to store your result, and set the counters to zero
            allocate(MoleculeCount(MoleculeKinds))
            MoleculeCount=0d0
            allocate(MoleculeCountTemp(MoleculeKinds))
            MoleculeCountTemp=0d0
            allocate(Virial(2,bucketnumbers))
            Virial=0d0
            allocate(VirialTemp(2,bucketnumbers))
            VirialTemp=0d0
            allocate(MoleculeNumber(bucketnumbers,MoleculeKinds))
            MoleculeNumber=0d0
            allocate(MoleculeNumberTemp(bucketnumbers,MoleculeKinds))
            MoleculeNumberTemp=0d0
            allocate(TotalDensity(bucketnumbers))
        !Trial for TotalSteps times
        do k=1,TotalSteps
            !Trial move: insertion or deletion, relaxation
            call random_number(direction)
            if(direction<0.01d0) then
                call Insertion1D(tp,chosen,direction,accept)
                if(accept) then
                    if(direction<0.5d0) then
                        !Delete the chosen molecule by flushing it with the last exist one before the trial move
                        LastOne=NExist(tp)+1
                        !Single-body
                        mtpv(tp).mol(chosen)=mtpv(tp).mol(LastOne)
                        !Many-body
                        do j=1,tp-1
                            forall(i=1:NExist(j))
                                intertpv(tp,j).molab(chosen,i)=intertpv(tp,j).molab(LastOne,i)
                            end forall
                        end do
                        forall(i=1:chosen-1)
                            intertpv(tp,tp).molab(chosen,i)=intertpv(tp,tp).molab(LastOne,i)
                        end forall
                        !molab(LastOne,chosen) is not needed as mol(LastOne) is deleted after the trial move
                        forall(i=chosen+1:NExist(tp))
                            intertpv(tp,tp).molab(i,chosen)=intertpv(tp,tp).molab(LastOne,i)
                        end forall
                        do j=tp+1,MoleculeKinds
                            forall(i=1:NExist(j))
                                intertpv(j,tp).molab(i,chosen)=intertpv(j,tp).molab(i,LastOne)
                            end forall
                        end do
                    else
                        call UpdateChosenVirial1D(tp,chosen)
                    end if
                end if
            else
                call Relax1D(i,j,direction,accept)
                if(accept) then
                    call UpdateChosenVirial1D(i,j)
                end if
            end if
            !Collect data
            MoleculeCountTemp=MoleculeCountTemp+NExist
            VirialTemp=VirialTemp+Virial1DInEachBin()
            MoleculeNumberTemp=MoleculeNumberTemp+MoleculeNumberInEachBin()
            !Follow evolving progress and add the analyzation
            if(mod(k,Evolve_FollowFreq)==0) then
                call showtime()
                write(*,*)'step =',k
                AverageEnergy=EnergyTotalAverageCurrentOutput1D()
                if(AverageEnergy<AverageGroundEnergy) then
                    AverageGroundEnergy=AverageEnergy
                end if
                write(*,*)'Save check point files'
                call WriteCheckPoint()
                call WriteVirialCheckPoint()
                !Add the analyzation
                MoleculeCount=MoleculeCount+MoleculeCountTemp/TotalSteps
                Virial=Virial+VirialTemp/TotalSteps
                MoleculeNumber=MoleculeNumber+MoleculeNumberTemp/TotalSteps
                !Output current ensemble average
                MoleculeCountTemp=MoleculeCount*TotalSteps/k
                open(unit=99,file='MoleculeCount.txt',status='replace')
                    write(99,*)MoleculeCountTemp
                close(99)
                VirialTemp=Virial*TotalSteps/k
                MoleculeNumberTemp=MoleculeNumber*TotalSteps/k
                call MoleculeNumber2Density1D(MoleculeNumberTemp,TotalDensity)
                call Virial2Pressure1D(VirialTemp,TotalDensity)
                !Reset the counters
                MoleculeCountTemp=0d0
                VirialTemp=0d0
                MoleculeNumberTemp=0d0
            end if
            !Adjust translational trial step
            if(mod(k,Evolve_AdjustStep)==0) then
                do i=1,MoleculeKinds
                    do j=1,3
                        if(TrialCount(j,i)>Evolve_adjust) then
                            if(TrialSuccessCount(j,i)<TrialCount(j,i)/2) then
                                TrialStep(j,i)=TrialStep(j,i)*0.95d0
                            else
                                TrialStep(j,i)=TrialStep(j,i)*1.05d0
                                if(TrialStep(j,i)>MaxTranStep(j)) then
                                    TrialStep(j,i)=MaxTranStep(j)
                                end if
                            end if
                            TrialSuccessCount(j,i)=0
                            TrialCount(j,i)=0
                        end if
                    end do
                end do
            end if
        end do
    end subroutine mVT1DA
    !For 2D
    subroutine mVT2DA(MoleculeCount,Virial,MoleculeNumber)
        real*8,allocatable,dimension(:),intent(out)::MoleculeCount,Virial
        real*8,allocatable,dimension(:,:),intent(out)::MoleculeNumber
        logical::accept
        integer::i,j,k,tp,chosen,LastOne
        real*8::direction
        real*8,allocatable,dimension(:)::MoleculeCountTemp,VirialTemp,TotalDensity
        real*8,allocatable,dimension(:,:)::MoleculeNumberTemp
        !Allocate memory to store your result, and set the counters to zero
            allocate(MoleculeCount(MoleculeKinds))
            MoleculeCount=0d0
            allocate(MoleculeCountTemp(MoleculeKinds))
            MoleculeCountTemp=0d0
            allocate(Virial(bucketnumbers))
            Virial=0d0
            allocate(VirialTemp(bucketnumbers))
            VirialTemp=0d0
            allocate(MoleculeNumber(bucketnumbers,MoleculeKinds))
            MoleculeNumber=0d0
            allocate(MoleculeNumberTemp(bucketnumbers,MoleculeKinds))
            MoleculeNumberTemp=0d0
            allocate(TotalDensity(bucketnumbers))
        !Trial for TotalSteps times
        do k=1,TotalSteps
            !Trial move: insertion or deletion, relaxation
            call random_number(direction)
            if(direction<0.01d0) then
                call Insertion2D(tp,chosen,direction,accept)
                if(accept) then
                    if(direction<0.5d0) then
                        !Delete the chosen molecule by flushing it with the last exist one before the trial move
                        LastOne=NExist(tp)+1
                        !Single-body
                        mtpv(tp).mol(chosen)=mtpv(tp).mol(LastOne)
                        !Many-body
                        do j=1,tp-1
                            forall(i=1:NExist(j))
                                intertpv(tp,j).molab(chosen,i)=intertpv(tp,j).molab(LastOne,i)
                            end forall
                        end do
                        forall(i=1:chosen-1)
                            intertpv(tp,tp).molab(chosen,i)=intertpv(tp,tp).molab(LastOne,i)
                        end forall
                        !molab(LastOne,chosen) is not needed as mol(LastOne) is deleted after the trial move
                        forall(i=chosen+1:NExist(tp))
                            intertpv(tp,tp).molab(i,chosen)=intertpv(tp,tp).molab(LastOne,i)
                        end forall
                        do j=tp+1,MoleculeKinds
                            forall(i=1:NExist(j))
                                intertpv(j,tp).molab(i,chosen)=intertpv(j,tp).molab(i,LastOne)
                            end forall
                        end do
                    else
                        call UpdateChosenVirial2D(tp,chosen)
                    end if
                end if
            else
                call Relax2D(i,j,direction,accept)
                if(accept) then
                    call UpdateChosenVirial2D(i,j)
                end if
            end if
            !Collect data
            MoleculeCountTemp=MoleculeCountTemp+NExist
            VirialTemp=VirialTemp+Virial2DInEachBin()
            MoleculeNumberTemp=MoleculeNumberTemp+MoleculeNumberInEachBin()
            !Follow evolving progress and add the analyzation
            if(mod(k,Evolve_FollowFreq)==0) then
                call showtime()
                write(*,*)'step =',k
                AverageEnergy=EnergyTotalAverageCurrentOutput2D()
                if(AverageEnergy<AverageGroundEnergy) then
                    AverageGroundEnergy=AverageEnergy
                end if
                write(*,*)'Save check point files'
                call WriteCheckPoint()
                call WriteVirialCheckPoint()
                !Add the analyzation
                MoleculeCount=MoleculeCount+MoleculeCountTemp/TotalSteps
                Virial=Virial+VirialTemp/TotalSteps
                MoleculeNumber=MoleculeNumber+MoleculeNumberTemp/TotalSteps
                !Output current ensemble average
                MoleculeCountTemp=MoleculeCount*TotalSteps/k
                open(unit=99,file='MoleculeCount.txt',status='replace')
                    write(99,*)MoleculeCountTemp
                close(99)
                VirialTemp=Virial*TotalSteps/k
                MoleculeNumberTemp=MoleculeNumber*TotalSteps/k
                call MoleculeNumber2Density2D(MoleculeNumberTemp,TotalDensity)
                call Virial2Pressure2D(VirialTemp,TotalDensity)
                !Reset the counters
                MoleculeCountTemp=0d0
                VirialTemp=0d0
                MoleculeNumberTemp=0d0
            end if
            !Adjust translational trial step
            if(mod(k,Evolve_AdjustStep)==0) then
                do i=1,MoleculeKinds
                    do j=1,3
                        if(TrialCount(j,i)>Evolve_adjust) then
                            if(TrialSuccessCount(j,i)<TrialCount(j,i)/2) then
                                TrialStep(j,i)=TrialStep(j,i)*0.95d0
                            else
                                TrialStep(j,i)=TrialStep(j,i)*1.05d0
                                if(TrialStep(j,i)>MaxTranStep(j)) then
                                    TrialStep(j,i)=MaxTranStep(j)
                                end if
                            end if
                            TrialSuccessCount(j,i)=0
                            TrialCount(j,i)=0
                        end if
                    end do
                end do
            end if
        end do
    end subroutine mVT2DA
    !For 3D
    subroutine mVT3DA(MoleculeCount)
        real*8,allocatable,dimension(:),intent(out)::MoleculeCount
        logical::accept
        integer::i,j,k
        real*8::direction
        real*8,allocatable,dimension(:)::MoleculeCountTemp
        !Allocate memory to store your result, and set the counters to zero
            allocate(MoleculeCount(MoleculeKinds))
            MoleculeCount=0d0
            allocate(MoleculeCountTemp(MoleculeKinds))
            MoleculeCountTemp=0d0
        !Trial for TotalSteps times
        do k=1,TotalSteps
            !Trial move: insertion or deletion, relaxation
            call random_number(direction)
            if(direction<0.01d0) then
                call Insertion3D(i,j,direction,accept)
            else
                call Relax3D(i,j,direction,accept)
            end if
            !Collect data
            MoleculeCountTemp=MoleculeCountTemp+NExist
            !Follow evolving progress and add the analyzation
            if(mod(k,Evolve_FollowFreq)==0) then
                call showtime()
                write(*,*)'step =',k
                AverageEnergy=EnergyTotalAverageCurrentOutput3D()
                if(AverageEnergy<AverageGroundEnergy) then
                    AverageGroundEnergy=AverageEnergy
                end if
                write(*,*)'Save check point files'
                call WriteCheckPoint()
                !Add the analyzation
                MoleculeCount=MoleculeCount+MoleculeCountTemp/TotalSteps
                !Output current ensemble average
                MoleculeCountTemp=MoleculeCount*TotalSteps/k
                open(unit=99,file='MoleculeCount.txt',status='replace')
                    write(99,*)MoleculeCountTemp
                close(99)
                !Reset the counter
                MoleculeCountTemp=0d0
            end if
            !Adjust translational trial step
            if(mod(k,Evolve_AdjustStep)==0) then
                do i=1,MoleculeKinds
                    do j=1,3
                        if(TrialCount(j,i)>Evolve_adjust) then
                            if(TrialSuccessCount(j,i)<TrialCount(j,i)/2) then
                                TrialStep(j,i)=TrialStep(j,i)*0.95d0
                            else
                                TrialStep(j,i)=TrialStep(j,i)*1.05d0
                                if(TrialStep(j,i)>MaxTranStep(j)) then
                                    TrialStep(j,i)=MaxTranStep(j)
                                end if
                            end if
                            TrialSuccessCount(j,i)=0
                            TrialCount(j,i)=0
                        end if
                    end do
                end do
            end if
        end do
    end subroutine mVT3DA
    
    !RVT
    !For 1D
    subroutine RVT1DA(MoleculeCount,Virial,MoleculeNumber)
        real*8,allocatable,dimension(:),intent(out)::MoleculeCount
        real*8,allocatable,dimension(:,:),intent(out)::Virial,MoleculeNumber
        logical::accept
        integer::i,j,k,rtp,r,p,chosenr,chosenp,LastOne
        integer,dimension(2)::r2,chosenr2
        real*8::direction
        real*8,allocatable,dimension(:)::MoleculeCountTemp,TotalDensity
        real*8,allocatable,dimension(:,:)::VirialTemp,MoleculeNumberTemp
        !Allocate memory to store your result, and set the counters to zero
            allocate(MoleculeCount(MoleculeKinds))
            MoleculeCount=0d0
            allocate(MoleculeCountTemp(MoleculeKinds))
            MoleculeCountTemp=0d0
            allocate(Virial(2,bucketnumbers))
            Virial=0d0
            allocate(VirialTemp(2,bucketnumbers))
            VirialTemp=0d0
            allocate(MoleculeNumber(bucketnumbers,MoleculeKinds))
            MoleculeNumber=0d0
            allocate(MoleculeNumberTemp(bucketnumbers,MoleculeKinds))
            MoleculeNumberTemp=0d0
            allocate(TotalDensity(bucketnumbers))
        !Trial for TotalSteps times
        do k=1,TotalSteps
            !Trial move: reaction, relaxation
            call random_number(direction)
            if(direction<0.01d0) then
                !Select a reaction
                call random_number(direction)
                rtp=ceiling(direction*ReactionKinds)
                if(ReactionInput(rtp).NReactant==1) then
                    call ReactAB1D(rtp,r,p,chosenr,chosenp,direction,accept)
                    if(accept) then
                        !Delete the reactant by flushing it with the last exist one before the trial move
                        LastOne=NExist(r)+1
                        !Single-body
                        mtpv(r).mol(chosenr)=mtpv(r).mol(LastOne)
                        !Many-body
                        do j=1,r-1
                            forall(i=1:NExist(j))
                                intertpv(r,j).molab(chosenr,i)=intertpv(r,j).molab(LastOne,i)
                            end forall
                        end do
                        forall(i=1:chosenr-1)
                            intertpv(r,r).molab(chosenr,i)=intertpv(r,r).molab(LastOne,i)
                        end forall
                        !molab(LastOne,chosenr) is not needed as mol(LastOne) is deleted after the trial move
                        forall(i=chosenr+1:NExist(r))
                            intertpv(r,r).molab(i,chosenr)=intertpv(r,r).molab(LastOne,i)
                        end forall
                        do j=r+1,MoleculeKinds
                            forall(i=1:NExist(j))
                                intertpv(j,r).molab(i,chosenr)=intertpv(j,r).molab(i,LastOne)
                            end forall
                        end do
                        !Compute the virial of the product
                        call UpdateChosenVirial1D(p,chosenp)
                    end if
                else
                    call ReactABC1D(rtp,r2,p,chosenr2,chosenp,direction,accept)
                    if(accept) then
                        if(direction<0.5d0) then
                            !Delete the reactants by flushing them with the last exist ones before the trial move
                            !Swaping more than 1 molecule needs to use the real NExistOld before the trial move
                            forall(r=1:2)
                                NExistOld(r2(r))=NExist(r2(r))+1
                            end forall
                            do r=1,2
                                !Single-body
                                mtpv(r2(r)).mol(chosenr2(r))=mtpv(r2(r)).mol(NExistOld(r2(r)))
                                !Many-body
                                do j=1,r2(r)-1
                                    forall(i=1:NExistOld(j))
                                        intertpv(r2(r),j).molab(chosenr2(r),i)=intertpv(r2(r),j).molab(NExistOld(r2(r)),i)
                                    end forall
                                end do
                                forall(i=1:chosenr2(r)-1)
                                    intertpv(r2(r),r2(r)).molab(chosenr2(r),i)=intertpv(r2(r),r2(r)).molab(NExistOld(r2(r)),i)
                                end forall
                                !Same type molab(NExistOld(r),chosenr) is not needed as mol(NExistOld(r)) is deleted
                                forall(i=chosenr2(r)+1:NExist(r2(r)))
                                    intertpv(r2(r),r2(r)).molab(i,chosenr2(r))=intertpv(r2(r),r2(r)).molab(NExistOld(r2(r)),i)
                                end forall
                                do j=r2(r)+1,MoleculeKinds
                                    forall(i=1:NExistOld(j))
                                        intertpv(j,r2(r)).molab(i,chosenr2(r))=intertpv(j,r2(r)).molab(i,NExistOld(r2(r)))
                                    end forall
                                end do
                            end do
                            !Get back the NExistOld after the trial move
                            forall(r=1:2)
                                NExistOld(r2(r))=NExist(r2(r))
                            end forall
                            !Compute the virial of the product
                            call UpdateChosenVirial1D(p,chosenp)
                        else
                            !Delete the chosen molecule by flushing it with the last exist one before the trial move
                            LastOne=NExist(p)+1
                            !Single-body
                            mtpv(p).mol(chosenp)=mtpv(p).mol(LastOne)
                            !Many-body
                            do j=1,p-1
                                forall(i=1:NExist(j))
                                    intertpv(p,j).molab(chosenp,i)=intertpv(p,j).molab(LastOne,i)
                                end forall
                            end do
                            forall(i=1:chosenp-1)
                                intertpv(p,p).molab(chosenp,i)=intertpv(p,p).molab(LastOne,i)
                            end forall
                            !molab(LastOne,chosen) is not needed as mol(LastOne) is deleted after the trial move
                            forall(i=chosenp+1:NExist(p))
                                intertpv(p,p).molab(i,chosenp)=intertpv(p,p).molab(LastOne,i)
                            end forall
                            do j=p+1,MoleculeKinds
                                forall(i=1:NExist(j))
                                    intertpv(j,p).molab(i,chosenp)=intertpv(j,p).molab(i,LastOne)
                                end forall
                            end do
                            !Compute the virial of the reactants
                            call UpdateChosenVirial1D(r2(1),chosenr2(1))
                            call UpdateChosenVirial1D(r2(2),chosenr2(2))
                        end if
                    end if
                end if
            else
                call Relax1D(i,j,direction,accept)
                if(accept) then
                    call UpdateChosenVirial1D(i,j)
                end if
            end if
            !Collect data
            MoleculeCountTemp=MoleculeCountTemp+NExist
            VirialTemp=VirialTemp+Virial1DInEachBin()
            MoleculeNumberTemp=MoleculeNumberTemp+MoleculeNumberInEachBin()
            !Follow evolving progress and add the analyzation
            if(mod(k,Evolve_FollowFreq)==0) then
                call showtime()
                write(*,*)'step =',k
                AverageEnergy=EnergyTotalAverageCurrentOutput1D()
                if(AverageEnergy<AverageGroundEnergy) then
                    AverageGroundEnergy=AverageEnergy
                end if
                write(*,*)'Save check point files'
                call WriteCheckPoint()
                call WriteVirialCheckPoint()
                !Add the analyzation
                MoleculeCount=MoleculeCount+MoleculeCountTemp/TotalSteps
                Virial=Virial+VirialTemp/TotalSteps
                MoleculeNumber=MoleculeNumber+MoleculeNumberTemp/TotalSteps
                !Output current ensemble average
                MoleculeCountTemp=MoleculeCount*TotalSteps/k
                open(unit=99,file='MoleculeCount.txt',status='replace')
                    write(99,*)MoleculeCountTemp
                close(99)
                VirialTemp=Virial*TotalSteps/k
                MoleculeNumberTemp=MoleculeNumber*TotalSteps/k
                call MoleculeNumber2Density1D(MoleculeNumberTemp,TotalDensity)
                call Virial2Pressure1D(VirialTemp,TotalDensity)
                !Reset the counters
                MoleculeCountTemp=0d0
                VirialTemp=0d0
                MoleculeNumberTemp=0d0
            end if
            !Adjust translational trial step
            if(mod(k,Evolve_AdjustStep)==0) then
                do i=1,MoleculeKinds
                    do j=1,3
                        if(TrialCount(j,i)>Evolve_adjust) then
                            if(TrialSuccessCount(j,i)<TrialCount(j,i)/2) then
                                TrialStep(j,i)=TrialStep(j,i)*0.95d0
                            else
                                TrialStep(j,i)=TrialStep(j,i)*1.05d0
                                if(TrialStep(j,i)>MaxTranStep(j)) then
                                    TrialStep(j,i)=MaxTranStep(j)
                                end if
                            end if
                            TrialSuccessCount(j,i)=0
                            TrialCount(j,i)=0
                        end if
                    end do
                end do
            end if
        end do
    end subroutine RVT1DA
    !For 2D
    subroutine RVT2DA(MoleculeCount,Virial,MoleculeNumber)
        real*8,allocatable,dimension(:),intent(out)::Virial,MoleculeCount
        real*8,allocatable,dimension(:,:),intent(out)::MoleculeNumber
        logical::accept
        integer::i,j,k,rtp,r,p,chosenr,chosenp,LastOne
        integer,dimension(2)::r2,chosenr2
        real*8::direction
        real*8,allocatable,dimension(:)::VirialTemp,MoleculeCountTemp,TotalDensity
        real*8,allocatable,dimension(:,:)::MoleculeNumberTemp
        !Allocate memory to store your result, and set the counters to zero
            allocate(MoleculeCount(MoleculeKinds))
            MoleculeCount=0d0
            allocate(MoleculeCountTemp(MoleculeKinds))
            MoleculeCountTemp=0d0
            allocate(Virial(bucketnumbers))
            Virial=0d0
            allocate(VirialTemp(bucketnumbers))
            VirialTemp=0d0
            allocate(MoleculeNumber(bucketnumbers,MoleculeKinds))
            MoleculeNumber=0d0
            allocate(MoleculeNumberTemp(bucketnumbers,MoleculeKinds))
            MoleculeNumberTemp=0d0
            allocate(TotalDensity(bucketnumbers))
        !Trial for TotalSteps times
        do k=1,TotalSteps
            !Trial move: reaction, relaxation
            call random_number(direction)
            if(direction<0.01d0) then
                !Select a reaction
                call random_number(direction)
                rtp=ceiling(direction*ReactionKinds)
                if(ReactionInput(rtp).NReactant==1) then
                    call ReactAB2D(rtp,r,p,chosenr,chosenp,direction,accept)
                    if(accept) then
                        !Delete the reactant by flushing it with the last exist one before the trial move
                        LastOne=NExist(r)+1
                        !Single-body
                        mtpv(r).mol(chosenr)=mtpv(r).mol(LastOne)
                        !Many-body
                        do j=1,r-1
                            forall(i=1:NExist(j))
                                intertpv(r,j).molab(chosenr,i)=intertpv(r,j).molab(LastOne,i)
                            end forall
                        end do
                        forall(i=1:chosenr-1)
                            intertpv(r,r).molab(chosenr,i)=intertpv(r,r).molab(LastOne,i)
                        end forall
                        !molab(LastOne,chosenr) is not needed as mol(LastOne) is deleted after the trial move
                        forall(i=chosenr+1:NExist(r))
                            intertpv(r,r).molab(i,chosenr)=intertpv(r,r).molab(LastOne,i)
                        end forall
                        do j=r+1,MoleculeKinds
                            forall(i=1:NExist(j))
                                intertpv(j,r).molab(i,chosenr)=intertpv(j,r).molab(i,LastOne)
                            end forall
                        end do
                        !Compute the virial of the product
                        call UpdateChosenVirial1D(p,chosenp)
                    end if
                else
                    call ReactABC2D(rtp,r2,p,chosenr2,chosenp,direction,accept)
                    if(accept) then
                        if(direction<0.5d0) then
                            !Delete the reactants by flushing them with the last exist ones before the trial move
                            !Swaping more than 1 molecule needs to use the real NExistOld before the trial move
                            forall(r=1:2)
                                NExistOld(r2(r))=NExist(r2(r))+1
                            end forall
                            do r=1,2
                                !Single-body
                                mtpv(r2(r)).mol(chosenr2(r))=mtpv(r2(r)).mol(NExistOld(r2(r)))
                                !Many-body
                                do j=1,r2(r)-1
                                    forall(i=1:NExistOld(j))
                                        intertpv(r2(r),j).molab(chosenr2(r),i)=intertpv(r2(r),j).molab(NExistOld(r2(r)),i)
                                    end forall
                                end do
                                forall(i=1:chosenr2(r)-1)
                                    intertpv(r2(r),r2(r)).molab(chosenr2(r),i)=intertpv(r2(r),r2(r)).molab(NExistOld(r2(r)),i)
                                end forall
                                !Same type molab(NExistOld(r),chosenr) is not needed as mol(NExistOld(r)) is deleted
                                forall(i=chosenr2(r)+1:NExist(r2(r)))
                                    intertpv(r2(r),r2(r)).molab(i,chosenr2(r))=intertpv(r2(r),r2(r)).molab(NExistOld(r2(r)),i)
                                end forall
                                do j=r2(r)+1,MoleculeKinds
                                    forall(i=1:NExistOld(j))
                                        intertpv(j,r2(r)).molab(i,chosenr2(r))=intertpv(j,r2(r)).molab(i,NExistOld(r2(r)))
                                    end forall
                                end do
                            end do
                            !Get back the NExistOld after the trial move
                            forall(r=1:2)
                                NExistOld(r2(r))=NExist(r2(r))
                            end forall
                            !Compute the virial of the product
                            call UpdateChosenVirial1D(p,chosenp)
                        else
                            !Delete the chosen molecule by flushing it with the last exist one before the trial move
                            LastOne=NExist(p)+1
                            !Single-body
                            mtpv(p).mol(chosenp)=mtpv(p).mol(LastOne)
                            !Many-body
                            do j=1,p-1
                                forall(i=1:NExist(j))
                                    intertpv(p,j).molab(chosenp,i)=intertpv(p,j).molab(LastOne,i)
                                end forall
                            end do
                            forall(i=1:chosenp-1)
                                intertpv(p,p).molab(chosenp,i)=intertpv(p,p).molab(LastOne,i)
                            end forall
                            !molab(LastOne,chosen) is not needed as mol(LastOne) is deleted after the trial move
                            forall(i=chosenp+1:NExist(p))
                                intertpv(p,p).molab(i,chosenp)=intertpv(p,p).molab(LastOne,i)
                            end forall
                            do j=p+1,MoleculeKinds
                                forall(i=1:NExist(j))
                                    intertpv(j,p).molab(i,chosenp)=intertpv(j,p).molab(i,LastOne)
                                end forall
                            end do
                            !Compute the virial of the reactants
                            call UpdateChosenVirial1D(r2(1),chosenr2(1))
                            call UpdateChosenVirial1D(r2(2),chosenr2(2))
                        end if
                    end if
                end if
            else
                call Relax2D(i,j,direction,accept)
                if(accept) then
                    call UpdateChosenVirial2D(i,j)
                end if
            end if
            !Collect data
            MoleculeCountTemp=MoleculeCountTemp+NExist
            VirialTemp=VirialTemp+Virial2DInEachBin()
            MoleculeNumberTemp=MoleculeNumberTemp+MoleculeNumberInEachBin()
            !Follow evolving progress and add the analyzation
            if(mod(k,Evolve_FollowFreq)==0) then
                call showtime()
                write(*,*)'step =',k
                AverageEnergy=EnergyTotalAverageCurrentOutput2D()
                if(AverageEnergy<AverageGroundEnergy) then
                    AverageGroundEnergy=AverageEnergy
                end if
                write(*,*)'Save check point files'
                call WriteCheckPoint()
                call WriteVirialCheckPoint()
                !Add the analyzation
                MoleculeCount=MoleculeCount+MoleculeCountTemp/TotalSteps
                Virial=Virial+VirialTemp/TotalSteps
                MoleculeNumber=MoleculeNumber+MoleculeNumberTemp/TotalSteps
                !Output current ensemble average
                MoleculeCountTemp=MoleculeCount*TotalSteps/k
                open(unit=99,file='MoleculeCount.txt',status='replace')
                    write(99,*)MoleculeCountTemp
                close(99)
                VirialTemp=Virial*TotalSteps/k
                MoleculeNumberTemp=MoleculeNumber*TotalSteps/k
                call MoleculeNumber2Density2D(MoleculeNumberTemp,TotalDensity)
                call Virial2Pressure2D(VirialTemp,TotalDensity)
                !Reset the counters
                MoleculeCountTemp=0d0
                VirialTemp=0d0
                MoleculeNumberTemp=0d0
            end if
            !Adjust translational trial step
            if(mod(k,Evolve_AdjustStep)==0) then
                do i=1,MoleculeKinds
                    do j=1,3
                        if(TrialCount(j,i)>Evolve_adjust) then
                            if(TrialSuccessCount(j,i)<TrialCount(j,i)/2) then
                                TrialStep(j,i)=TrialStep(j,i)*0.95d0
                            else
                                TrialStep(j,i)=TrialStep(j,i)*1.05d0
                                if(TrialStep(j,i)>MaxTranStep(j)) then
                                    TrialStep(j,i)=MaxTranStep(j)
                                end if
                            end if
                            TrialSuccessCount(j,i)=0
                            TrialCount(j,i)=0
                        end if
                    end do
                end do
            end if
        end do
    end subroutine RVT2DA
    !For 3D
    subroutine RVT3DA(MoleculeCount)
        real*8,allocatable,dimension(:),intent(out)::MoleculeCount
        logical::accept
        integer::i,j,k,chosenp
        integer,dimension(2)::r,chosenr
        real*8::direction
        real*8,allocatable,dimension(:)::MoleculeCountTemp
        !Allocate memory to store your result, and set the counters to zero
            allocate(MoleculeCount(MoleculeKinds))
            MoleculeCount=0d0
            allocate(MoleculeCountTemp(MoleculeKinds))
            MoleculeCountTemp=0d0
        !Trial for TotalSteps times
        do k=1,TotalSteps
            !Trial move: reaction, relaxation
            call random_number(direction)
            if(direction<0.01d0) then
                !Select a reaction
                call random_number(direction)
                i=ceiling(direction*ReactionKinds)
                if(ReactionInput(i).NReactant==1) then
                    call ReactAB3D(i,r(1),r(2),chosenr(1),chosenr(2),direction,accept)
                else
                    call ReactABC3D(i,r,j,chosenr,chosenp,direction,accept)
                end if
            else
                call Relax3D(i,j,direction,accept)
            end if
            !Collect data
            MoleculeCountTemp=MoleculeCountTemp+NExist
            !Follow evolving progress and add the analyzation
            if(mod(k,Evolve_FollowFreq)==0) then
                call showtime()
                write(*,*)'step =',k
                AverageEnergy=EnergyTotalAverageCurrentOutput3D()
                if(AverageEnergy<AverageGroundEnergy) then
                    AverageGroundEnergy=AverageEnergy
                end if
                write(*,*)'Save check point files'
                call WriteCheckPoint()
                !Add the analyzation
                MoleculeCount=MoleculeCount+MoleculeCountTemp/TotalSteps
                !Output current ensemble average
                MoleculeCountTemp=MoleculeCount*TotalSteps/k
                open(unit=99,file='MoleculeCount.txt',status='replace')
                    write(99,*)MoleculeCountTemp
                close(99)
                !Reset the counter
                MoleculeCountTemp=0d0
            end if
            !Adjust translational trial step
            if(mod(k,Evolve_AdjustStep)==0) then
                do i=1,MoleculeKinds
                    do j=1,3
                        if(TrialCount(j,i)>Evolve_adjust) then
                            if(TrialSuccessCount(j,i)<TrialCount(j,i)/2) then
                                TrialStep(j,i)=TrialStep(j,i)*0.95d0
                            else
                                TrialStep(j,i)=TrialStep(j,i)*1.05d0
                                if(TrialStep(j,i)>MaxTranStep(j)) then
                                    TrialStep(j,i)=MaxTranStep(j)
                                end if
                            end if
                            TrialSuccessCount(j,i)=0
                            TrialCount(j,i)=0
                        end if
                    end do
                end do
            end if
        end do
    end subroutine RVT3DA

    !RPT
    subroutine RPTA(Concentration,Widom)
        real*8,allocatable,dimension(:),intent(out)::Concentration,Widom
        logical::accept
        integer::i,j,k,chosenp
        integer,dimension(2)::r,chosenr
        real*8::direction
        real*8,allocatable,dimension(:)::ConcentrationTemp,WidomTemp
        !Allocate memory to store your result, and set the counters to zero
            allocate(Concentration(MoleculeKinds))
            Concentration=0d0
            allocate(ConcentrationTemp(MoleculeKinds))
            ConcentrationTemp=0d0
            allocate(Widom(MoleculeKinds))
            Widom=0d0
            allocate(WidomTemp(MoleculeKinds))
            WidomTemp=0d0
        !Trial for TotalSteps times
        do k=1,TotalSteps
            !Trial move: volume change, reaction, relaxation
            call random_number(direction)
            if(direction<0.01d0) then
                call VolumeChange(accept)
            else if(direction<0.02d0) then
                !Select a reaction
                call random_number(direction)
                i=ceiling(direction*ReactionKinds)
                if(ReactionInput(i).NReactant==1) then
                    call ReactAB3D(i,r(1),r(2),chosenr(1),chosenr(2),direction,accept)
                else
                    call ReactABC3D(i,r,j,chosenr,chosenp,direction,accept)
                end if
            else
                call Relax3D(i,j,direction,accept)
            end if
            !Collect data
            ConcentrationTemp=ConcentrationTemp/volume+NExist
            WidomTemp=WidomTemp+Widom3D()
            !Follow evolving progress and add the analyzation
            if(mod(k,Evolve_FollowFreq)==0) then
                call showtime()
                write(*,*)'step =',k
                AverageEnergy=EnergyTotalAverageCurrentOutput3D()
                if(AverageEnergy<AverageGroundEnergy) then
                    AverageGroundEnergy=AverageEnergy
                end if
                write(*,*)'Save check point files'
                call WriteCheckPoint()
                !Add the analyzation
                Concentration=Concentration+ConcentrationTemp/TotalSteps
                Widom=Widom+WidomTemp/TotalSteps
                !Output current ensemble average
                ConcentrationTemp=Concentration*TotalSteps/k*(AInAU*AInAU*AInAU*1d27)/NAvogadro
                open(unit=99,file='Concentration.txt',status='replace')
                    write(99,*)ConcentrationTemp
                close(99)
                forall(i=1:MoleculeKinds)
                    WidomTemp(i)=-temperature*Log(Widom(i)*TotalSteps/k/ConcentrationTemp(i)&
                        *(MoleculeInput(i).mass*temperature/pim2)**1.5d0)/kJmolInAU
                end forall
                open(unit=99,file='ChemicalPotential.txt')
                    write(99,*)WidomTemp
                close(99)
                !Reset the counter
                ConcentrationTemp=0d0
            end if
            !Adjust translational trial step
            if(mod(k,Evolve_AdjustStep)==0) then
                do i=1,MoleculeKinds
                    do j=1,3
                        if(TrialCount(j,i)>Evolve_adjust) then
                            if(TrialSuccessCount(j,i)<TrialCount(j,i)/2) then
                                TrialStep(j,i)=TrialStep(j,i)*0.95d0
                            else
                                TrialStep(j,i)=TrialStep(j,i)*1.05d0
                                if(TrialStep(j,i)>MaxTranStep(j)) then
                                    TrialStep(j,i)=MaxTranStep(j)
                                end if
                            end if
                            TrialSuccessCount(j,i)=0
                            TrialCount(j,i)=0
                        end if
                    end do
                end do
            end if
        end do
    end subroutine RPTA    
!-------------------- End ---------------------

end module Evolve