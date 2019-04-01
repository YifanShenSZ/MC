!This is the module to account for configurational energy
!Currently it only has the nonbonded part (intermolecular interaction and outer field)
!Flag: probably will add vibration some day
module ConfigurationalEnergy
    use General
    use MCBasis
    use Nonbonded
    implicit none

!Configurational energy module only variable
    real*8,allocatable,dimension(:)::ConfigurationalEnergy_ESSBDependent2DComponent,ConfigurationalEnergy_WallLJComponent
    real*8,allocatable,dimension(:,:)::ConfigurationalEnergy_ESComponent,ConfigurationalEnergy_LJComponent

contains
!The initializer for ConfigurationalEnergy module
subroutine InitializeConfigurationalEnergy()
    if(PeriodicDims==2) then
        allocate(ConfigurationalEnergy_ESSBDependent2DComponent(MoleculeKinds))
    end if
    allocate(ConfigurationalEnergy_WallLJComponent(MoleculeKinds))
    allocate(ConfigurationalEnergy_ESComponent(MoleculeKinds,MoleculeKinds))
    allocate(ConfigurationalEnergy_LJComponent(MoleculeKinds,MoleculeKinds))
    call InitializeNonbonded()
end subroutine InitializeConfigurationalEnergy

!Update the chosen molecule's interaction and return its energy
!For 1D
function EnergyChosenUpdated1D(tp,chosen)
    real*8::EnergyChosenUpdated1D
    integer,intent(in)::tp,chosen
    integer::i,j
    !Single-body
    !Currently WallLJ is the only 1D position dependent single-body term, so we do not refer to Total
    mtp(tp).mol(chosen).WallLJ=WallLJPotential1D(tp,mtp(tp).mol(chosen).Coordinate)
    EnergyChosenUpdated1D=mtp(tp).mol(chosen).WallLJ+ESSBIndependent(tp)
    !Many-body
    !Interaction between smaller type molecules and chosen molecule
    do j=1,tp-1
        do i=1,NExist(j)
            intertp(tp,j).molab(chosen,i).Total=NonbondedManyBody1D(tp,j,mtp(tp).mol(chosen).Coordinate,mtp(j).mol(i).Coordinate,&
                intertp(tp,j).molab(chosen,i).ESMB,intertp(tp,j).molab(chosen,i).LJ)
            EnergyChosenUpdated1D=EnergyChosenUpdated1D+intertp(tp,j).molab(chosen,i).Total
        end do
    end do
    !Interaction between same type smaller order molecules and chosen molecule
    do i=1,chosen-1
        intertp(tp,tp).molab(chosen,i).Total=NonbondedManyBody1D(tp,tp,mtp(tp).mol(chosen).Coordinate,mtp(tp).mol(i).Coordinate,&
            intertp(tp,tp).molab(chosen,i).ESMB,intertp(tp,tp).molab(chosen,i).LJ)
        EnergyChosenUpdated1D=EnergyChosenUpdated1D+intertp(tp,tp).molab(chosen,i).Total
    end do
    !Interaction between same type larger order molecules and chosen molecule
    do i=chosen+1,NExist(tp)
        intertp(tp,tp).molab(i,chosen).Total=NonbondedManyBody1D(tp,tp,mtp(tp).mol(i).Coordinate,mtp(tp).mol(chosen).Coordinate,&
            intertp(tp,tp).molab(i,chosen).ESMB,intertp(tp,tp).molab(i,chosen).LJ)
        EnergyChosenUpdated1D=EnergyChosenUpdated1D+intertp(tp,tp).molab(i,chosen).Total
    end do
    !Interaction between larger type molecules and chosen molecule
    do j=tp+1,MoleculeKinds
        do i=1,NExist(j)
            intertp(j,tp).molab(i,chosen).Total=NonbondedManyBody1D(j,tp,mtp(j).mol(i).Coordinate,mtp(tp).mol(chosen).Coordinate,&
                intertp(j,tp).molab(i,chosen).ESMB,intertp(j,tp).molab(i,chosen).LJ)
            EnergyChosenUpdated1D=EnergyChosenUpdated1D+intertp(j,tp).molab(i,chosen).Total
        end do
    end do
end function EnergyChosenUpdated1D
!For 2D
function EnergyChosenUpdated2D(tp,chosen)
    real*8::EnergyChosenUpdated2D
    integer,intent(in)::tp,chosen
    integer::i,j
    !Single-body
    mtp(tp).mol(chosen).ESSBDependent2D=EwaldSummationSingleBodyPositionDependent(tp,mtp(tp).mol(chosen).Coordinate(3,:))
    mtp(tp).mol(chosen).WallLJ=WallLJPotential2D(tp,mtp(tp).mol(chosen).Coordinate)
    mtp(tp).mol(chosen).Total=mtp(tp).mol(chosen).ESSBDependent2D+mtp(tp).mol(chosen).WallLJ
    EnergyChosenUpdated2D=mtp(tp).mol(chosen).Total+ESSBIndependent(tp)
    !Many-body
    !Interaction between smaller type molecules and chosen molecule
    do j=1,tp-1
        do i=1,NExist(j)
            intertp(tp,j).molab(chosen,i).Total=NonbondedManyBody2D(tp,j,mtp(tp).mol(chosen).Coordinate,mtp(j).mol(i).Coordinate,&
                intertp(tp,j).molab(chosen,i).ESMB,intertp(tp,j).molab(chosen,i).LJ)
            EnergyChosenUpdated2D=EnergyChosenUpdated2D+intertp(tp,j).molab(chosen,i).Total
        end do
    end do
    !Interaction between same type smaller order molecules and chosen molecule
    do i=1,chosen-1
        intertp(tp,tp).molab(chosen,i).Total=NonbondedManyBody2D(tp,tp,mtp(tp).mol(chosen).Coordinate,mtp(tp).mol(i).Coordinate,&
            intertp(tp,tp).molab(chosen,i).ESMB,intertp(tp,tp).molab(chosen,i).LJ)
        EnergyChosenUpdated2D=EnergyChosenUpdated2D+intertp(tp,tp).molab(chosen,i).Total
    end do
    !Interaction between same type larger order molecules and chosen molecule
    do i=chosen+1,NExist(tp)
        intertp(tp,tp).molab(i,chosen).Total=NonbondedManyBody2D(tp,tp,mtp(tp).mol(i).Coordinate,mtp(tp).mol(chosen).Coordinate,&
            intertp(tp,tp).molab(i,chosen).ESMB,intertp(tp,tp).molab(i,chosen).LJ)
        EnergyChosenUpdated2D=EnergyChosenUpdated2D+intertp(tp,tp).molab(i,chosen).Total
    end do
    !Interaction between larger type molecules and chosen molecule
    do j=tp+1,MoleculeKinds
        do i=1,NExist(j)
            intertp(j,tp).molab(i,chosen).Total=NonbondedManyBody2D(j,tp,mtp(j).mol(i).Coordinate,mtp(tp).mol(chosen).Coordinate,&
                intertp(j,tp).molab(i,chosen).ESMB,intertp(j,tp).molab(i,chosen).LJ)
            EnergyChosenUpdated2D=EnergyChosenUpdated2D+intertp(j,tp).molab(i,chosen).Total
        end do
    end do
end function EnergyChosenUpdated2D
!For 3D
function EnergyChosenUpdated3D(tp,chosen)
    real*8::EnergyChosenUpdated3D
    integer,intent(in)::tp,chosen
    integer::i,j
    !Single-body
    !Currently no position dependent single-body interaction for 3D
    EnergyChosenUpdated3D=ESSBIndependent(tp)
    !Many-body
    !Interaction between smaller type molecules and chosen molecule
    do j=1,tp-1
        do i=1,NExist(j)
            intertp(tp,j).molab(chosen,i).Total=NonbondedManyBody3D(tp,j,mtp(tp).mol(chosen).Coordinate,mtp(j).mol(i).Coordinate,&
                intertp(tp,j).molab(chosen,i).ESMB,intertp(tp,j).molab(chosen,i).LJ)
            EnergyChosenUpdated3D=EnergyChosenUpdated3D+intertp(tp,j).molab(chosen,i).Total
        end do
    end do
    !Interaction between same type smaller order molecules and chosen molecule
    do i=1,chosen-1
        intertp(tp,tp).molab(chosen,i).Total=NonbondedManyBody3D(tp,tp,mtp(tp).mol(chosen).Coordinate,mtp(tp).mol(i).Coordinate,&
            intertp(tp,tp).molab(chosen,i).ESMB,intertp(tp,tp).molab(chosen,i).LJ)
        EnergyChosenUpdated3D=EnergyChosenUpdated3D+intertp(tp,tp).molab(chosen,i).Total
    end do
    !Interaction between same type larger order molecules and chosen molecule
    do i=chosen+1,NExist(tp)
        intertp(tp,tp).molab(i,chosen).Total=NonbondedManyBody3D(tp,tp,mtp(tp).mol(i).Coordinate,mtp(tp).mol(chosen).Coordinate,&
            intertp(tp,tp).molab(i,chosen).ESMB,intertp(tp,tp).molab(i,chosen).LJ)
        EnergyChosenUpdated3D=EnergyChosenUpdated3D+intertp(tp,tp).molab(i,chosen).Total
    end do
    !Interaction between larger type molecules and chosen molecule
    do j=tp+1,MoleculeKinds
        do i=1,NExist(j)
            intertp(j,tp).molab(i,chosen).Total=NonbondedManyBody3D(j,tp,mtp(j).mol(i).Coordinate,mtp(tp).mol(chosen).Coordinate,&
                intertp(j,tp).molab(i,chosen).ESMB,intertp(j,tp).molab(i,chosen).LJ)
            EnergyChosenUpdated3D=EnergyChosenUpdated3D+intertp(j,tp).molab(i,chosen).Total
        end do
    end do
end function EnergyChosenUpdated3D

!Current energy of the chosen molecule
!For 1D
function EnergyChosenCurrent1D(tp,chosen)
    real*8::EnergyChosenCurrent1D
    integer,intent(in)::tp,chosen
    integer::i,j
    !Single-body
    !Currently WallLJ is the only 1D position dependent single-body term, so we do not refer to Total
    EnergyChosenCurrent1D=mtpold(tp).mol(chosen).WallLJ+ESSBIndependent(tp)
    !Many-body
    !Interaction between smaller type molecules and chosen molecule
    do j=1,tp-1
        do i=1,NExistOld(j)
            EnergyChosenCurrent1D=EnergyChosenCurrent1D+intertpold(tp,j).molab(chosen,i).Total
        end do
    end do
    !Interaction between same type smaller order molecules and chosen molecule
    do i=1,chosen-1
        EnergyChosenCurrent1D=EnergyChosenCurrent1D+intertpold(tp,tp).molab(chosen,i).Total
    end do
    !Interaction between same type larger order molecules and chosen molecule
    do i=chosen+1,NExistOld(tp)
        EnergyChosenCurrent1D=EnergyChosenCurrent1D+intertpold(tp,tp).molab(i,chosen).Total
    end do
    !Interaction between larger type molecules and chosen molecule
    do j=tp+1,MoleculeKinds
        do i=1,NExistOld(j)
            EnergyChosenCurrent1D=EnergyChosenCurrent1D+intertpold(j,tp).molab(i,chosen).Total
        end do
    end do
end function EnergyChosenCurrent1D
!For 2D
function EnergyChosenCurrent2D(tp,chosen)
    real*8::EnergyChosenCurrent2D
    integer,intent(in)::tp,chosen
    integer::i,j
    !Single-body
    EnergyChosenCurrent2D=mtpold(tp).mol(chosen).Total+ESSBIndependent(tp)
    !Many-body
    !Interaction between smaller type molecules and chosen molecule
    do j=1,tp-1
        do i=1,NExistOld(j)
            EnergyChosenCurrent2D=EnergyChosenCurrent2D+intertpold(tp,j).molab(chosen,i).Total
        end do
    end do
    !Interaction between same type smaller order molecules and chosen molecule
    do i=1,chosen-1
        EnergyChosenCurrent2D=EnergyChosenCurrent2D+intertpold(tp,tp).molab(chosen,i).Total
    end do
    !Interaction between same type larger order molecules and chosen molecule
    do i=chosen+1,NExistOld(tp)
        EnergyChosenCurrent2D=EnergyChosenCurrent2D+intertpold(tp,tp).molab(i,chosen).Total
    end do
    !Interaction between larger type molecules and chosen molecule
    do j=tp+1,MoleculeKinds
        do i=1,NExistOld(j)
            EnergyChosenCurrent2D=EnergyChosenCurrent2D+intertpold(j,tp).molab(i,chosen).Total
        end do
    end do
end function EnergyChosenCurrent2D
!For 3D
function EnergyChosenCurrent3D(tp,chosen)
    real*8::EnergyChosenCurrent3D
    integer,intent(in)::tp,chosen
    integer::i,j
    !Single-body
    !Currently no position dependent single-body interaction for 3D
    EnergyChosenCurrent3D=ESSBIndependent(tp)
    !Many-body
    !Interaction between smaller type molecules and chosen molecule
    do j=1,tp-1
        do i=1,NExistOld(j)
            EnergyChosenCurrent3D=EnergyChosenCurrent3D+intertpold(tp,j).molab(chosen,i).Total
        end do
    end do
    !Interaction between same type smaller order molecules and chosen molecule
    do i=1,chosen-1
        EnergyChosenCurrent3D=EnergyChosenCurrent3D+intertpold(tp,tp).molab(chosen,i).Total
    end do
    !Interaction between same type larger order molecules and chosen molecule
    do i=chosen+1,NExistOld(tp)
        EnergyChosenCurrent3D=EnergyChosenCurrent3D+intertpold(tp,tp).molab(i,chosen).Total
    end do
    !Interaction between larger type molecules and chosen molecule
    do j=tp+1,MoleculeKinds
        do i=1,NExistOld(j)
            EnergyChosenCurrent3D=EnergyChosenCurrent3D+intertpold(j,tp).molab(i,chosen).Total
        end do
    end do
end function EnergyChosenCurrent3D

!Update all molecules' interaction and return the average energy
!For 1D
function EnergyTotalAverageUpdated1D()
    real*8::EnergyTotalAverageUpdated1D
    integer::i,j,ii,jj
    !Initialization: set the counter with the single-body position independent term
    EnergyTotalAverageUpdated1D=dot_product(NExist,ESSBIndependent)
    do j=1,MoleculeKinds
        do i=1,NExist(j)
            !Single-body
            !Currently WallLJ is the only 1D position dependent single-body term, so we do not refer to Total
            mtp(j).mol(i).WallLJ=WallLJPotential1D(j,mtp(j).mol(i).Coordinate)
            EnergyTotalAverageUpdated1D=EnergyTotalAverageUpdated1D+mtp(j).mol(i).WallLJ
            !Many-body
            do ii=i+1,NExist(j)
                intertp(j,j).molab(ii,i).Total=NonbondedManyBody1D(j,j,mtp(j).mol(ii).Coordinate,mtp(j).mol(i).Coordinate,&
                    intertp(j,j).molab(ii,i).ESMB,intertp(j,j).molab(ii,i).LJ)
                EnergyTotalAverageUpdated1D=EnergyTotalAverageUpdated1D+intertp(j,j).molab(ii,i).Total
            end do
            do jj=j+1,MoleculeKinds
                do ii=1,NExist(jj)
                    intertp(jj,j).molab(ii,i).Total=NonbondedManyBody1D(jj,j,mtp(jj).mol(ii).Coordinate,mtp(j).mol(i).Coordinate,&
                        intertp(jj,j).molab(ii,i).ESMB,intertp(jj,j).molab(ii,i).LJ)
                    EnergyTotalAverageUpdated1D=EnergyTotalAverageUpdated1D+intertp(jj,j).molab(ii,i).Total
                end do
            end do
        end do
    end do
    !Take average
    EnergyTotalAverageUpdated1D=EnergyTotalAverageUpdated1D/NTotal
end function EnergyTotalAverageUpdated1D
!For 2D
function EnergyTotalAverageUpdated2D()
    real*8::EnergyTotalAverageUpdated2D
    integer::i,j,ii,jj
    !Initialization: set the counter with the single-body position independent term
    EnergyTotalAverageUpdated2D=dot_product(NExist,ESSBIndependent)
    do j=1,MoleculeKinds
        do i=1,NExist(j)
            !Single-body
            mtp(j).mol(i).ESSBDependent2D=EwaldSummationSingleBodyPositionDependent(j,mtp(j).mol(i).Coordinate(3,:))
            mtp(j).mol(i).WallLJ=WallLJPotential2D(j,mtp(j).mol(i).Coordinate)
            mtp(j).mol(i).Total=mtp(j).mol(i).ESSBDependent2D+mtp(j).mol(i).WallLJ
            EnergyTotalAverageUpdated2D=EnergyTotalAverageUpdated2D+mtp(j).mol(i).Total
            !Many-body
            do ii=i+1,NExist(j)
                intertp(j,j).molab(ii,i).Total=NonbondedManyBody2D(j,j,mtp(j).mol(ii).Coordinate,mtp(j).mol(i).Coordinate,&
                    intertp(j,j).molab(ii,i).ESMB,intertp(j,j).molab(ii,i).LJ)
                EnergyTotalAverageUpdated2D=EnergyTotalAverageUpdated2D+intertp(j,j).molab(ii,i).Total
            end do
            do jj=j+1,MoleculeKinds
                do ii=1,NExist(jj)
                    intertp(jj,j).molab(ii,i).Total=NonbondedManyBody2D(jj,j,mtp(jj).mol(ii).Coordinate,mtp(j).mol(i).Coordinate,&
                        intertp(jj,j).molab(ii,i).ESMB,intertp(jj,j).molab(ii,i).LJ)
                    EnergyTotalAverageUpdated2D=EnergyTotalAverageUpdated2D+intertp(jj,j).molab(ii,i).Total
                end do
            end do
        end do
    end do
    !Take average
    EnergyTotalAverageUpdated2D=EnergyTotalAverageUpdated2D/NTotal
end function EnergyTotalAverageUpdated2D
!For 3D
function EnergyTotalAverageUpdated3D()
    real*8::EnergyTotalAverageUpdated3D
    integer::i,j,ii,jj
    !Initialization: set the counter with the single-body position independent term
    EnergyTotalAverageUpdated3D=dot_product(NExist,ESSBIndependent)
    do j=1,MoleculeKinds
        do i=1,NExist(j)
            !Currently no position dependent single-body interaction for 3D
            !Many-body
            do ii=i+1,NExist(j)
                intertp(j,j).molab(ii,i).Total=NonbondedManyBody3D(j,j,mtp(j).mol(ii).Coordinate,mtp(j).mol(i).Coordinate,&
                    intertp(j,j).molab(ii,i).ESMB,intertp(j,j).molab(ii,i).LJ)
                EnergyTotalAverageUpdated3D=EnergyTotalAverageUpdated3D+intertp(j,j).molab(ii,i).Total
            end do
            do jj=j+1,MoleculeKinds
                do ii=1,NExist(jj)
                    intertp(jj,j).molab(ii,i).Total=NonbondedManyBody3D(jj,j,mtp(jj).mol(ii).Coordinate,mtp(j).mol(i).Coordinate,&
                        intertp(jj,j).molab(ii,i).ESMB,intertp(jj,j).molab(ii,i).LJ)
                    EnergyTotalAverageUpdated3D=EnergyTotalAverageUpdated3D+intertp(jj,j).molab(ii,i).Total
                end do
            end do
        end do
    end do
    !Take average
    EnergyTotalAverageUpdated3D=EnergyTotalAverageUpdated3D/NTotal
end function EnergyTotalAverageUpdated3D

!Current average energy of the whole system
!For 1D
function EnergyTotalAverageCurrent1D()
    real*8::EnergyTotalAverageCurrent1D
    integer::i,j,ii,jj
    !Initialization: set the counter with the single-body position independent term
    EnergyTotalAverageCurrent1D=dot_product(NExist,ESSBIndependent)
    do j=1,MoleculeKinds
        do i=1,NExist(j)
            !Single-body
            !Currently WallLJ is the only 1D position dependent single-body term, so we do not refer to Total
            EnergyTotalAverageCurrent1D=EnergyTotalAverageCurrent1D+mtpold(j).mol(i).WallLJ
            !Many-body
            do ii=i+1,NExist(j)
                EnergyTotalAverageCurrent1D=EnergyTotalAverageCurrent1D+intertpold(j,j).molab(ii,i).Total
            end do
            do jj=j+1,MoleculeKinds
                do ii=1,NExist(jj)
                    EnergyTotalAverageCurrent1D=EnergyTotalAverageCurrent1D+intertpold(jj,j).molab(ii,i).Total
                end do
            end do
        end do
    end do
    !Take average
    EnergyTotalAverageCurrent1D=EnergyTotalAverageCurrent1D/NTotal
end function EnergyTotalAverageCurrent1D
!For 2D
function EnergyTotalAverageCurrent2D()
    real*8::EnergyTotalAverageCurrent2D
    integer::i,j,ii,jj
    !Initialization: set the counter with the single-body position independent term
    EnergyTotalAverageCurrent2D=dot_product(NExist,ESSBIndependent)
    do j=1,MoleculeKinds
        do i=1,NExist(j)
            !Single-body
            EnergyTotalAverageCurrent2D=EnergyTotalAverageCurrent2D+mtpold(j).mol(i).Total
            !Many-body
            do ii=i+1,NExist(j)
                EnergyTotalAverageCurrent2D=EnergyTotalAverageCurrent2D+intertpold(j,j).molab(ii,i).Total
            end do
            do jj=j+1,MoleculeKinds
                do ii=1,NExist(jj)
                    EnergyTotalAverageCurrent2D=EnergyTotalAverageCurrent2D+intertpold(jj,j).molab(ii,i).Total
                end do
            end do
        end do
    end do
    !Take average
    EnergyTotalAverageCurrent2D=EnergyTotalAverageCurrent2D/NTotal
end function EnergyTotalAverageCurrent2D
!For 3D
function EnergyTotalAverageCurrent3D()
    real*8::EnergyTotalAverageCurrent3D
    integer::i,j,ii,jj
    !Initialization: set the counter with the single-body position independent term
    EnergyTotalAverageCurrent3D=dot_product(NExist,ESSBIndependent)
    do j=1,MoleculeKinds
        do i=1,NExist(j)
            !Currently no position dependent single-body interaction for 3D
            !Many-body
            do ii=i+1,NExist(j)
                EnergyTotalAverageCurrent3D=EnergyTotalAverageCurrent3D+intertpold(j,j).molab(ii,i).Total
            end do
            do jj=j+1,MoleculeKinds
                do ii=1,NExist(jj)
                    EnergyTotalAverageCurrent3D=EnergyTotalAverageCurrent3D+intertpold(jj,j).molab(ii,i).Total
                end do
            end do
        end do
    end do
    !Take average
    EnergyTotalAverageCurrent3D=EnergyTotalAverageCurrent3D/NTotal
end function EnergyTotalAverageCurrent3D

!Current average energy of the whole system, and output every component
!For 1D
function EnergyTotalAverageCurrentOutput1D()
    real*8::EnergyTotalAverageCurrentOutput1D
    integer::i,j,ii,jj
    !Initialization: set the counters to zero
    EnergyTotalAverageCurrentOutput1D=0d0
    ConfigurationalEnergy_ESComponent=0d0
    ConfigurationalEnergy_WallLJComponent=0d0
    ConfigurationalEnergy_LJComponent=0d0
    !Sum up total
    do j=1,MoleculeKinds
        do i=1,NExist(j)
            !Single-type
                !Single-body
                ConfigurationalEnergy_WallLJComponent(j)=ConfigurationalEnergy_WallLJComponent(j)+mtp(j).mol(i).WallLJ
                !Many-body
                do ii=i+1,NExist(j)
                    ConfigurationalEnergy_ESComponent(j,j)=ConfigurationalEnergy_ESComponent(j,j)+intertp(j,j).molab(ii,i).ESMB
                    ConfigurationalEnergy_LJComponent(j,j)=ConfigurationalEnergy_LJComponent(j,j)+intertp(j,j).molab(ii,i).LJ
                end do
            !Many-type
            !Many-body
            do jj=j+1,MoleculeKinds
                do ii=1,NExist(jj)
                    ConfigurationalEnergy_ESComponent(jj,j)=ConfigurationalEnergy_ESComponent(jj,j)+intertp(jj,j).molab(ii,i).ESMB
                    ConfigurationalEnergy_LJComponent(jj,j)=ConfigurationalEnergy_LJComponent(jj,j)+intertp(jj,j).molab(ii,i).LJ
                end do
            end do
        end do
        !Single-type
        !Add single-body term to electro static interaction
        ConfigurationalEnergy_ESComponent(j,j)=ConfigurationalEnergy_ESComponent(j,j)+NExist(j)*ESSBIndependent(j)
        EnergyTotalAverageCurrentOutput1D=EnergyTotalAverageCurrentOutput1D+ConfigurationalEnergy_WallLJComponent(j)+ConfigurationalEnergy_ESComponent(j,j)+ConfigurationalEnergy_LJComponent(j,j)
        !Many-type
        do jj=j+1,MoleculeKinds
            EnergyTotalAverageCurrentOutput1D=EnergyTotalAverageCurrentOutput1D+ConfigurationalEnergy_ESComponent(jj,j)+ConfigurationalEnergy_LJComponent(jj,j)
            !Sum of typewise electro static energy /= total electro static energy,
            !as there are multiple counting on single-body interaction.
            !To avoid this, add single-body term to electro static interaction after contribute to total
            ConfigurationalEnergy_ESComponent(jj,j)=ConfigurationalEnergy_ESComponent(jj,j)+NExist(j)*ESSBIndependent(j)+NExist(jj)*ESSBIndependent(jj)
        end do
    end do
    !Take average
    EnergyTotalAverageCurrentOutput1D=EnergyTotalAverageCurrentOutput1D/NTotal
    !Output
    write(*,'(A32,F15.6,A6)')'Average configurational energy:',EnergyTotalAverageCurrentOutput1D/kJmolInAU,'kJ/mol'
    write(*,'(A25)')'Single-body interaction:'
    do j=1,MoleculeKinds
        if(NExist(j)>0) then
            write(*,'(A10,I2,A20,F15.6,A6)')'Type:',j,'Wall Lennard-Jones:',ConfigurationalEnergy_WallLJComponent(j)/NExist(j)/kJmolInAU,'kJ/mol'
        end if
    end do
    write(*,'(A23)')'Many-body interaction:'
    do j=1,MoleculeKinds
        if(NExist(j)>0) then
            do jj=1,j-1
                if(NExist(jj)>0) then
                    write(*,'(A10,I2,A3,I2,A16,F15.6,A6)')'Type:',j,'by',jj,'Electro static:',&
                        ConfigurationalEnergy_ESComponent(j,jj)/NExist(j)/kJmolInAU,'kJ/mol'
                    write(*,'(A10,I2,A3,I2,A16,F15.6,A6)')'Type:',j,'by',jj,'Lennard-Jones:',&
                        ConfigurationalEnergy_LJComponent(j,jj)/NExist(j)/kJmolInAU,'kJ/mol'
                end if
            end do
            write(*,'(A10,I2,A3,I2,A16,F15.6,A6)')'Type:',j,'by',j,'Electro static:',&
                ConfigurationalEnergy_ESComponent(j,j)/NExist(j)/kJmolInAU,'kJ/mol'
            write(*,'(A10,I2,A3,I2,A16,F15.6,A6)')'Type:',j,'by',j,'Lennard-Jones:',&
                ConfigurationalEnergy_LJComponent(j,j)/NExist(j)/kJmolInAU,'kJ/mol'
            do jj=j+1,MoleculeKinds
                if(NExist(jj)>0) then
                    write(*,'(A10,I2,A3,I2,A16,F15.6,A6)')'Type:',j,'by',jj,'Electro static:',&
                        ConfigurationalEnergy_ESComponent(jj,j)/NExist(j)/kJmolInAU,'kJ/mol'
                    write(*,'(A10,I2,A3,I2,A16,F15.6,A6)')'Type:',j,'by',jj,'Lennard-Jones:',&
                        ConfigurationalEnergy_LJComponent(jj,j)/NExist(j)/kJmolInAU,'kJ/mol'
                end if
            end do
        end if
    end do
end function EnergyTotalAverageCurrentOutput1D
!For 2D
function EnergyTotalAverageCurrentOutput2D()
    real*8::EnergyTotalAverageCurrentOutput2D
    integer::i,j,ii,jj
    !Initialization: set the counters to zero
    EnergyTotalAverageCurrentOutput2D=0d0
    ConfigurationalEnergy_ESSBDependent2DComponent=0d0
    ConfigurationalEnergy_ESComponent=0d0
    ConfigurationalEnergy_WallLJComponent=0d0
    ConfigurationalEnergy_LJComponent=0d0
    !Sum up total
    do j=1,MoleculeKinds
        do i=1,NExist(j)
            !Single-type
                !Single-body
                ConfigurationalEnergy_ESSBDependent2DComponent(j)=ConfigurationalEnergy_ESSBDependent2DComponent(j)+mtp(j).mol(i).ESSBDependent2D
                ConfigurationalEnergy_WallLJComponent(j)=ConfigurationalEnergy_WallLJComponent(j)+mtp(j).mol(i).WallLJ
                !Many-body
                do ii=i+1,NExist(j)
                    ConfigurationalEnergy_ESComponent(j,j)=ConfigurationalEnergy_ESComponent(j,j)+intertp(j,j).molab(ii,i).ESMB
                    ConfigurationalEnergy_LJComponent(j,j)=ConfigurationalEnergy_LJComponent(j,j)+intertp(j,j).molab(ii,i).LJ
                end do
            !Many-type
            !Many-body
            do jj=j+1,MoleculeKinds
                do ii=1,NExist(jj)
                    ConfigurationalEnergy_ESComponent(jj,j)=ConfigurationalEnergy_ESComponent(jj,j)+intertp(jj,j).molab(ii,i).ESMB
                    ConfigurationalEnergy_LJComponent(jj,j)=ConfigurationalEnergy_LJComponent(jj,j)+intertp(jj,j).molab(ii,i).LJ
                end do
            end do
        end do
        !Single-type
        !Add single-body term to electro static interaction
        ConfigurationalEnergy_ESComponent(j,j)=ConfigurationalEnergy_ESComponent(j,j)+NExist(j)*ESSBIndependent(j)+ConfigurationalEnergy_ESSBDependent2DComponent(j)
        EnergyTotalAverageCurrentOutput2D=EnergyTotalAverageCurrentOutput2D+ConfigurationalEnergy_WallLJComponent(j)+ConfigurationalEnergy_ESComponent(j,j)+ConfigurationalEnergy_LJComponent(j,j)
        !Many-type
        do jj=j+1,MoleculeKinds
            EnergyTotalAverageCurrentOutput2D=EnergyTotalAverageCurrentOutput2D+ConfigurationalEnergy_ESComponent(jj,j)+ConfigurationalEnergy_LJComponent(jj,j)
            !Sum of typewise electro static energy /= total electro static energy,
            !as there are multiple counting on single-body interaction.
            !To avoid this, add single-body term to electro static interaction after contribute to total
            ConfigurationalEnergy_ESComponent(jj,j)=ConfigurationalEnergy_ESComponent(jj,j)+NExist(j)*ESSBIndependent(j)+ConfigurationalEnergy_ESSBDependent2DComponent(j)&
                +NExist(jj)*ESSBIndependent(jj)+ConfigurationalEnergy_ESSBDependent2DComponent(jj)
        end do
    end do
    !Take average
    EnergyTotalAverageCurrentOutput2D=EnergyTotalAverageCurrentOutput2D/NTotal
    !Output
    write(*,'(A32,F15.6,A6)')'Average configurational energy:',EnergyTotalAverageCurrentOutput2D/kJmolInAU,'kJ/mol'
    write(*,'(A25)')'Single-body interaction:'
    do j=1,MoleculeKinds
        if(NExist(j)>0) then
            write(*,'(A10,I2,A20,F15.6,A6)')'Type:',j,'Wall Lennard-Jones:',ConfigurationalEnergy_WallLJComponent(j)/NExist(j)/kJmolInAU,'kJ/mol'
        end if
    end do
    write(*,'(A23)')'Many-body interaction:'
    do j=1,MoleculeKinds
        if(NExist(j)>0) then
            do jj=1,j-1
                if(NExist(jj)>0) then
                    write(*,'(A10,I2,A3,I2,A16,F15.6,A6)')'Type:',j,'by',jj,'Electro static:',&
                        ConfigurationalEnergy_ESComponent(j,jj)/NExist(j)/kJmolInAU,'kJ/mol'
                    write(*,'(A10,I2,A3,I2,A16,F15.6,A6)')'Type:',j,'by',jj,'Lennard-Jones:',&
                        ConfigurationalEnergy_LJComponent(j,jj)/NExist(j)/kJmolInAU,'kJ/mol'
                end if
            end do
            write(*,'(A10,I2,A3,I2,A16,F15.6,A6)')'Type:',j,'by',j,'Electro static:',&
                ConfigurationalEnergy_ESComponent(j,j)/NExist(j)/kJmolInAU,'kJ/mol'
            write(*,'(A10,I2,A3,I2,A16,F15.6,A6)')'Type:',j,'by',j,'Lennard-Jones:',&
                ConfigurationalEnergy_LJComponent(j,j)/NExist(j)/kJmolInAU,'kJ/mol'
            do jj=j+1,MoleculeKinds
                if(NExist(jj)>0) then
                    write(*,'(A10,I2,A3,I2,A16,F15.6,A6)')'Type:',j,'by',jj,'Electro static:',&
                        ConfigurationalEnergy_ESComponent(jj,j)/NExist(j)/kJmolInAU,'kJ/mol'
                    write(*,'(A10,I2,A3,I2,A16,F15.6,A6)')'Type:',j,'by',jj,'Lennard-Jones:',&
                        ConfigurationalEnergy_LJComponent(jj,j)/NExist(j)/kJmolInAU,'kJ/mol'
                end if
            end do
        end if
    end do
end function EnergyTotalAverageCurrentOutput2D
!For 3D
function EnergyTotalAverageCurrentOutput3D()
    real*8::EnergyTotalAverageCurrentOutput3D
    integer::i,j,ii,jj
    !Initialization: set the counters to zero
    EnergyTotalAverageCurrentOutput3D=0d0
    ConfigurationalEnergy_ESComponent=0d0
    ConfigurationalEnergy_WallLJComponent=0d0
    ConfigurationalEnergy_LJComponent=0d0
    !Sum up total
    do j=1,MoleculeKinds
        do i=1,NExist(j)
            !Single-type
            !Currently no position dependent single-body interaction for 3D
            !Many-body
            do ii=i+1,NExist(j)
                ConfigurationalEnergy_ESComponent(j,j)=ConfigurationalEnergy_ESComponent(j,j)+intertp(j,j).molab(ii,i).ESMB
                ConfigurationalEnergy_LJComponent(j,j)=ConfigurationalEnergy_LJComponent(j,j)+intertp(j,j).molab(ii,i).LJ
            end do
            !Many-type
            !Currently no position dependent single-body interaction for 3D
            !Many-body
            do jj=j+1,MoleculeKinds
                do ii=1,NExist(jj)
                    ConfigurationalEnergy_ESComponent(jj,j)=ConfigurationalEnergy_ESComponent(jj,j)+intertp(jj,j).molab(ii,i).ESMB
                    ConfigurationalEnergy_LJComponent(jj,j)=ConfigurationalEnergy_LJComponent(jj,j)+intertp(jj,j).molab(ii,i).LJ
                end do
            end do
        end do
        !Single-type
        !Add single-body term to electro static interaction
        ConfigurationalEnergy_ESComponent(j,j)=ConfigurationalEnergy_ESComponent(j,j)+NExist(j)*ESSBIndependent(j)
        EnergyTotalAverageCurrentOutput3D=EnergyTotalAverageCurrentOutput3D+ConfigurationalEnergy_ESComponent(j,j)+ConfigurationalEnergy_LJComponent(j,j)
        !Many-type
        do jj=j+1,MoleculeKinds
            EnergyTotalAverageCurrentOutput3D=EnergyTotalAverageCurrentOutput3D+ConfigurationalEnergy_ESComponent(jj,j)+ConfigurationalEnergy_LJComponent(jj,j)
            !Sum of typewise electro static energy /= total electro static energy,
            !as there are multiple counting on single-body interaction.
            !To avoid this, add single-body term to electro static interaction after contribute to total
            ConfigurationalEnergy_ESComponent(jj,j)=ConfigurationalEnergy_ESComponent(jj,j)+NExist(j)*ESSBIndependent(j)+NExist(jj)*ESSBIndependent(jj)
        end do
    end do
    !Take average
    EnergyTotalAverageCurrentOutput3D=EnergyTotalAverageCurrentOutput3D/NTotal
    !Output
    write(*,'(A32,F15.6,A6)')'Average configurational energy:',EnergyTotalAverageCurrentOutput3D/kJmolInAU,'kJ/mol'
    write(*,'(A23)')'Many-body interaction:'
    do j=1,MoleculeKinds
        if(NExist(j)>0) then
            do jj=1,j-1
                if(NExist(jj)>0) then
                    write(*,'(A10,I2,A3,I2,A16,F15.6,A6)')'Type:',j,'by',jj,'Electro static:',&
                        ConfigurationalEnergy_ESComponent(j,jj)/NExist(j)/kJmolInAU,'kJ/mol'
                    write(*,'(A10,I2,A3,I2,A16,F15.6,A6)')'Type:',j,'by',jj,'Lennard-Jones:',&
                        ConfigurationalEnergy_LJComponent(j,jj)/NExist(j)/kJmolInAU,'kJ/mol'
                end if
            end do
            write(*,'(A10,I2,A3,I2,A16,F15.6,A6)')'Type:',j,'by',j,'Electro static:',&
                ConfigurationalEnergy_ESComponent(j,j)/NExist(j)/kJmolInAU,'kJ/mol'
            write(*,'(A10,I2,A3,I2,A16,F15.6,A6)')'Type:',j,'by',j,'Lennard-Jones:',&
                ConfigurationalEnergy_LJComponent(j,j)/NExist(j)/kJmolInAU,'kJ/mol'
            do jj=j+1,MoleculeKinds
                if(NExist(jj)>0) then
                    write(*,'(A10,I2,A3,I2,A16,F15.6,A6)')'Type:',j,'by',jj,'Electro static:',&
                        ConfigurationalEnergy_ESComponent(jj,j)/NExist(j)/kJmolInAU,'kJ/mol'
                    write(*,'(A10,I2,A3,I2,A16,F15.6,A6)')'Type:',j,'by',jj,'Lennard-Jones:',&
                        ConfigurationalEnergy_LJComponent(jj,j)/NExist(j)/kJmolInAU,'kJ/mol'
                end if
            end do
        end if
    end do
end function EnergyTotalAverageCurrentOutput3D

!---------- Updating the interaction status and computing the energy difference can be done in a same loop ----------
    !For relaxation
    !For 1D
    function DeltaRelax1D(tp,chosen)
        real*8::DeltaRelax1D
        integer,intent(in)::tp,chosen
        integer::i,j
        !Single-body
        !Currently WallLJ is the only 1D position dependent single-body term, so we do not refer to Total
        mtp(tp).mol(chosen).WallLJ=WallLJPotential1D(tp,mtp(tp).mol(chosen).Coordinate)
        DeltaRelax1D=mtp(tp).mol(chosen).WallLJ-mtpold(tp).mol(chosen).WallLJ
        !Many-body
        !Interaction between smaller type molecules and chosen molecule
        do j=1,tp-1
            do i=1,NExist(j)
                intertp(tp,j).molab(chosen,i).Total=NonbondedManyBody1D(tp,j,mtp(tp).mol(chosen).Coordinate,mtp(j).mol(i).Coordinate,&
                    intertp(tp,j).molab(chosen,i).ESMB,intertp(tp,j).molab(chosen,i).LJ)
                DeltaRelax1D=DeltaRelax1D+intertp(tp,j).molab(chosen,i).Total-intertpold(tp,j).molab(chosen,i).Total
            end do
        end do
        !Interaction between same type smaller order molecules and chosen molecule
        do i=1,chosen-1
            intertp(tp,tp).molab(chosen,i).Total=NonbondedManyBody1D(tp,tp,mtp(tp).mol(chosen).Coordinate,mtp(tp).mol(i).Coordinate,&
                intertp(tp,tp).molab(chosen,i).ESMB,intertp(tp,tp).molab(chosen,i).LJ)
            DeltaRelax1D=DeltaRelax1D+intertp(tp,tp).molab(chosen,i).Total-intertpold(tp,tp).molab(chosen,i).Total
        end do
        !Interaction between same type larger order molecules and chosen molecule
        do i=chosen+1,NExist(tp)
            intertp(tp,tp).molab(i,chosen).Total=NonbondedManyBody1D(tp,tp,mtp(tp).mol(i).Coordinate,mtp(tp).mol(chosen).Coordinate,&
                intertp(tp,tp).molab(i,chosen).ESMB,intertp(tp,tp).molab(i,chosen).LJ)
            DeltaRelax1D=DeltaRelax1D+intertp(tp,tp).molab(i,chosen).Total-intertpold(tp,tp).molab(i,chosen).Total
        end do
        !Interaction between larger type molecules and chosen molecule
        do j=tp+1,MoleculeKinds
            do i=1,NExist(j)
                intertp(j,tp).molab(i,chosen).Total=NonbondedManyBody1D(j,tp,mtp(j).mol(i).Coordinate,mtp(tp).mol(chosen).Coordinate,&
                    intertp(j,tp).molab(i,chosen).ESMB,intertp(j,tp).molab(i,chosen).LJ)
                DeltaRelax1D=DeltaRelax1D+intertp(j,tp).molab(i,chosen).Total-intertpold(j,tp).molab(i,chosen).Total
            end do
        end do
    end function DeltaRelax1D
    !For 2D
    function DeltaRelax2D(tp,chosen)
        real*8::DeltaRelax2D
        integer,intent(in)::tp,chosen
        integer::i,j
        !Single-body
        mtp(tp).mol(chosen).ESSBDependent2D=EwaldSummationSingleBodyPositionDependent(tp,mtp(tp).mol(chosen).Coordinate(3,:))
        mtp(tp).mol(chosen).WallLJ=WallLJPotential2D(tp,mtp(tp).mol(chosen).Coordinate)
        mtp(tp).mol(chosen).Total=mtp(tp).mol(chosen).ESSBDependent2D+mtp(tp).mol(chosen).WallLJ
        DeltaRelax2D=mtp(tp).mol(chosen).Total-mtpold(tp).mol(chosen).Total
        !Many-body
        !Interaction between smaller type molecules and chosen molecule
        do j=1,tp-1
            do i=1,NExist(j)
                intertp(tp,j).molab(chosen,i).Total=NonbondedManyBody2D(tp,j,mtp(tp).mol(chosen).Coordinate,mtp(j).mol(i).Coordinate,&
                    intertp(tp,j).molab(chosen,i).ESMB,intertp(tp,j).molab(chosen,i).LJ)
                DeltaRelax2D=DeltaRelax2D+intertp(tp,j).molab(chosen,i).Total-intertpold(tp,j).molab(chosen,i).Total
            end do
        end do
        !Interaction between same type smaller order molecules and chosen molecule
        do i=1,chosen-1
            intertp(tp,tp).molab(chosen,i).Total=NonbondedManyBody2D(tp,tp,mtp(tp).mol(chosen).Coordinate,mtp(tp).mol(i).Coordinate,&
                intertp(tp,tp).molab(chosen,i).ESMB,intertp(tp,tp).molab(chosen,i).LJ)
            DeltaRelax2D=DeltaRelax2D+intertp(tp,tp).molab(chosen,i).Total-intertpold(tp,tp).molab(chosen,i).Total
        end do
        !Interaction between same type larger order molecules and chosen molecule
        do i=chosen+1,NExist(tp)
            intertp(tp,tp).molab(i,chosen).Total=NonbondedManyBody2D(tp,tp,mtp(tp).mol(i).Coordinate,mtp(tp).mol(chosen).Coordinate,&
                intertp(tp,tp).molab(i,chosen).ESMB,intertp(tp,tp).molab(i,chosen).LJ)
            DeltaRelax2D=DeltaRelax2D+intertp(tp,tp).molab(i,chosen).Total-intertpold(tp,tp).molab(i,chosen).Total
        end do
        !Interaction between larger type molecules and chosen molecule
        do j=tp+1,MoleculeKinds
            do i=1,NExist(j)
                intertp(j,tp).molab(i,chosen).Total=NonbondedManyBody2D(j,tp,mtp(j).mol(i).Coordinate,mtp(tp).mol(chosen).Coordinate,&
                    intertp(j,tp).molab(i,chosen).ESMB,intertp(j,tp).molab(i,chosen).LJ)
                DeltaRelax2D=DeltaRelax2D+intertp(j,tp).molab(i,chosen).Total-intertpold(j,tp).molab(i,chosen).Total
            end do
        end do
    end function DeltaRelax2D
    !For 3D
    function DeltaRelax3D(tp,chosen)
        real*8::DeltaRelax3D
        integer,intent(in)::tp,chosen
        integer::i,j
        !Currently no position dependent single-body interaction for 3D
        !Initialization: set the counter to zero
        DeltaRelax3D=0d0
        !Many-body
        !Interaction between smaller type molecules and chosen molecule
        do j=1,tp-1
            do i=1,NExist(j)
                intertp(tp,j).molab(chosen,i).Total=NonbondedManyBody3D(tp,j,mtp(tp).mol(chosen).Coordinate,mtp(j).mol(i).Coordinate,&
                    intertp(tp,j).molab(chosen,i).ESMB,intertp(tp,j).molab(chosen,i).LJ)
                DeltaRelax3D=DeltaRelax3D+intertp(tp,j).molab(chosen,i).Total-intertpold(tp,j).molab(chosen,i).Total
            end do
        end do
        !Interaction between same type smaller order molecules and chosen molecule
        do i=1,chosen-1
            intertp(tp,tp).molab(chosen,i).Total=NonbondedManyBody3D(tp,tp,mtp(tp).mol(chosen).Coordinate,mtp(tp).mol(i).Coordinate,&
                intertp(tp,tp).molab(chosen,i).ESMB,intertp(tp,tp).molab(chosen,i).LJ)
            DeltaRelax3D=DeltaRelax3D+intertp(tp,tp).molab(chosen,i).Total-intertpold(tp,tp).molab(chosen,i).Total
        end do
        !Interaction between same type larger order molecules and chosen molecule
        do i=chosen+1,NExist(tp)
            intertp(tp,tp).molab(i,chosen).Total=NonbondedManyBody3D(tp,tp,mtp(tp).mol(i).Coordinate,mtp(tp).mol(chosen).Coordinate,&
                intertp(tp,tp).molab(i,chosen).ESMB,intertp(tp,tp).molab(i,chosen).LJ)
            DeltaRelax3D=DeltaRelax3D+intertp(tp,tp).molab(i,chosen).Total-intertpold(tp,tp).molab(i,chosen).Total
        end do
        !Interaction between larger type molecules and chosen molecule
        do j=tp+1,MoleculeKinds
            do i=1,NExist(j)
                intertp(j,tp).molab(i,chosen).Total=NonbondedManyBody3D(j,tp,mtp(j).mol(i).Coordinate,mtp(tp).mol(chosen).Coordinate,&
                    intertp(j,tp).molab(i,chosen).ESMB,intertp(j,tp).molab(i,chosen).LJ)
                DeltaRelax3D=DeltaRelax3D+intertp(j,tp).molab(i,chosen).Total-intertpold(j,tp).molab(i,chosen).Total
            end do
        end do
    end function DeltaRelax3D
    
    !For volume change
    function DeltaVolumeChange()
        real*8::DeltaVolumeChange
        integer::i,j,ii,jj
        !Initialization: set the counter with the single-body position independent term
        DeltaVolumeChange=dot_product(NExist,ESSBIndependent-ESSBIndependentOld)
        do j=1,MoleculeKinds
            do i=1,NExist(j)
                !Single-body: currently no position dependent single-body interaction for 3D
                !Many-body
                do ii=i+1,NExist(j)
                    intertp(j,j).molab(ii,i).Total=NonbondedManyBody3D(j,j,mtp(j).mol(ii).Coordinate,mtp(j).mol(i).Coordinate,&
                        intertp(j,j).molab(ii,i).ESMB,intertp(j,j).molab(ii,i).LJ)
                    DeltaVolumeChange=DeltaVolumeChange+intertp(j,j).molab(ii,i).Total-intertpold(j,j).molab(ii,i).Total
                end do
                do jj=j+1,MoleculeKinds
                    do ii=1,NExist(jj)
                        intertp(jj,j).molab(ii,i).Total=NonbondedManyBody3D(jj,j,mtp(jj).mol(ii).Coordinate,mtp(j).mol(i).Coordinate,&
                            intertp(jj,j).molab(ii,i).ESMB,intertp(jj,j).molab(ii,i).LJ)
                        DeltaVolumeChange=DeltaVolumeChange+intertp(jj,j).molab(ii,i).Total-intertpold(jj,j).molab(ii,i).Total
                    end do
                end do
            end do
        end do
    end function DeltaVolumeChange
!---------------------------------------------------- End -----------------------------------------------------------

end module ConfigurationalEnergy