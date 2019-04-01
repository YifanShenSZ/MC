!Contains all required trial moves for available ensembles
module Trial
    use General
    use MCBasis
    use ConfigurationalEnergy
    implicit none

contains
!Perturb a translational or rotational degree of freedom
!The chosen-th type tp molecule is perturbed
!If direction<0.5, translational trial; else rotational trial
!For 1D
subroutine Relax1D(tp,chosen,direction,accept)
    integer,intent(out)::tp,chosen
    real*8,intent(out)::direction
    logical,intent(out)::accept
    integer::i,j
    real*8::rnd,del,step
    !Randomly choose a molecule
    call random_number(rnd)
    chosen=ceiling(rnd*NTotal)
    tp=MoleculeKinds
    do i=1,MoleculeKinds-1
        chosen=chosen-NExist(i)
        if(chosen<1) then
            chosen=chosen+NExist(i)
            tp=i
            exit
        end if
    end do
    !Determine the trial to make
    call random_number(direction)
    !Translational trial
    if(direction<0.5d0) then
        !Perturb ith dimension
            call random_number(rnd)
            i=ceiling(rnd*3d0)
            TrialCount(i,tp)=TrialCount(i,tp)+1
        call random_number(rnd)
        step=TrialStep(i,tp)*(rnd*2d0-1d0)
        mtp(tp).mol(chosen).Coordinate(i,:)=mtp(tp).mol(chosen).Coordinate(i,:)+step
        !Periodic boundary condition
        !For periodic dimension, modification standard is the centre of mass
        !For non-periodic dimensions, check all atoms whether there exists one outside the wall
            if(i==1) then
                !The 1st direction displacement of the centre of mass
                rnd=mtp(tp).mol(chosen).Coordinate(1,1)-MoleculeInput(tp).RefConfig(1,1)
                if(rnd>HalfLength(1)) then
                    mtp(tp).mol(chosen).Coordinate(1,:)=mtp(tp).mol(chosen).Coordinate(1,:)-BoxSize(1)
                else if(rnd<-HalfLength(1)) then
                    mtp(tp).mol(chosen).Coordinate(1,:)=mtp(tp).mol(chosen).Coordinate(1,:)+BoxSize(1)
                end if
            else
                !xy plane is confined by the cylindrical wall
                do j=1,MoleculeInput(tp).NAtoms
                    if(dot_product(mtp(tp).mol(chosen).Coordinate(2:3,j),mtp(tp).mol(chosen).Coordinate(2:3,j))>radiussq) then
                        mtp(tp).mol(chosen).Coordinate(i,:)=mtp(tp).mol(chosen).Coordinate(i,:)-step
                        accept=0
                        return
                    end if
                end do
            end if
    !Rotational trial
    else
        !Form the unit quaternion
            do while(1)
                call random_number(QuaternionSupporter1)
                QuaternionSupporter1=QuaternionSupporter1*2d0-1d0
                if(dot_product(QuaternionSupporter1,QuaternionSupporter1)<1d0) then
                    exit
                end if
            end do
            do while(1)
                call random_number(QuaternionSupporter2)
                QuaternionSupporter2=QuaternionSupporter2*2d0-1d0
                if(dot_product(QuaternionSupporter2,QuaternionSupporter2)<1d0) then
                    exit
                end if
            end do
            RotationQuaternion(1:2)=QuaternionSupporter1
            RotationQuaternion(3:4)=QuaternionSupporter2&
                *Sqrt((1d0-dot_product(QuaternionSupporter1,QuaternionSupporter1))&
                /dot_product(QuaternionSupporter2,QuaternionSupporter2))
        !The position of the centre of mass
        DisplacementVector=mtp(tp).mol(chosen).Coordinate(:,1)-MoleculeInput(tp).RefConfig(:,1)
        do i=1,MoleculeInput(tp).NAtoms
            !Transform to the centre of mass frame
            mtp(tp).mol(chosen).Coordinate(:,i)=mtp(tp).mol(chosen).Coordinate(:,i)-DisplacementVector
            !Perform the rotation
            mtp(tp).mol(chosen).Coordinate(:,i)=rotate(RotationQuaternion,mtp(tp).mol(chosen).Coordinate(:,i))
            !Transform back to the simulation box frame
            mtp(tp).mol(chosen).Coordinate(:,i)=mtp(tp).mol(chosen).Coordinate(:,i)+DisplacementVector
            !xy plane is confined by the cylindrical wall, check if the atom is outside
            if(dot_product(mtp(tp).mol(chosen).Coordinate(2:3,i),mtp(tp).mol(chosen).Coordinate(2:3,i))>radiussq) then
                mtp(tp).mol(chosen).Coordinate=mtpold(tp).mol(chosen).Coordinate
                accept=0
                return
            end if
        end do
    end if
    !Energy difference
    del=DeltaRelax1D(tp,chosen)
    !Accept or reject this trial
    accept=del<ImpossibleDel
    if(accept.and.del>0d0) then
        call random_number(rnd)
        accept=rnd<exp(-del/temperature)
    end if
    if(accept) then
        !old=new
        !Single-body
        mtpold(tp).mol(chosen)=mtp(tp).mol(chosen)
        mtpold(tp).mol(chosen).WallLJ=mtp(tp).mol(chosen).WallLJ
        !Many-body
        do j=1,tp-1
            forall(i=1:NExist(j))
                intertpold(tp,j).molab(chosen,i)=intertp(tp,j).molab(chosen,i)
            end forall
        end do
        forall(i=1:chosen-1)
            intertpold(tp,tp).molab(chosen,i)=intertp(tp,tp).molab(chosen,i)
        end forall
        forall(i=chosen+1:NExist(tp))
            intertpold(tp,tp).molab(i,chosen)=intertp(tp,tp).molab(i,chosen)
        end forall
        do j=tp+1,MoleculeKinds
            forall(i=1:NExist(j))
                intertpold(j,tp).molab(i,chosen)=intertp(j,tp).molab(i,chosen)
            end forall
        end do
        !Count translational trial success
        if(direction<0.5d0) then
            TrialSuccessCount(i,tp)=TrialSuccessCount(i,tp)+1
        end if
    else
        !new=old
        !Single-body
        mtp(tp).mol(chosen)=mtpold(tp).mol(chosen)
        !Many-body
        do j=1,tp-1
            forall(i=1:NExist(j))
                intertp(tp,j).molab(chosen,i)=intertpold(tp,j).molab(chosen,i)
            end forall
        end do
        forall(i=1:chosen-1)
            intertp(tp,tp).molab(chosen,i)=intertpold(tp,tp).molab(chosen,i)
        end forall
        forall(i=chosen+1:NExist(tp))
            intertp(tp,tp).molab(i,chosen)=intertpold(tp,tp).molab(i,chosen)
        end forall
        do j=tp+1,MoleculeKinds
            forall(i=1:NExist(j))
                intertp(j,tp).molab(i,chosen)=intertpold(j,tp).molab(i,chosen)
            end forall
        end do
    end if
end subroutine Relax1D
!For 2D
subroutine Relax2D(tp,chosen,direction,accept)
    integer,intent(out)::tp,chosen
    real*8,intent(out)::direction
    logical,intent(out)::accept
    integer::i,j
    real*8::rnd,del,step
    !Randomly choose a molecule
    call random_number(rnd)
    chosen=ceiling(rnd*NTotal)
    tp=MoleculeKinds
    do i=1,MoleculeKinds-1
        chosen=chosen-NExist(i)
        if(chosen<1) then
            chosen=chosen+NExist(i)
            tp=i
            exit
        end if
    end do
    !Determine the trial to make
    call random_number(direction)
    !Translational trial
    if(direction<0.5d0) then
        !Perturb ith dimension
            call random_number(rnd)
            i=ceiling(rnd*3d0)
            TrialCount(i,tp)=TrialCount(i,tp)+1
        call random_number(rnd)
        step=TrialStep(i,tp)*(rnd*2d0-1d0)
        mtp(tp).mol(chosen).Coordinate(i,:)=mtp(tp).mol(chosen).Coordinate(i,:)+step
        !Periodic boundary condition
        !For non-periodic dimension, check all atoms whether there exists one outside the wall
        !For periodic dimensions, modification standard is the centre of mass
        if(i==3) then
            !z direction is confined by the slit wall
            do j=1,MoleculeInput(tp).NAtoms
                if(abs(mtp(tp).mol(chosen).Coordinate(3,j))>HalfLength(3)) then
                    mtp(tp).mol(chosen).Coordinate(3,:)=mtp(tp).mol(chosen).Coordinate(3,:)-step
                    accept=0
                    return
                end if
            end do
        else
            !The i-th direction displacement of the centre of mass
            rnd=mtp(tp).mol(chosen).Coordinate(i,1)-MoleculeInput(tp).RefConfig(i,1)
            if(rnd>HalfLength(i)) then
                mtp(tp).mol(chosen).Coordinate(i,:)=mtp(tp).mol(chosen).Coordinate(i,:)-BoxSize(i)
            else if(rnd<-HalfLength(i)) then
                mtp(tp).mol(chosen).Coordinate(i,:)=mtp(tp).mol(chosen).Coordinate(i,:)+BoxSize(i)
            end if
        end if
    !Rotational trial
    else
        !Form the unit quaternion
            do while(1)
                call random_number(QuaternionSupporter1)
                QuaternionSupporter1=QuaternionSupporter1*2d0-1d0
                if(dot_product(QuaternionSupporter1,QuaternionSupporter1)<1d0) then
                    exit
                end if
            end do
            do while(1)
                call random_number(QuaternionSupporter2)
                QuaternionSupporter2=QuaternionSupporter2*2d0-1d0
                if(dot_product(QuaternionSupporter2,QuaternionSupporter2)<1d0) then
                    exit
                end if
            end do
            RotationQuaternion(1:2)=QuaternionSupporter1
            RotationQuaternion(3:4)=QuaternionSupporter2&
                *Sqrt((1d0-dot_product(QuaternionSupporter1,QuaternionSupporter1))&
                /dot_product(QuaternionSupporter2,QuaternionSupporter2))
        !The position of the centre of mass
        DisplacementVector=mtp(tp).mol(chosen).Coordinate(:,1)-MoleculeInput(tp).RefConfig(:,1)
        do i=1,MoleculeInput(tp).NAtoms
            !Transform to the centre of mass frame
            mtp(tp).mol(chosen).Coordinate(:,i)=mtp(tp).mol(chosen).Coordinate(:,i)-DisplacementVector
            !Perform the rotation
            mtp(tp).mol(chosen).Coordinate(:,i)=rotate(RotationQuaternion,mtp(tp).mol(chosen).Coordinate(:,i))
            !Transform back to the simulation box frame
            mtp(tp).mol(chosen).Coordinate(:,i)=mtp(tp).mol(chosen).Coordinate(:,i)+DisplacementVector
            !z direction is confined by the slit wall, check if the atom is outside
            if(abs(mtp(tp).mol(chosen).Coordinate(3,i))>HalfLength(3)) then
                mtp(tp).mol(chosen).Coordinate=mtpold(tp).mol(chosen).Coordinate
                accept=0
                return
            end if
        end do
    end if
    !Energy difference
    del=DeltaRelax2D(tp,chosen)
    !Accept or reject this trial
    accept=del<ImpossibleDel
    if(accept.and.del>0d0) then
        call random_number(rnd)
        accept=rnd<exp(-del/temperature)
    end if
    if(accept) then
        !old=new
        !Single-body
        mtpold(tp).mol(chosen)=mtp(tp).mol(chosen)
        !Many-body
        do j=1,tp-1
            forall(i=1:NExist(j))
                intertpold(tp,j).molab(chosen,i)=intertp(tp,j).molab(chosen,i)
            end forall
        end do
        forall(i=1:chosen-1)
            intertpold(tp,tp).molab(chosen,i)=intertp(tp,tp).molab(chosen,i)
        end forall
        forall(i=chosen+1:NExist(tp))
            intertpold(tp,tp).molab(i,chosen)=intertp(tp,tp).molab(i,chosen)
        end forall
        do j=tp+1,MoleculeKinds
            forall(i=1:NExist(j))
                intertpold(j,tp).molab(i,chosen)=intertp(j,tp).molab(i,chosen)
            end forall
        end do
        !Count translational trial success
        if(direction<0.5d0) then
            TrialSuccessCount(i,tp)=TrialSuccessCount(i,tp)+1
        end if
    else
        !new=old
        !Single-body
        mtp(tp).mol(chosen)=mtpold(tp).mol(chosen)
        !Many-body
        do j=1,tp-1
            forall(i=1:NExist(j))
                intertp(tp,j).molab(chosen,i)=intertpold(tp,j).molab(chosen,i)
            end forall
        end do
        forall(i=1:chosen-1)
            intertp(tp,tp).molab(chosen,i)=intertpold(tp,tp).molab(chosen,i)
        end forall
        forall(i=chosen+1:NExist(tp))
            intertp(tp,tp).molab(i,chosen)=intertpold(tp,tp).molab(i,chosen)
        end forall
        do j=tp+1,MoleculeKinds
            forall(i=1:NExist(j))
                intertp(j,tp).molab(i,chosen)=intertpold(j,tp).molab(i,chosen)
            end forall
        end do
    end if
end subroutine Relax2D
!For 3D
subroutine Relax3D(tp,chosen,direction,accept)
    integer,intent(out)::tp,chosen
    real*8,intent(out)::direction
    logical,intent(out)::accept
    integer::i,j
    real*8::rnd,del,step
    !Randomly choose a molecule
    call random_number(rnd)
    chosen=ceiling(rnd*NTotal)
    tp=MoleculeKinds
    do i=1,MoleculeKinds-1
        chosen=chosen-NExist(i)
        if(chosen<1) then
            chosen=chosen+NExist(i)
            tp=i
            exit
        end if
    end do
    !Determine the trial to make
    call random_number(direction)
    !Translational trial
    if(direction<0.5d0) then
        !Perturb ith dimension
            call random_number(rnd)
            i=ceiling(rnd*3d0)
            TrialCount(i,tp)=TrialCount(i,tp)+1
        call random_number(rnd)
        step=TrialStep(i,tp)*(rnd*2d0-1d0)
        mtp(tp).mol(chosen).Coordinate(i,:)=mtp(tp).mol(chosen).Coordinate(i,:)+step
        !Periodic boundary condition
        !For periodic dimensions, modification standard is the centre of mass
        !The i-th direction displacement of the centre of mass
        rnd=mtp(tp).mol(chosen).Coordinate(i,1)-MoleculeInput(tp).RefConfig(i,1)
        if(rnd>HalfLength(i)) then
            mtp(tp).mol(chosen).Coordinate(i,:)=mtp(tp).mol(chosen).Coordinate(i,:)-BoxSize(i)
        else if(rnd<-HalfLength(i)) then
            mtp(tp).mol(chosen).Coordinate(i,:)=mtp(tp).mol(chosen).Coordinate(i,:)+BoxSize(i)
        end if
    !Rotational trial
    else
        !Form the unit quaternion
            do while(1)
                call random_number(QuaternionSupporter1)
                QuaternionSupporter1=QuaternionSupporter1*2d0-1d0
                if(dot_product(QuaternionSupporter1,QuaternionSupporter1)<1d0) then
                    exit
                end if
            end do
            do while(1)
                call random_number(QuaternionSupporter2)
                QuaternionSupporter2=QuaternionSupporter2*2d0-1d0
                if(dot_product(QuaternionSupporter2,QuaternionSupporter2)<1d0) then
                    exit
                end if
            end do
            RotationQuaternion(1:2)=QuaternionSupporter1
            RotationQuaternion(3:4)=QuaternionSupporter2&
                *Sqrt((1d0-dot_product(QuaternionSupporter1,QuaternionSupporter1))&
                /dot_product(QuaternionSupporter2,QuaternionSupporter2))
        !The position of the centre of mass
        DisplacementVector=mtp(tp).mol(chosen).Coordinate(:,1)-MoleculeInput(tp).RefConfig(:,1)
        do i=1,MoleculeInput(tp).NAtoms
            !Transform to the centre of mass frame
            mtp(tp).mol(chosen).Coordinate(:,i)=mtp(tp).mol(chosen).Coordinate(:,i)-DisplacementVector
            !Perform the rotation
            mtp(tp).mol(chosen).Coordinate(:,i)=rotate(RotationQuaternion,mtp(tp).mol(chosen).Coordinate(:,i))
            !Transform back to the simulation box frame
            mtp(tp).mol(chosen).Coordinate(:,i)=mtp(tp).mol(chosen).Coordinate(:,i)+DisplacementVector
        end do
    end if
    !Energy difference
    del=DeltaRelax3D(tp,chosen)
    !Accept or reject this trial
    accept=del<ImpossibleDel
    if(accept.and.del>0d0) then
        call random_number(rnd)
        accept=rnd<exp(-del/temperature)
    end if
    if(accept) then
        !old=new
        !Single-body
        mtpold(tp).mol(chosen)=mtp(tp).mol(chosen)
        !Many-body
        do j=1,tp-1
            forall(i=1:NExist(j))
                intertpold(tp,j).molab(chosen,i)=intertp(tp,j).molab(chosen,i)
            end forall
        end do
        forall(i=1:chosen-1)
            intertpold(tp,tp).molab(chosen,i)=intertp(tp,tp).molab(chosen,i)
        end forall
        forall(i=chosen+1:NExist(tp))
            intertpold(tp,tp).molab(i,chosen)=intertp(tp,tp).molab(i,chosen)
        end forall
        do j=tp+1,MoleculeKinds
            forall(i=1:NExist(j))
                intertpold(j,tp).molab(i,chosen)=intertp(j,tp).molab(i,chosen)
            end forall
        end do
        !Count translational trial success
        if(direction<0.5d0) then
            TrialSuccessCount(i,tp)=TrialSuccessCount(i,tp)+1
        end if
    else
        !new=old
        !Single-body
        mtp(tp).mol(chosen)=mtpold(tp).mol(chosen)
        !Many-body
        do j=1,tp-1
            forall(i=1:NExist(j))
                intertp(tp,j).molab(chosen,i)=intertpold(tp,j).molab(chosen,i)
            end forall
        end do
        forall(i=1:chosen-1)
            intertp(tp,tp).molab(chosen,i)=intertpold(tp,tp).molab(chosen,i)
        end forall
        forall(i=chosen+1:NExist(tp))
            intertp(tp,tp).molab(i,chosen)=intertpold(tp,tp).molab(i,chosen)
        end forall
        do j=tp+1,MoleculeKinds
            forall(i=1:NExist(j))
                intertp(j,tp).molab(i,chosen)=intertpold(j,tp).molab(i,chosen)
            end forall
        end do
    end if
end subroutine Relax3D

!Perturb the box size, for P type ensemble only
!Specialized for 3 dimensionally periodic system (only such system could undergo P type ensemble)
subroutine VolumeChange(accept)
    logical,intent(out)::accept
    integer::i,j,ii,jj
    real*8::rnd,del,rate,lengthrate
    !lnV trial
    call random_number(rnd)
    lengthrate=1d0+(rnd*2d0-1d0)*(Norm2(TrialStep(:,1))/Norm2(HalfLength))
    rate=lengthrate*lengthrate*lengthrate
    !Update box size
    HalfLength=HalfLength*lengthrate
    BoxSize=BoxSize*lengthrate
    volume=volume*rate
    !Reinitialize Ewald summation
    wavenumber=wavenumber/lengthrate
    pim2dV=pim2dV/rate
    pim4dV=pim4dV/rate
    call EwaldSummationSingleBodyPositionIndependent3D()
    call EwaldSummationDumpFactor3D()
    !Update coordinates
    do j=1,MoleculeKinds
        do i=1,NExist(j)
            !Update coordinates
            DisplacementVector=(mtp(j).mol(i).Coordinate(:,1)-MoleculeInput(j).RefConfig(:,1))*(lengthrate-1d0)
            forall(ii=1:MoleculeInput(j).NAtoms)
                mtp(j).mol(i).Coordinate(:,ii)=mtp(j).mol(i).Coordinate(:,ii)+DisplacementVector
            end forall
        end do
    end do
    !Energy difference and volume change special term
    del=DeltaVolumeChange()+pressure*(1d0-1d0/rate)*volume-temperature*(NTotal+1)*Log(rate)
    !Accept or reject this trial
    accept=del<ImpossibleDel
    if(accept.and.del>0d0) then
        call random_number(rnd)
        accept=rnd<exp(-del/temperature)
    end if
    if(accept) then
        !This program does not have such isobaric ensemble that requires translational partition function
        !old=new
        !Position independent term
        ESSBIndependentOld=ESSBIndependent
        do j=1,MoleculeKinds
            do i=1,NExist(j)
                !Single-body
                mtpold(j).mol(i)=mtp(j).mol(i)
                !Many-body
                forall(ii=i+1:NExist(j))
                    intertpold(j,j).molab(ii,i)=intertp(j,j).molab(ii,i)
                end forall
                do jj=j+1,MoleculeKinds
                    forall(ii=1:NExist(jj))
                        intertpold(jj,j).molab(ii,i)=intertp(jj,j).molab(ii,i)
                    end forall
                end do
            end do
        end do
    else
        !new=old
        !Box size
        HalfLength=HalfLength/lengthrate
        BoxSize=BoxSize/lengthrate
        volume=volume/rate
        !Volume dependent parameters in Ewald summation
        wavenumber=wavenumber*lengthrate
        pim2dV=pim2dV*rate
        pim4dV=pim4dV*rate
        chi=chi*lengthrate
        chisq=chi*chi
        minus4chisq=-4d0*chisq
        !Position independent term
        ESSBIndependent=ESSBIndependentOld
        do j=1,MoleculeKinds
            do i=1,NExist(j)
                !Single-body
                mtp(j).mol(i)=mtpold(j).mol(i)
                !Many-body
                forall(ii=i+1:NExist(j))
                    intertp(j,j).molab(ii,i)=intertpold(j,j).molab(ii,i)
                end forall
                do jj=j+1,MoleculeKinds
                    forall(ii=1:NExist(jj))
                        intertp(jj,j).molab(ii,i)=intertpold(jj,j).molab(ii,i)
                    end forall
                end do
            end do
        end do
    end if
end subroutine VolumeChange

!Delete or insert a molecule
!For deletion, the chosen-th type tp molecule is deleted
!For insertion, a type tp molecule is inserted
!If direction<0.5, deletion trial; else insertion trial
!For 1D
subroutine Insertion1D(tp,chosen,direction,accept)
    integer,intent(out)::tp,chosen
    real*8,intent(out)::direction
    logical,intent(out)::accept
    integer::i,j
    real*8::rnd,del
    !Randomly choose a type of molecule
    call random_number(rnd)
    tp=MoleculeKinds
    do i=1,MoleculeKinds-1
        if(rnd<ExChoose(i)) then
            tp=i
            exit
        end if
    end do
    !Determine the trial to make
    call random_number(direction)
    !Deletion
    if(direction<0.5d0) then
        !Check the availability of the molecule to be deleted
        if(NExist(tp)<1) then
            accept=0
            return
        end if
        !Randomly choose a type tp molecule
        call random_number(rnd)
        chosen=ceiling(rnd*NExist(tp))
        !Energy difference
        del=-EnergyChosenCurrent1D(tp,chosen)
        !Accept or reject this trial
        call random_number(rnd)
        accept=rnd<NExist(tp)/TranslationalPartitionFunction(tp)*exp(-(ChemicalPotential(tp)+del)/temperature)
    !Insertion
    else
        !Check the availability of workspace for the molecule to be inserted
        if(NExist(tp)==MaxPossibleN(tp)) then
            if(MaxPossibleNReachedWarning(tp)) then
                write(*,'(A33,I2,A17)')'Maximum possible number for type',tp,'molecule reached'
                write(*,'(A61)')'Please increase it before run grand canonical ensemble again'
                MaxPossibleNReachedWarning(tp)=0
            end if
            accept=0
            return
        end if
        chosen=NExist(tp)+1
        !Insert the new molecule randomly and uniformly
        call InsertMolecule1D(DisplacementVector,RotationQuaternion,.true.)
        mtp(tp).mol(chosen).Coordinate=struct(DisplacementVector,RotationQuaternion,&
            MoleculeInput(tp).RefConfig,MoleculeInput(tp).NAtoms)
        !xy plane is confined by the cylindrical wall, check if the atom is outside
        do i=1,MoleculeInput(tp).NAtoms
            if(dot_product(mtp(tp).mol(chosen).Coordinate(2:3,i),mtp(tp).mol(chosen).Coordinate(2:3,i))>radiussq) then
                accept=0
                return
            end if
        end do
        !Energy difference
        del=EnergyChosenUpdated1D(tp,chosen)
        !Accept or reject this trial
        call random_number(rnd)
        accept=rnd<TranslationalPartitionFunction(tp)/(NExist(tp)+1)*exp((ChemicalPotential(tp)-del)/temperature)
    end if
    !Nothing is changed, so nothing needs to be done if rejected
    if(accept) then
        if(direction<0.5d0) then
            NExist(tp)=NExist(tp)-1
            NTotal=NTotal-1
            !Delete the chosen molecule by flushing it with the last exist one, old=new by the way
            !Single-body
            mtp(tp).mol(chosen)=mtpold(tp).mol(NExistOld(tp))
            mtpold(tp).mol(chosen)=mtp(tp).mol(chosen)
            !Many-body
            do j=1,tp-1
                forall(i=1:NExistOld(j))
                    intertp(tp,j).molab(chosen,i)=intertpold(tp,j).molab(NExistOld(tp),i)
                    intertpold(tp,j).molab(chosen,i)=intertp(tp,j).molab(chosen,i)
                end forall
            end do
            forall(i=1:chosen-1)
                intertp(tp,tp).molab(chosen,i)=intertpold(tp,tp).molab(NExistOld(tp),i)
                intertpold(tp,tp).molab(chosen,i)=intertp(tp,tp).molab(chosen,i)
            end forall
            !molab(NExistOld(tp),chosen) is not needed as mol(NExistOld(tp)) is deleted
            forall(i=chosen+1:NExist(tp))
                intertp(tp,tp).molab(i,chosen)=intertpold(tp,tp).molab(NExistOld(tp),i)
                intertpold(tp,tp).molab(i,chosen)=intertp(tp,tp).molab(i,chosen)
            end forall
            do j=tp+1,MoleculeKinds
                forall(i=1:NExistOld(j))
                    intertp(j,tp).molab(i,chosen)=intertpold(j,tp).molab(i,NExistOld(tp))
                    intertpold(j,tp).molab(i,chosen)=intertp(j,tp).molab(i,chosen)
                end forall
            end do
            NExistOld(tp)=NExist(tp)
        else
            NExist(tp)=NExist(tp)+1
            NTotal=NTotal+1
            !old=new
            NExistOld(tp)=NExist(tp)
            !Single-body
            mtpold(tp).mol(chosen)=mtp(tp).mol(chosen)
            !Many-body
            do j=1,tp-1
                forall(i=1:NExist(j))
                    intertpold(tp,j).molab(chosen,i)=intertp(tp,j).molab(chosen,i)
                end forall
            end do
            !The chosen molecule is the last one so no need for the latter loop
            forall(i=1:chosen-1)
                intertpold(tp,tp).molab(chosen,i)=intertp(tp,tp).molab(chosen,i)
            end forall
            do j=tp+1,MoleculeKinds
                forall(i=1:NExist(j))
                    intertpold(j,tp).molab(i,chosen)=intertp(j,tp).molab(i,chosen)
                end forall
            end do
        end if
    end if
end subroutine Insertion1D
!For 2D
subroutine Insertion2D(tp,chosen,direction,accept)
    integer,intent(out)::tp,chosen
    real*8,intent(out)::direction
    logical,intent(out)::accept
    integer::i,j
    real*8::rnd,del
    !Randomly choose a type of molecule
    call random_number(rnd)
    tp=MoleculeKinds
    do i=1,MoleculeKinds-1
        if(rnd<ExChoose(i)) then
            tp=i
            exit
        end if
    end do
    !Determine the trial to make
    call random_number(direction)
    !Deletion
    if(direction<0.5d0) then
        !Check the availability of the molecule to be deleted
        if(NExist(tp)<1) then
            accept=0
            return
        end if
        !Randomly choose a type tp molecule
        call random_number(rnd)
        chosen=ceiling(rnd*NExist(tp))
        !Energy difference
        del=-EnergyChosenCurrent2D(tp,chosen)
        !Accept or reject this trial
        call random_number(rnd)
        accept=rnd<NExist(tp)/TranslationalPartitionFunction(tp)*exp(-(ChemicalPotential(tp)+del)/temperature)
    !Insertion
    else
        !Check the availability of workspace for the molecule to be inserted
        if(NExist(tp)==MaxPossibleN(tp)) then
            if(MaxPossibleNReachedWarning(tp)) then
                write(*,'(A33,I2,A17)')'Maximum possible number for type',tp,'molecule reached'
                write(*,'(A61)')'Please increase it before run grand canonical ensemble again'
                MaxPossibleNReachedWarning(tp)=0
            end if
            accept=0
            return
        end if
        chosen=NExist(tp)+1
        !Insert the new molecule randomly and uniformly
        call InsertMolecule2D(DisplacementVector,RotationQuaternion,.true.)
        mtp(tp).mol(chosen).Coordinate=struct(DisplacementVector,RotationQuaternion,&
            MoleculeInput(tp).RefConfig,MoleculeInput(tp).NAtoms)
        !z direction is confined by the slit wall, check if the atom is outside
        do i=1,MoleculeInput(tp).NAtoms
            if(abs(mtp(tp).mol(chosen).Coordinate(3,i))>HalfLength(3)) then
                accept=0
                return
            end if
        end do
        !Energy difference
        del=EnergyChosenUpdated2D(tp,chosen)
        !Accept or reject this trial
        call random_number(rnd)
        accept=rnd<TranslationalPartitionFunction(tp)/(NExist(tp)+1)*exp((ChemicalPotential(tp)-del)/temperature)
    end if
    !Nothing is changed, so nothing needs to be done if rejected
    if(accept) then
        if(direction<0.5d0) then
            NExist(tp)=NExist(tp)-1
            NTotal=NTotal-1
            !Delete the chosen molecule by flushing it with the last exist one, old=new by the way
            !Single-body
            mtp(tp).mol(chosen)=mtpold(tp).mol(NExistOld(tp))
            mtpold(tp).mol(chosen)=mtp(tp).mol(chosen)
            !Many-body
            do j=1,tp-1
                forall(i=1:NExistOld(j))
                    intertp(tp,j).molab(chosen,i)=intertpold(tp,j).molab(NExistOld(tp),i)
                    intertpold(tp,j).molab(chosen,i)=intertp(tp,j).molab(chosen,i)
                end forall
            end do
            forall(i=1:chosen-1)
                intertp(tp,tp).molab(chosen,i)=intertpold(tp,tp).molab(NExistOld(tp),i)
                intertpold(tp,tp).molab(chosen,i)=intertp(tp,tp).molab(chosen,i)
            end forall
            !molab(NExistOld(tp),chosen) is not needed as mol(NExistOld(tp)) is deleted
            forall(i=chosen+1:NExist(tp))
                intertp(tp,tp).molab(i,chosen)=intertpold(tp,tp).molab(NExistOld(tp),i)
                intertpold(tp,tp).molab(i,chosen)=intertp(tp,tp).molab(i,chosen)
            end forall
            do j=tp+1,MoleculeKinds
                forall(i=1:NExistOld(j))
                    intertp(j,tp).molab(i,chosen)=intertpold(j,tp).molab(i,NExistOld(tp))
                    intertpold(j,tp).molab(i,chosen)=intertp(j,tp).molab(i,chosen)
                end forall
            end do
            NExistOld(tp)=NExist(tp)
        else
            NExist(tp)=NExist(tp)+1
            NTotal=NTotal+1
            !old=new
            NExistOld(tp)=NExist(tp)
            !Single-body
            mtpold(tp).mol(chosen)=mtp(tp).mol(chosen)
            !Many-body
            do j=1,tp-1
                forall(i=1:NExist(j))
                    intertpold(tp,j).molab(chosen,i)=intertp(tp,j).molab(chosen,i)
                end forall
            end do
            !The chosen molecule is the last one so no need for the latter loop
            forall(i=1:chosen-1)
                intertpold(tp,tp).molab(chosen,i)=intertp(tp,tp).molab(chosen,i)
            end forall
            do j=tp+1,MoleculeKinds
                forall(i=1:NExist(j))
                    intertpold(j,tp).molab(i,chosen)=intertp(j,tp).molab(i,chosen)
                end forall
            end do
        end if
    end if
end subroutine Insertion2D
!For 3D
subroutine Insertion3D(tp,chosen,direction,accept)
    integer,intent(out)::tp,chosen
    real*8,intent(out)::direction
    logical,intent(out)::accept
    integer::i,j
    real*8::rnd,del
    !Randomly choose a type of molecule
    call random_number(rnd)
    tp=MoleculeKinds
    do i=1,MoleculeKinds-1
        if(rnd<ExChoose(i)) then
            tp=i
            exit
        end if
    end do
    !Determine the trial to make
    call random_number(direction)
    !Deletion
    if(direction<0.5d0) then
        !Check the availability of the molecule to be deleted
        if(NExist(tp)<1) then
            accept=0
            return
        end if
        !Randomly choose a type tp molecule
        call random_number(rnd)
        chosen=ceiling(rnd*NExist(tp))
        !Energy difference
        del=-EnergyChosenCurrent3D(tp,chosen)
        !Accept or reject this trial
        call random_number(rnd)
        accept=rnd<NExist(tp)/TranslationalPartitionFunction(tp)*exp(-(ChemicalPotential(tp)+del)/temperature)
    !Insertion
    else
        !Check the availability of workspace for the molecule to be inserted
        if(NExist(tp)==MaxPossibleN(tp)) then
            if(MaxPossibleNReachedWarning(tp)) then
                write(*,'(A33,I2,A17)')'Maximum possible number for type',tp,'molecule reached'
                write(*,'(A61)')'Please increase it before run grand canonical ensemble again'
                MaxPossibleNReachedWarning(tp)=0
            end if
            accept=0
            return
        end if
        chosen=NExist(tp)+1
        !Insert the new molecule randomly and uniformly
        call InsertMolecule3D(DisplacementVector,RotationQuaternion)
        mtp(tp).mol(chosen).Coordinate=struct(DisplacementVector,RotationQuaternion,&
            MoleculeInput(tp).RefConfig,MoleculeInput(tp).NAtoms)
        !Energy difference
        del=EnergyChosenUpdated3D(tp,chosen)
        !Accept or reject this trial
        call random_number(rnd)
        accept=rnd<TranslationalPartitionFunction(tp)/(NExist(tp)+1)*exp((ChemicalPotential(tp)-del)/temperature)
    end if
    !Nothing is changed, so nothing needs to be done if rejected
    if(accept) then
        if(direction<0.5d0) then
            NExist(tp)=NExist(tp)-1
            NTotal=NTotal-1
            !Delete the chosen molecule by flushing it with the last exist one, old=new by the way
            !Single-body
            mtp(tp).mol(chosen)=mtpold(tp).mol(NExistOld(tp))
            mtpold(tp).mol(chosen)=mtp(tp).mol(chosen)
            !Many-body
            do j=1,tp-1
                forall(i=1:NExistOld(j))
                    intertp(tp,j).molab(chosen,i)=intertpold(tp,j).molab(NExistOld(tp),i)
                    intertpold(tp,j).molab(chosen,i)=intertp(tp,j).molab(chosen,i)
                end forall
            end do
            forall(i=1:chosen-1)
                intertp(tp,tp).molab(chosen,i)=intertpold(tp,tp).molab(NExistOld(tp),i)
                intertpold(tp,tp).molab(chosen,i)=intertp(tp,tp).molab(chosen,i)
            end forall
            !molab(NExistOld(tp),chosen) is not needed as mol(NExistOld(tp)) is deleted
            forall(i=chosen+1:NExist(tp))
                intertp(tp,tp).molab(i,chosen)=intertpold(tp,tp).molab(NExistOld(tp),i)
                intertpold(tp,tp).molab(i,chosen)=intertp(tp,tp).molab(i,chosen)
            end forall
            do j=tp+1,MoleculeKinds
                forall(i=1:NExistOld(j))
                    intertp(j,tp).molab(i,chosen)=intertpold(j,tp).molab(i,NExistOld(tp))
                    intertpold(j,tp).molab(i,chosen)=intertp(j,tp).molab(i,chosen)
                end forall
            end do
            NExistOld(tp)=NExist(tp)
        else
            NExist(tp)=NExist(tp)+1
            NTotal=NTotal+1
            !old=new
            NExistOld(tp)=NExist(tp)
            !Single-body
            mtpold(tp).mol(chosen)=mtp(tp).mol(chosen)
            !Many-body
            do j=1,tp-1
                forall(i=1:NExist(j))
                    intertpold(tp,j).molab(chosen,i)=intertp(tp,j).molab(chosen,i)
                end forall
            end do
            !The chosen molecule is the last one so no need for the latter loop
            forall(i=1:chosen-1)
                intertpold(tp,tp).molab(chosen,i)=intertp(tp,tp).molab(chosen,i)
            end forall
            do j=tp+1,MoleculeKinds
                forall(i=1:NExist(j))
                    intertpold(j,tp).molab(i,chosen)=intertp(j,tp).molab(i,chosen)
                end forall
            end do
        end if
    end if
end subroutine Insertion3D

!Perform a single molecule forward or backward reaction A=B, for R type ensemble only
!Type rtp A=B reaction, chosenr-th type r molecule turns into chosenp-th type p molecule
!If direction<0.5, forward reaction; else backward reaction
!For 1D
subroutine ReactAB1D(rtp,r,p,chosenr,chosenp,direction,accept)
    integer,intent(in)::rtp
    integer,intent(out)::r,p,chosenr,chosenp
    real*8,intent(out)::direction
    logical,intent(out)::accept
    integer::i,j
    real*8::rnd,del,temp
    !Determine the trial to make
    call random_number(direction)
    !Forward
    if(direction<0.5d0) then
        !Check the availability of the reactant
        if(NExist(ReactionInput(rtp).ReactantSequence(1))<1) then
            accept=0
            return
        end if
        r=ReactionInput(rtp).ReactantSequence(1)
        p=ReactionInput(rtp).ProductSequence(1)
    !Backward
    else
        !Check the availability of the product
        if(NExist(ReactionInput(rtp).ProductSequence(1))<1) then
            accept=0
            return
        end if
        r=ReactionInput(rtp).ProductSequence(1)
        p=ReactionInput(rtp).ReactantSequence(1)
    end if
    !Choose a reactant molecule
    call random_number(rnd)
    chosenr=ceiling(rnd*NExist(r))
    !Delete it by flushing the chosen one with the last exist one
    !Here only flush the coordinate and the existance number
    !Single-body and many-body interactions will be flushed once accepted
    mtp(r).mol(chosenr).Coordinate=mtp(r).mol(NExist(r)).Coordinate
    NExist(r)=NExist(r)-1
    !Compute its energy
    del=-EnergyChosenCurrent1D(r,chosenr)
    !Insert product at the cavity left by the reactant
        NExist(p)=NExist(p)+1
        chosenp=NExist(p)
        !Compute the reactant's displacement vector relative to its reference configuration
        DisplacementVector=mtpold(r).mol(chosenr).Coordinate(:,1)-MoleculeInput(r).RefConfig(:,1)
        !Form the unit quaternion of the product orientation
            do while(1)
                call random_number(QuaternionSupporter1)
                QuaternionSupporter1=QuaternionSupporter1*2d0-1d0
                if(dot_product(QuaternionSupporter1,QuaternionSupporter1)<1d0) then
                    exit
                end if
            end do
            do while(1)
                call random_number(QuaternionSupporter2)
                QuaternionSupporter2=QuaternionSupporter2*2d0-1d0
                if(dot_product(QuaternionSupporter2,QuaternionSupporter2)<1d0) then
                    exit
                end if
            end do
            RotationQuaternion(1:2)=QuaternionSupporter1
            RotationQuaternion(3:4)=QuaternionSupporter2&
                *Sqrt((1d0-dot_product(QuaternionSupporter1,QuaternionSupporter1))&
                /dot_product(QuaternionSupporter2,QuaternionSupporter2))
        mtp(p).mol(chosenp).Coordinate=struct(DisplacementVector,RotationQuaternion,&
            MoleculeInput(p).RefConfig,MoleculeInput(p).NAtoms)
        !xy plane is confined by the cylindrical wall, check if the atom is outside
        do i=1,MoleculeInput(p).NAtoms
            if(dot_product(mtp(p).mol(chosenp).Coordinate(2:3,i),mtp(p).mol(chosenp).Coordinate(2:3,i))>radiussq) then
                accept=0
                !new=old
                NExist(r)=NExistOld(r)
                NExist(p)=NExistOld(p)
                !Coordinate only, no need for single-body interaction
                mtp(r).mol(chosenr).Coordinate=mtpold(r).mol(chosenr).Coordinate
                return
            end if
        end do
    !Compute the intermolecule energy of the product
    del=del+EnergyChosenUpdated1D(p,chosenp)
    !Add the ideal gas Helmholtz free energy change
    if(direction<0.5d0) then
        del=del+ReactionInput(rtp).deltaAig
    else
        del=del-ReactionInput(rtp).deltaAig
    end if
    !Accept or reject this trial
    del=del-temperature*Log(Real(NExistOld(r)/NExist(p)))
    accept=del<ImpossibleDel
    if(accept.and.del>0d0) then
        call random_number(rnd)
        accept=rnd<exp(-del/temperature)
    end if
    if(accept) then
        !Complete flushing the chosen reactant molecule with the last exist reactant molecule
        !Flushing single-body interaction, old(r)=new(r) by the way
        mtp(r).mol(chosenr).WallLJ=mtpold(r).mol(NExistOld(r)).WallLJ
        mtpold(r).mol(chosenr)=mtp(r).mol(chosenr)
        !Flushing many-body interaction, old(r)=new(r) by the way
        do j=1,r-1
            forall(i=1:NExistOld(j))
                intertp(r,j).molab(chosenr,i)=intertpold(r,j).molab(NExistOld(r),i)
                intertpold(r,j).molab(chosenr,i)=intertp(r,j).molab(chosenr,i)
            end forall
        end do
        forall(i=1:chosenr-1)
            intertp(r,r).molab(chosenr,i)=intertpold(r,r).molab(NExistOld(r),i)
            intertpold(r,r).molab(chosenr,i)=intertp(r,r).molab(chosenr,i)
        end forall
        !molab(NExistOld(r),chosenr) is not needed as mol(NExistOld(r)) is deleted
        forall(i=chosenr+1:NExist(r))
            intertp(r,r).molab(i,chosenr)=intertpold(r,r).molab(NExistOld(r),i)
            intertpold(r,r).molab(i,chosenr)=intertp(r,r).molab(i,chosenr)
        end forall
        do j=r+1,MoleculeKinds
            forall(i=1:NExistOld(j))
                intertp(j,r).molab(i,chosenr)=intertpold(j,r).molab(i,NExistOld(r))
                intertpold(j,r).molab(i,chosenr)=intertp(j,r).molab(i,chosenr)
            end forall
        end do
        !Complete old(r)=new(r)
        NExistOld(r)=NExist(r)
        !old(p)=new(p)
        NExistOld(p)=NExist(p)
        !Single-body
        mtpold(p).mol(chosenp)=mtp(p).mol(chosenp)
        !Many-body
        do j=1,p-1
            forall(i=1:NExist(j))
                intertpold(p,j).molab(chosenp,i)=intertp(p,j).molab(chosenp,i)
            end forall
        end do
        !chosenp=NExist(p) so no need for the latter loop
        forall(i=1:chosenp-1)
            intertpold(p,p).molab(chosenp,i)=intertp(p,p).molab(chosenp,i)
        end forall
        do j=p+1,MoleculeKinds
            forall(i=1:NExist(j))
                intertpold(j,p).molab(i,chosenp)=intertp(j,p).molab(i,chosenp)
            end forall
        end do
    else
        !new=old
        NExist(r)=NExistOld(r)
        NExist(p)=NExistOld(p)
        !Coordinate only, no need for single-body interaction
        mtp(r).mol(chosenr).Coordinate=mtpold(r).mol(chosenr).Coordinate
        !No need for many-body interaction
    end if
end subroutine ReactAB1D
!For 2D
subroutine ReactAB2D(rtp,r,p,chosenr,chosenp,direction,accept)
    integer,intent(in)::rtp
    integer,intent(out)::r,p,chosenr,chosenp
    real*8,intent(out)::direction
    logical,intent(out)::accept
    integer::i,j
    real*8::rnd,del,temp
    !Determine the trial to make
    call random_number(direction)
    !Forward
    if(direction<0.5d0) then
        !Check the availability of the reactant
        if(NExist(ReactionInput(rtp).ReactantSequence(1))<1) then
            accept=0
            return
        end if
        r=ReactionInput(rtp).ReactantSequence(1)
        p=ReactionInput(rtp).ProductSequence(1)
    !Backward
    else
        !Check the availability of the product
        if(NExist(ReactionInput(rtp).ProductSequence(1))<1) then
            accept=0
            return
        end if
        r=ReactionInput(rtp).ProductSequence(1)
        p=ReactionInput(rtp).ReactantSequence(1)
    end if
    !Choose a reactant molecule
    call random_number(rnd)
    chosenr=ceiling(rnd*NExist(r))
    !Delete it by flushing the chosen one with the last exist one
    !Here only flush the coordinate and the existance number
    !Single-body and many-body interactions will be flushed once accepted
    mtp(r).mol(chosenr).Coordinate=mtp(r).mol(NExist(r)).Coordinate
    NExist(r)=NExist(r)-1
    !Compute its energy
    del=-EnergyChosenCurrent2D(r,chosenr)
    !Insert product at the cavity left by the reactant
        NExist(p)=NExist(p)+1
        chosenp=NExist(p)
        !Compute the reactant's displacement vector relative to its reference configuration
        DisplacementVector=mtpold(r).mol(chosenr).Coordinate(:,1)-MoleculeInput(r).RefConfig(:,1)
        !Form the unit quaternion of the product orientation
            do while(1)
                call random_number(QuaternionSupporter1)
                QuaternionSupporter1=QuaternionSupporter1*2d0-1d0
                if(dot_product(QuaternionSupporter1,QuaternionSupporter1)<1d0) then
                    exit
                end if
            end do
            do while(1)
                call random_number(QuaternionSupporter2)
                QuaternionSupporter2=QuaternionSupporter2*2d0-1d0
                if(dot_product(QuaternionSupporter2,QuaternionSupporter2)<1d0) then
                    exit
                end if
            end do
            RotationQuaternion(1:2)=QuaternionSupporter1
            RotationQuaternion(3:4)=QuaternionSupporter2&
                *Sqrt((1d0-dot_product(QuaternionSupporter1,QuaternionSupporter1))&
                /dot_product(QuaternionSupporter2,QuaternionSupporter2))
        mtp(p).mol(chosenp).Coordinate=struct(DisplacementVector,RotationQuaternion,&
            MoleculeInput(p).RefConfig,MoleculeInput(p).NAtoms)
        !z direction is confined by the slit wall, check if the atom is outside
        do i=1,MoleculeInput(p).NAtoms
            if(abs(mtp(p).mol(chosenp).Coordinate(3,i))>HalfLength(3)) then
                accept=0
                !new=old
                NExist(r)=NExistOld(r)
                NExist(p)=NExistOld(p)
                !Coordinate only, no need for single-body interaction
                mtp(r).mol(chosenr).Coordinate=mtpold(r).mol(chosenr).Coordinate
                return
            end if
        end do
    !Compute the intermolecule energy of the product
    del=del+EnergyChosenUpdated2D(p,chosenp)
    !Add the ideal gas Helmholtz free energy change
    if(direction<0.5d0) then
        del=del+ReactionInput(rtp).deltaAig
    else
        del=del-ReactionInput(rtp).deltaAig
    end if
    !Accept or reject this trial
    del=del-temperature*Log(Real(NExistOld(r)/NExist(p)))
    accept=del<ImpossibleDel
    if(accept.and.del>0d0) then
        call random_number(rnd)
        accept=rnd<exp(-del/temperature)
    end if
    if(accept) then
        !Complete flushing the chosen reactant molecule with the last exist reactant molecule
        !Flushing single-body interaction, old(r)=new(r) by the way
        mtp(r).mol(chosenr).Total=mtpold(r).mol(NExistOld(r)).Total
        mtp(r).mol(chosenr).ESSBDependent2D=mtpold(r).mol(NExistOld(r)).ESSBDependent2D
        mtp(r).mol(chosenr).WallLJ=mtpold(r).mol(NExistOld(r)).WallLJ
        mtpold(r).mol(chosenr)=mtp(r).mol(chosenr)
        !Flushing many-body interaction, old(r)=new(r) by the way
        do j=1,r-1
            forall(i=1:NExistOld(j))
                intertp(r,j).molab(chosenr,i)=intertpold(r,j).molab(NExistOld(r),i)
                intertpold(r,j).molab(chosenr,i)=intertp(r,j).molab(chosenr,i)
            end forall
        end do
        forall(i=1:chosenr-1)
            intertp(r,r).molab(chosenr,i)=intertpold(r,r).molab(NExistOld(r),i)
            intertpold(r,r).molab(chosenr,i)=intertp(r,r).molab(chosenr,i)
        end forall
        !molab(NExistOld(r),chosenr) is not needed as mol(NExistOld(r)) is deleted
        forall(i=chosenr+1:NExist(r))
            intertp(r,r).molab(i,chosenr)=intertpold(r,r).molab(NExistOld(r),i)
            intertpold(r,r).molab(i,chosenr)=intertp(r,r).molab(i,chosenr)
        end forall
        do j=r+1,MoleculeKinds
            forall(i=1:NExistOld(j))
                intertp(j,r).molab(i,chosenr)=intertpold(j,r).molab(i,NExistOld(r))
                intertpold(j,r).molab(i,chosenr)=intertp(j,r).molab(i,chosenr)
            end forall
        end do
        !Complete old(r)=new(r)
        NExistOld(r)=NExist(r)
        !old(p)=new(p)
        NExistOld(p)=NExist(p)
        !Single-body
        mtpold(p).mol(chosenp)=mtp(p).mol(chosenp)
        !Many-body
        do j=1,p-1
            forall(i=1:NExist(j))
                intertpold(p,j).molab(chosenp,i)=intertp(p,j).molab(chosenp,i)
            end forall
        end do
        !chosenp=NExist(p) so no need for the latter loop
        forall(i=1:chosenp-1)
            intertpold(p,p).molab(chosenp,i)=intertp(p,p).molab(chosenp,i)
        end forall
        do j=p+1,MoleculeKinds
            forall(i=1:NExist(j))
                intertpold(j,p).molab(i,chosenp)=intertp(j,p).molab(i,chosenp)
            end forall
        end do
    else
        !new=old
        NExist(r)=NExistOld(r)
        NExist(p)=NExistOld(p)
        !Coordinate only, no need for single-body interaction
        mtp(r).mol(chosenr).Coordinate=mtpold(r).mol(chosenr).Coordinate
        !No need for many-body interaction
    end if
end subroutine ReactAB2D
!For 3D
subroutine ReactAB3D(rtp,r,p,chosenr,chosenp,direction,accept)
    integer,intent(in)::rtp
    integer,intent(out)::r,p,chosenr,chosenp
    real*8,intent(out)::direction
    logical,intent(out)::accept
    integer::i,j
    real*8::rnd,del,temp
    !Determine the trial to make
    call random_number(direction)
    !Forward
    if(direction<0.5d0) then
        !Check the availability of the reactant
        if(NExist(ReactionInput(rtp).ReactantSequence(1))<1) then
            accept=0
            return
        end if
        r=ReactionInput(rtp).ReactantSequence(1)
        p=ReactionInput(rtp).ProductSequence(1)
    !Backward
    else
        !Check the availability of the product
        if(NExist(ReactionInput(rtp).ProductSequence(1))<1) then
            accept=0
            return
        end if
        r=ReactionInput(rtp).ProductSequence(1)
        p=ReactionInput(rtp).ReactantSequence(1)
    end if
    !Choose a reactant molecule
    call random_number(rnd)
    chosenr=ceiling(rnd*NExist(r))
    !Delete it by flushing the chosen one with the last exist one
    !Here only flush the coordinate and the existance number
    !Single-body and many-body interactions will be flushed once accepted
    mtp(r).mol(chosenr).Coordinate=mtp(r).mol(NExist(r)).Coordinate
    NExist(r)=NExist(r)-1
    !Compute its energy
    del=-EnergyChosenCurrent3D(r,chosenr)
    !Insert product at the cavity left by the reactant
        NExist(p)=NExist(p)+1
        chosenp=NExist(p)
        !Compute the reactant's displacement vector relative to its reference configuration
        DisplacementVector=mtpold(r).mol(chosenr).Coordinate(:,1)-MoleculeInput(r).RefConfig(:,1)
        !Form the unit quaternion of the product orientation
            do while(1)
                call random_number(QuaternionSupporter1)
                QuaternionSupporter1=QuaternionSupporter1*2d0-1d0
                if(dot_product(QuaternionSupporter1,QuaternionSupporter1)<1d0) then
                    exit
                end if
            end do
            do while(1)
                call random_number(QuaternionSupporter2)
                QuaternionSupporter2=QuaternionSupporter2*2d0-1d0
                if(dot_product(QuaternionSupporter2,QuaternionSupporter2)<1d0) then
                    exit
                end if
            end do
            RotationQuaternion(1:2)=QuaternionSupporter1
            RotationQuaternion(3:4)=QuaternionSupporter2&
                *Sqrt((1d0-dot_product(QuaternionSupporter1,QuaternionSupporter1))&
                /dot_product(QuaternionSupporter2,QuaternionSupporter2))
        mtp(p).mol(chosenp).Coordinate=struct(DisplacementVector,RotationQuaternion,&
            MoleculeInput(p).RefConfig,MoleculeInput(p).NAtoms)
    !Compute the intermolecule energy of the product
    del=del+EnergyChosenUpdated3D(p,chosenp)
    !Add the ideal gas Helmholtz free energy change
    if(direction<0.5d0) then
        del=del+ReactionInput(rtp).deltaAig
    else
        del=del-ReactionInput(rtp).deltaAig
    end if
    !Accept or reject this trial
    del=del-temperature*Log(Real(NExistOld(r)/NExist(p)))
    accept=del<ImpossibleDel
    if(accept.and.del>0d0) then
        call random_number(rnd)
        accept=rnd<exp(-del/temperature)
    end if
    if(accept) then
        !Complete flushing the chosen reactant molecule with the last exist reactant molecule
        !Currently no position dependent single-body interaction for 3D, old(r)=new(r) only
        mtpold(r).mol(chosenr).Coordinate=mtp(r).mol(chosenr).Coordinate
        !Flushing many-body interaction, old(r)=new(r) by the way
        do j=1,r-1
            forall(i=1:NExistOld(j))
                intertp(r,j).molab(chosenr,i)=intertpold(r,j).molab(NExistOld(r),i)
                intertpold(r,j).molab(chosenr,i)=intertp(r,j).molab(chosenr,i)
            end forall
        end do
        forall(i=1:chosenr-1)
            intertp(r,r).molab(chosenr,i)=intertpold(r,r).molab(NExistOld(r),i)
            intertpold(r,r).molab(chosenr,i)=intertp(r,r).molab(chosenr,i)
        end forall
        !molab(NExistOld(r),chosenr) is not needed as mol(NExistOld(r)) is deleted
        forall(i=chosenr+1:NExist(r))
            intertp(r,r).molab(i,chosenr)=intertpold(r,r).molab(NExistOld(r),i)
            intertpold(r,r).molab(i,chosenr)=intertp(r,r).molab(i,chosenr)
        end forall
        do j=r+1,MoleculeKinds
            forall(i=1:NExistOld(j))
                intertp(j,r).molab(i,chosenr)=intertpold(j,r).molab(i,NExistOld(r))
                intertpold(j,r).molab(i,chosenr)=intertp(j,r).molab(i,chosenr)
            end forall
        end do
        !Complete old(r)=new(r)
        NExistOld(r)=NExist(r)
        !old(p)=new(p)
        NExistOld(p)=NExist(p)
        !Currently no position dependent single-body interaction for 3D, coordinate only
        mtpold(p).mol(chosenp).Coordinate=mtp(p).mol(chosenp).Coordinate
        !Many-body
        do j=1,p-1
            forall(i=1:NExist(j))
                intertpold(p,j).molab(chosenp,i)=intertp(p,j).molab(chosenp,i)
            end forall
        end do
        !chosenp=NExist(p) so no need for the latter loop
        forall(i=1:chosenp-1)
            intertpold(p,p).molab(chosenp,i)=intertp(p,p).molab(chosenp,i)
        end forall
        do j=p+1,MoleculeKinds
            forall(i=1:NExist(j))
                intertpold(j,p).molab(i,chosenp)=intertp(j,p).molab(i,chosenp)
            end forall
        end do
    else
        !new=old
        NExist(r)=NExistOld(r)
        NExist(p)=NExistOld(p)
        !Coordinate only, no need for single-body interaction
        mtp(r).mol(chosenr).Coordinate=mtpold(r).mol(chosenr).Coordinate
        !No need for many-body interaction
    end if
end subroutine ReactAB3D

!Perform a single molecule forward or backward reaction A+B=C, for R type ensemble only
!B has to be a small molecule for detailed balance issue
!Type rtp A+B=C reaction
!For forward, chosenr(1)-th type r(1) and chosenr(2)-th type r(2) molecule turn into chosenp-th type p molecule
!For backward, chosenp-th type p molecule turns into chosenr(1)-th type r(1) and chosenr(2)-th type r(2) molecule
!If direction<0.5, forward reaction; else backward reaction
!For 1D
subroutine ReactABC1D(rtp,r,p,chosenr,chosenp,direction,accept)
    integer,intent(in)::rtp
    integer,dimension(2),intent(out)::r,chosenr
    integer,intent(out)::p,chosenp
    real*8,intent(out)::direction
    logical,intent(out)::accept
    integer::i,j,k,ii,jj
    real*8::rnd,del
    real*8,dimension(2)::choose
    !To make notation simpler and reduce addressing
    r=ReactionInput(rtp).ReactantSequence
    p=ReactionInput(rtp).ProductSequence(1)
    !Determine the trial to make
    call random_number(direction)
    !Forward
    if(direction<0.5d0) then
        !Check the availability of the reactant
        if(NExist(r(1))<1.or.NExist(r(2))<1) then
            accept=0
            return
        end if
        !Choose two reactants and delete them by flushing the chosen ones with the last exist ones
        !Here only flush the coordinate and the existance number
        !Single-body and many-body interactions will be flushed once accepted
        call random_number(choose)
        forall(i=1:2)
            chosenr(i)=ceiling(choose(i)*NExist(r(i)))
            mtp(r(i)).mol(chosenr(i)).Coordinate=mtpold(r(i)).mol(NExist(r(i))).Coordinate
            NExist(r(i))=NExist(r(i))-1
        end forall
        !Compute their energy
        del=-EnergyChosenCurrent1D(r(1),chosenr(1))-EnergyChosenCurrent1D(r(2),chosenr(2))&
            !By input rule, the smaller molecule has the smaller type number
            +intertpold(r(1),r(2)).molab(chosenr(1),chosenr(2)).Total
        !Insert product at the cavity left by reactant A
        NExist(p)=NExist(p)+1
        chosenp=NExist(p)
        !The position of the centre of mass of the larger reactant
        DisplacementVector=mtpold(r(1)).mol(chosenr(1)).Coordinate(:,1)-MoleculeInput(r(1)).RefConfig(:,1)
        !Form the unit quaternion of the product orientation
            do while(1)
                call random_number(QuaternionSupporter1)
                QuaternionSupporter1=QuaternionSupporter1*2d0-1d0
                if(dot_product(QuaternionSupporter1,QuaternionSupporter1)<1d0) then
                    exit
                end if
            end do
            do while(1)
                call random_number(QuaternionSupporter2)
                QuaternionSupporter2=QuaternionSupporter2*2d0-1d0
                if(dot_product(QuaternionSupporter2,QuaternionSupporter2)<1d0) then
                    exit
                end if
            end do
            RotationQuaternion(1:2)=QuaternionSupporter1
            RotationQuaternion(3:4)=QuaternionSupporter2&
                *Sqrt((1d0-dot_product(QuaternionSupporter1,QuaternionSupporter1))&
                /dot_product(QuaternionSupporter2,QuaternionSupporter2))
        mtp(p).mol(chosenp).Coordinate=struct(DisplacementVector,RotationQuaternion,&
            MoleculeInput(p).RefConfig,MoleculeInput(p).NAtoms)
        !xy plane is confined by the cylindrical wall, check if the atom is outside
        do i=1,MoleculeInput(p).NAtoms
            if(dot_product(mtp(p).mol(chosenp).Coordinate(2:3,i),mtp(p).mol(chosenp).Coordinate(2:3,i))>radiussq) then
                accept=0
                !new=old
                NExist(r(1))=NExistOld(r(1))
                NExist(r(2))=NExistOld(r(2))
                NExist(p)=NExistOld(p)
                !Coordinate only, no need for single-body interaction
                forall(j=1:2)
                    mtp(r(j)).mol(chosenr(j)).Coordinate=mtpold(r(j)).mol(chosenr(j)).Coordinate
                end forall
                return
            end if
        end do
        del=del+EnergyChosenUpdated1D(p,chosenp)+ReactionInput(rtp).deltaAig&
            -temperature*Log(Real(NExistOld(r(1))*NExistOld(r(2))/NExist(p)))
    !Backward
    else
        !Check the availability of the product
        if(NExist(p)<1) then
            accept=0
            return
        end if
        !Choose a product
        call random_number(rnd)
        chosenp=ceiling(rnd*NExist(p))
        !Delete it by flushing the chosen one with the last exist one
        !Here only flush the coordinate and the existance number
        !Single-body and many-body interactions will be flushed once accepted
        mtp(p).mol(chosenp).Coordinate=mtpold(p).mol(NExist(p)).Coordinate
        NExist(p)=NExist(p)-1
        !Compute its energy
        del=-EnergyChosenCurrent1D(p,chosenp)
        !Insert reactant A at the cavity left by product
        forall(i=1:2)
            chosenr(i)=NExist(r(i))+1
        end forall
        !The position of the centre of mass of the product
        DisplacementVector=mtpold(p).mol(chosenp).Coordinate(:,1)-MoleculeInput(p).RefConfig(:,1)
        !Form the unit quaternion of the reactant orientation
            do while(1)
                call random_number(QuaternionSupporter1)
                QuaternionSupporter1=QuaternionSupporter1*2d0-1d0
                if(dot_product(QuaternionSupporter1,QuaternionSupporter1)<1d0) then
                    exit
                end if
            end do
            do while(1)
                call random_number(QuaternionSupporter2)
                QuaternionSupporter2=QuaternionSupporter2*2d0-1d0
                if(dot_product(QuaternionSupporter2,QuaternionSupporter2)<1d0) then
                    exit
                end if
            end do
            RotationQuaternion(1:2)=QuaternionSupporter1
            RotationQuaternion(3:4)=QuaternionSupporter2&
                *Sqrt((1d0-dot_product(QuaternionSupporter1,QuaternionSupporter1))&
                /dot_product(QuaternionSupporter2,QuaternionSupporter2))
        mtp(r(1)).mol(chosenr(1)).Coordinate=struct(DisplacementVector,RotationQuaternion,&
            MoleculeInput(r(1)).RefConfig,MoleculeInput(r(1)).NAtoms)
        !xy plane is confined by the cylindrical wall, check if the atom is outside
        do i=1,MoleculeInput(r(1)).NAtoms
            if(dot_product(mtp(r(1)).mol(chosenr(1)).Coordinate(2:3,i),mtp(r(1)).mol(chosenr(1)).Coordinate(2:3,i))>radiussq) then
                accept=0
                !new=old
                NExist(p)=NExistOld(p)
                !Coordinate only, no need for single-body interaction
                mtp(p).mol(chosenp).Coordinate=mtpold(p).mol(chosenp).Coordinate
                return
            end if
        end do
        NExist(r(1))=NExist(r(1))+1
        del=del+EnergyChosenUpdated1D(r(1),chosenr(1))
        !Insert B randomly and uniformly
        call InsertMolecule1D(DisplacementVector,RotationQuaternion,.true.)
        mtp(r(2)).mol(chosenr(2)).Coordinate=struct(DisplacementVector,RotationQuaternion,&
            MoleculeInput(r(2)).RefConfig,MoleculeInput(r(2)).NAtoms)
        !xy plane is confined by the cylindrical wall, check if the atom is outside
        do i=1,MoleculeInput(r(2)).NAtoms
            if(dot_product(mtp(r(2)).mol(chosenr(2)).Coordinate(2:3,i),mtp(r(2)).mol(chosenr(2)).Coordinate(2:3,i))>radiussq) then
                accept=0
                !new=old
                NExist(r(1))=NExistOld(r(1))
                NExist(p)=NExistOld(p)
                !Coordinate only, no need for single-body interaction
                mtp(p).mol(chosenp).Coordinate=mtpold(p).mol(chosenp).Coordinate
                return
            end if
        end do
        NExist(r(2))=NExist(r(2))+1
        del=del+EnergyChosenUpdated1D(r(2),chosenr(2))-ReactionInput(rtp).deltaAig&
            -temperature*Log(Real(NExistOld(p)/NExist(r(1))/NExist(r(2))))
    end if
    !Accept or reject this trial
    accept=del<ImpossibleDel
    if(accept.and.del>0d0) then
        call random_number(rnd)
        accept=rnd<exp(-del/temperature)
    end if
    if(accept) then
        if(direction<0.5d0) then
            NTotal=NTotal-1
            !Complete flushing the chosen reactant molecules with the last exist reactant molecules
            do k=1,2
                !Flushing single-body interaction, old(r_k)=new(r_k) by the way
                mtp(r(k)).mol(chosenr(k)).WallLJ=mtpold(r(k)).mol(NExistOld(r(k))).WallLJ
                mtpold(r(k)).mol(chosenr(k)).Coordinate=mtp(r(k)).mol(chosenr(k)).Coordinate
                mtpold(r(k)).mol(chosenr(k)).WallLJ=mtp(r(k)).mol(chosenr(k)).WallLJ
                !Flushing many-body interaction, old(r_k)=new(r_k) by the way
                do j=1,r(k)-1
                    forall(i=1:NExistOld(j))
                        intertp(r(k),j).molab(chosenr(k),i)=intertpold(r(k),j).molab(NExistOld(r(k)),i)
                        intertpold(r(k),j).molab(chosenr(k),i)=intertp(r(k),j).molab(chosenr(k),i)
                    end forall
                end do
                forall(i=1:chosenr(k)-1)
                    intertp(r(k),r(k)).molab(chosenr(k),i)=intertpold(r(k),r(k)).molab(NExistOld(r(k)),i)
                    intertpold(r(k),r(k)).molab(chosenr(k),i)=intertp(r(k),r(k)).molab(chosenr(k),i)
                end forall
                !Same type molab(NExistOld(r),chosenr) is not needed as mol(NExistOld(r)) is deleted
                forall(i=chosenr(k)+1:NExist(r(k)))
                    intertp(r(k),r(k)).molab(i,chosenr(k))=intertpold(r(k),r(k)).molab(NExistOld(r(k)),i)
                    intertpold(r(k),r(k)).molab(i,chosenr(k))=intertp(r(k),r(k)).molab(i,chosenr(k))
                end forall
                do j=r(k)+1,MoleculeKinds
                    forall(i=1:NExistOld(j))
                        intertp(j,r(k)).molab(i,chosenr(k))=intertpold(j,r(k)).molab(i,NExistOld(r(k)))
                        intertpold(j,r(k)).molab(i,chosenr(k))=intertp(j,r(k)).molab(i,chosenr(k))
                    end forall
                end do
                !Complete old(r)=new(r)
                NExistOld(r(k))=NExist(r(k))
            end do
            !old(p)=new(p)
            NExistOld(p)=NExist(p)
            !Single-body
            mtpold(p).mol(chosenp).Coordinate=mtp(p).mol(chosenp).Coordinate
            mtpold(p).mol(chosenp).WallLJ=mtp(p).mol(chosenp).WallLJ
            !Many-body
            do j=1,p-1
                forall(i=1:NExist(j))
                    intertpold(p,j).molab(chosenp,i)=intertp(p,j).molab(chosenp,i)
                end forall
            end do
            !chosenp=NExist(p) so no need for the latter loop
            forall(i=1:chosenp-1)
                intertpold(p,p).molab(chosenp,i)=intertp(p,p).molab(chosenp,i)
            end forall
            do j=p+1,MoleculeKinds
                forall(i=1:NExist(j))
                    intertpold(j,p).molab(i,chosenp)=intertp(j,p).molab(i,chosenp)
                end forall
            end do
        else
            NTotal=NTotal+1
            !Complete flushing the chosen product molecule with the last exist product molecule
            !Flushing single-body interaction, old(p)=new(p) by the way
            mtp(p).mol(chosenp).WallLJ=mtpold(p).mol(NExistOld(p)).WallLJ
            mtpold(p).mol(chosenp).Coordinate=mtp(p).mol(chosenp).Coordinate
            mtpold(p).mol(chosenp).WallLJ=mtp(p).mol(chosenp).WallLJ
            !Flushing many-body interaction, old(p)=new(p) by the way
            do j=1,p-1
                forall(i=1:NExistOld(j))
                    intertp(p,j).molab(chosenp,i)=intertpold(p,j).molab(NExistOld(p),i)
                    intertpold(p,j).molab(chosenp,i)=intertp(p,j).molab(chosenp,i)
                end forall
            end do
            forall(i=1:chosenp-1)
                intertp(p,p).molab(chosenp,i)=intertpold(p,p).molab(NExistOld(p),i)
                intertpold(p,p).molab(chosenp,i)=intertp(p,p).molab(chosenp,i)
            end forall
            forall(i=chosenp+1:NExistOld(p))
                intertp(p,p).molab(i,chosenp)=intertpold(p,p).molab(NExistOld(p),i)
                intertpold(p,p).molab(i,chosenp)=intertp(p,p).molab(i,chosenp)
            end forall
            do j=p+1,MoleculeKinds
                forall(i=1:NExistOld(j))
                    intertp(j,p).molab(i,chosenp)=intertpold(j,p).molab(i,NExistOld(p))
                    intertpold(j,p).molab(i,chosenp)=intertp(j,p).molab(i,chosenp)
                end forall
            end do
            !Complete old(p)=new(p)
            NExistOld(p)=NExist(p)
            !old(r)=new(r)
            do k=1,2
                NExistOld(r(k))=NExist(r(k))
                !Single-body
                mtpold(r(k)).mol(chosenr(k)).Coordinate=mtp(r(k)).mol(chosenr(k)).Coordinate
                mtpold(r(k)).mol(chosenr(k)).WallLJ=mtp(r(k)).mol(chosenr(k)).WallLJ
                !Many-body
                do j=1,r(k)-1
                    forall(i=1:NExist(j))
                        intertpold(r(k),j).molab(chosenr(k),i)=intertp(r(k),j).molab(chosenr(k),i)
                    end forall
                end do
                forall(i=1:chosenr(k)-1)
                    intertpold(r(k),r(k)).molab(chosenr(k),i)=intertp(r(k),r(k)).molab(chosenr(k),i)
                end forall
                forall(i=chosenr(k)+1:NExist(r(k)))
                    intertpold(r(k),r(k)).molab(i,chosenr(k))=intertp(r(k),r(k)).molab(i,chosenr(k))
                end forall
                do j=r(k)+1,MoleculeKinds
                    forall(i=1:NExist(j))
                        intertpold(j,r(k)).molab(i,chosenr(k))=intertp(j,r(k)).molab(i,chosenr(k))
                    end forall
                end do
            end do
        end if
    else
        NExist(r(1))=NExistOld(r(1))
        NExist(r(2))=NExistOld(r(2))
        NExist(p)=NExistOld(p)
        !Coordinate only, no need for single-body interaction
        if(direction<0.5d0) then
            forall(i=1:2)
                mtp(r(i)).mol(chosenr(i)).Coordinate=mtpold(r(i)).mol(chosenr(i)).Coordinate
            end forall
        else
            mtp(p).mol(chosenp).Coordinate=mtpold(p).mol(chosenp).Coordinate
        end if
        !No need for many-body interaction
    end if
end subroutine ReactABC1D
!For 2D
subroutine ReactABC2D(rtp,r,p,chosenr,chosenp,direction,accept)
    integer,intent(in)::rtp
    integer,dimension(2),intent(out)::r,chosenr
    integer,intent(out)::p,chosenp
    real*8,intent(out)::direction
    logical,intent(out)::accept
    integer::i,j,k,ii,jj
    real*8::rnd,del
    real*8,dimension(2)::choose
    !To make notation simpler and reduce addressing
    r=ReactionInput(rtp).ReactantSequence
    p=ReactionInput(rtp).ProductSequence(1)
    !Determine the trial to make
    call random_number(direction)
    !Forward
    if(direction<0.5d0) then
        !Check the availability of the reactant
        if(NExist(r(1))<1.or.NExist(r(2))<1) then
            accept=0
            return
        end if
        !Choose two reactants and delete them by flushing the chosen ones with the last exist ones
        !Here only flush the coordinate and the existance number
        !Single-body and many-body interactions will be flushed once accepted
        call random_number(choose)
        forall(i=1:2)
            chosenr(i)=ceiling(choose(i)*NExist(r(i)))
            mtp(r(i)).mol(chosenr(i)).Coordinate=mtpold(r(i)).mol(NExist(r(i))).Coordinate
            NExist(r(i))=NExist(r(i))-1
        end forall
        !Compute their energy
        del=-EnergyChosenCurrent2D(r(1),chosenr(1))-EnergyChosenCurrent2D(r(2),chosenr(2))&
            !By input rule, the smaller molecule has the smaller type number
            +intertpold(r(1),r(2)).molab(chosenr(1),chosenr(2)).Total
        !Insert product at the cavity left by reactant A
        NExist(p)=NExist(p)+1
        chosenp=NExist(p)
        !The position of the centre of mass of the larger reactant
        DisplacementVector=mtpold(r(1)).mol(chosenr(1)).Coordinate(:,1)-MoleculeInput(r(1)).RefConfig(:,1)
        !Form the unit quaternion of the product orientation
            do while(1)
                call random_number(QuaternionSupporter1)
                QuaternionSupporter1=QuaternionSupporter1*2d0-1d0
                if(dot_product(QuaternionSupporter1,QuaternionSupporter1)<1d0) then
                    exit
                end if
            end do
            do while(1)
                call random_number(QuaternionSupporter2)
                QuaternionSupporter2=QuaternionSupporter2*2d0-1d0
                if(dot_product(QuaternionSupporter2,QuaternionSupporter2)<1d0) then
                    exit
                end if
            end do
            RotationQuaternion(1:2)=QuaternionSupporter1
            RotationQuaternion(3:4)=QuaternionSupporter2&
                *Sqrt((1d0-dot_product(QuaternionSupporter1,QuaternionSupporter1))&
                /dot_product(QuaternionSupporter2,QuaternionSupporter2))
        mtp(p).mol(chosenp).Coordinate=struct(DisplacementVector,RotationQuaternion,&
            MoleculeInput(p).RefConfig,MoleculeInput(p).NAtoms)
        !z direction is confined by the slit wall, check if the atom is outside
            do i=1,MoleculeInput(p).NAtoms
                if(abs(mtp(p).mol(chosenp).Coordinate(3,i))>HalfLength(3)) then
                    accept=0
                    !new=old
                    NExist(r(1))=NExistOld(r(1))
                    NExist(r(2))=NExistOld(r(2))
                    NExist(p)=NExistOld(p)
                    !Coordinate only, no need for single-body interaction
                    forall(j=1:2)
                        mtp(r(j)).mol(chosenr(j)).Coordinate=mtpold(r(j)).mol(chosenr(j)).Coordinate
                    end forall
                    return
                end if
            end do
        del=del+EnergyChosenUpdated2D(p,chosenp)+ReactionInput(rtp).deltaAig&
            -temperature*Log(Real(NExistOld(r(1))*NExistOld(r(2))/NExist(p)))
    !Backward
    else
        !Check the availability of the product
        if(NExist(p)<1) then
            accept=0
            return
        end if
        !Choose a product
        call random_number(rnd)
        chosenp=ceiling(rnd*NExist(p))
        !Delete it by flushing the chosen one with the last exist one
        !Here only flush the coordinate and the existance number
        !Single-body and many-body interactions will be flushed once accepted
        mtp(p).mol(chosenp).Coordinate=mtpold(p).mol(NExist(p)).Coordinate
        NExist(p)=NExist(p)-1
        !Compute its energy
        del=-EnergyChosenCurrent2D(p,chosenp)
        !Insert reactant A at the cavity left by product
        forall(i=1:2)
            chosenr(i)=NExist(r(i))+1
        end forall
        !The position of the centre of mass of the product
        DisplacementVector=mtpold(p).mol(chosenp).Coordinate(:,1)-MoleculeInput(p).RefConfig(:,1)
        !Form the unit quaternion of the reactant orientation
            do while(1)
                call random_number(QuaternionSupporter1)
                QuaternionSupporter1=QuaternionSupporter1*2d0-1d0
                if(dot_product(QuaternionSupporter1,QuaternionSupporter1)<1d0) then
                    exit
                end if
            end do
            do while(1)
                call random_number(QuaternionSupporter2)
                QuaternionSupporter2=QuaternionSupporter2*2d0-1d0
                if(dot_product(QuaternionSupporter2,QuaternionSupporter2)<1d0) then
                    exit
                end if
            end do
            RotationQuaternion(1:2)=QuaternionSupporter1
            RotationQuaternion(3:4)=QuaternionSupporter2&
                *Sqrt((1d0-dot_product(QuaternionSupporter1,QuaternionSupporter1))&
                /dot_product(QuaternionSupporter2,QuaternionSupporter2))
        mtp(r(1)).mol(chosenr(1)).Coordinate=struct(DisplacementVector,RotationQuaternion,&
            MoleculeInput(r(1)).RefConfig,MoleculeInput(r(1)).NAtoms)
        !z direction is confined by the slit wall, check if the atom is outside
        do i=1,MoleculeInput(r(1)).NAtoms
            if(abs(mtp(r(1)).mol(chosenr(1)).Coordinate(3,i))>HalfLength(3)) then
                accept=0
                !new=old
                NExist(p)=NExistOld(p)
                !Coordinate only, no need for single-body interaction
                mtp(p).mol(chosenp).Coordinate=mtpold(p).mol(chosenp).Coordinate
                return
            end if
        end do
        NExist(r(1))=NExist(r(1))+1
        del=del+EnergyChosenUpdated2D(r(1),chosenr(1))
        !Insert B randomly and uniformly
        call InsertMolecule2D(DisplacementVector,RotationQuaternion,.true.)
        mtp(r(2)).mol(chosenr(2)).Coordinate=struct(DisplacementVector,RotationQuaternion,&
            MoleculeInput(r(2)).RefConfig,MoleculeInput(r(2)).NAtoms)
        !z direction is confined by the slit wall, check if the atom is outside
        do i=1,MoleculeInput(r(2)).NAtoms
            if(abs(mtp(r(2)).mol(chosenr(2)).Coordinate(3,i))>HalfLength(3)) then
                accept=0
                !new=old
                NExist(r(1))=NExistOld(r(1))
                NExist(p)=NExistOld(p)
                !Coordinate only, no need for single-body interaction
                mtp(p).mol(chosenp).Coordinate=mtpold(p).mol(chosenp).Coordinate
                return
            end if
        end do
        NExist(r(2))=NExist(r(2))+1
        del=del+EnergyChosenUpdated2D(r(2),chosenr(2))-ReactionInput(rtp).deltaAig&
            -temperature*Log(Real(NExistOld(p)/NExist(r(1))/NExist(r(2))))
    end if
    !Accept or reject this trial
    accept=del<ImpossibleDel
    if(accept.and.del>0d0) then
        call random_number(rnd)
        accept=rnd<exp(-del/temperature)
    end if
    if(accept) then
        if(direction<0.5d0) then
            NTotal=NTotal-1
            !Complete flushing the chosen reactant molecules with the last exist reactant molecules
            do k=1,2
                !Flushing single-body interaction, old(r_k)=new(r_k) by the way
                mtp(r(k)).mol(chosenr(k)).Total=mtpold(r(k)).mol(NExistOld(r(k))).Total
                mtp(r(k)).mol(chosenr(k)).ESSBDependent2D=mtpold(r(k)).mol(NExistOld(r(k))).ESSBDependent2D
                mtp(r(k)).mol(chosenr(k)).WallLJ=mtpold(r(k)).mol(NExistOld(r(k))).WallLJ
                mtpold(r(k)).mol(chosenr(k))=mtp(r(k)).mol(chosenr(k))
                !Flushing many-body interaction, old(r_k)=new(r_k) by the way
                do j=1,r(k)-1
                    forall(i=1:NExistOld(j))
                        intertp(r(k),j).molab(chosenr(k),i)=intertpold(r(k),j).molab(NExistOld(r(k)),i)
                        intertpold(r(k),j).molab(chosenr(k),i)=intertp(r(k),j).molab(chosenr(k),i)
                    end forall
                end do
                forall(i=1:chosenr(k)-1)
                    intertp(r(k),r(k)).molab(chosenr(k),i)=intertpold(r(k),r(k)).molab(NExistOld(r(k)),i)
                    intertpold(r(k),r(k)).molab(chosenr(k),i)=intertp(r(k),r(k)).molab(chosenr(k),i)
                end forall
                !Same type molab(NExistOld(r),chosenr) is not needed as mol(NExistOld(r)) is deleted
                forall(i=chosenr(k)+1:NExist(r(k)))
                    intertp(r(k),r(k)).molab(i,chosenr(k))=intertpold(r(k),r(k)).molab(NExistOld(r(k)),i)
                    intertpold(r(k),r(k)).molab(i,chosenr(k))=intertp(r(k),r(k)).molab(i,chosenr(k))
                end forall
                do j=r(k)+1,MoleculeKinds
                    forall(i=1:NExistOld(j))
                        intertp(j,r(k)).molab(i,chosenr(k))=intertpold(j,r(k)).molab(i,NExistOld(r(k)))
                        intertpold(j,r(k)).molab(i,chosenr(k))=intertp(j,r(k)).molab(i,chosenr(k))
                    end forall
                end do
                !Complete old(r)=new(r)
                NExistOld(r(k))=NExist(r(k))
            end do
            !old(p)=new(p)
            NExistOld(p)=NExist(p)
            !Single-body
            mtpold(p).mol(chosenp)=mtp(p).mol(chosenp)
            !Many-body
            do j=1,p-1
                forall(i=1:NExist(j))
                    intertpold(p,j).molab(chosenp,i)=intertp(p,j).molab(chosenp,i)
                end forall
            end do
            !chosenp=NExist(p) so no need for the latter loop
            forall(i=1:chosenp-1)
                intertpold(p,p).molab(chosenp,i)=intertp(p,p).molab(chosenp,i)
            end forall
            do j=p+1,MoleculeKinds
                forall(i=1:NExist(j))
                    intertpold(j,p).molab(i,chosenp)=intertp(j,p).molab(i,chosenp)
                end forall
            end do
        else
            NTotal=NTotal+1
            !Complete flushing the chosen product molecule with the last exist product molecule
            !Flushing single-body interaction, old(p)=new(p) by the way
            mtp(p).mol(chosenp).Total=mtpold(p).mol(NExistOld(p)).Total
            mtp(p).mol(chosenp).ESSBDependent2D=mtpold(p).mol(NExistOld(p)).ESSBDependent2D
            mtp(p).mol(chosenp).WallLJ=mtpold(p).mol(NExistOld(p)).WallLJ
            mtpold(p).mol(chosenp)=mtp(p).mol(chosenp)
            !Flushing many-body interaction, old(p)=new(p) by the way
            do j=1,p-1
                forall(i=1:NExistOld(j))
                    intertp(p,j).molab(chosenp,i)=intertpold(p,j).molab(NExistOld(p),i)
                    intertpold(p,j).molab(chosenp,i)=intertp(p,j).molab(chosenp,i)
                end forall
            end do
            forall(i=1:chosenp-1)
                intertp(p,p).molab(chosenp,i)=intertpold(p,p).molab(NExistOld(p),i)
                intertpold(p,p).molab(chosenp,i)=intertp(p,p).molab(chosenp,i)
            end forall
            forall(i=chosenp+1:NExistOld(p))
                intertp(p,p).molab(i,chosenp)=intertpold(p,p).molab(NExistOld(p),i)
                intertpold(p,p).molab(i,chosenp)=intertp(p,p).molab(i,chosenp)
            end forall
            do j=p+1,MoleculeKinds
                forall(i=1:NExistOld(j))
                    intertp(j,p).molab(i,chosenp)=intertpold(j,p).molab(i,NExistOld(p))
                    intertpold(j,p).molab(i,chosenp)=intertp(j,p).molab(i,chosenp)
                end forall
            end do
            !Complete old(p)=new(p)
            NExistOld(p)=NExist(p)
            !old(r)=new(r)
            do k=1,2
                NExistOld(r(k))=NExist(r(k))
                !Single-body
                mtpold(r(k)).mol(chosenr(k))=mtp(r(k)).mol(chosenr(k))
                !Many-body
                do j=1,r(k)-1
                    forall(i=1:NExist(j))
                        intertpold(r(k),j).molab(chosenr(k),i)=intertp(r(k),j).molab(chosenr(k),i)
                    end forall
                end do
                forall(i=1:chosenr(k)-1)
                    intertpold(r(k),r(k)).molab(chosenr(k),i)=intertp(r(k),r(k)).molab(chosenr(k),i)
                end forall
                forall(i=chosenr(k)+1:NExist(r(k)))
                    intertpold(r(k),r(k)).molab(i,chosenr(k))=intertp(r(k),r(k)).molab(i,chosenr(k))
                end forall
                do j=r(k)+1,MoleculeKinds
                    forall(i=1:NExist(j))
                        intertpold(j,r(k)).molab(i,chosenr(k))=intertp(j,r(k)).molab(i,chosenr(k))
                    end forall
                end do
            end do
        end if
    else
        NExist(r(1))=NExistOld(r(1))
        NExist(r(2))=NExistOld(r(2))
        NExist(p)=NExistOld(p)
        !Coordinate only, no need for single-body interaction
        if(direction<0.5d0) then
            forall(i=1:2)
                mtp(r(i)).mol(chosenr(i)).Coordinate=mtpold(r(i)).mol(chosenr(i)).Coordinate
            end forall
        else
            mtp(p).mol(chosenp).Coordinate=mtpold(p).mol(chosenp).Coordinate
        end if
        !No need for many-body interaction
    end if
end subroutine ReactABC2D
!For 3D
subroutine ReactABC3D(rtp,r,p,chosenr,chosenp,direction,accept)
    integer,intent(in)::rtp
    integer,dimension(2),intent(out)::r,chosenr
    integer,intent(out)::p,chosenp
    real*8,intent(out)::direction
    logical,intent(out)::accept
    integer::i,j,k,ii,jj
    real*8::rnd,del
    real*8,dimension(2)::choose
    !To make notation simpler and reduce addressing
    r=ReactionInput(rtp).ReactantSequence
    p=ReactionInput(rtp).ProductSequence(1)
    !Determine the trial to make
    call random_number(direction)
    !Forward
    if(direction<0.5d0) then
        !Check the availability of the reactant
        if(NExist(r(1))<1.or.NExist(r(2))<1) then
            accept=0
            return
        end if
        !Choose two reactants and delete them by flushing the chosen ones with the last exist ones
        !Here only flush the coordinate and the existance number
        !Single-body and many-body interactions will be flushed once accepted
        call random_number(choose)
        forall(i=1:2)
            chosenr(i)=ceiling(choose(i)*NExist(r(i)))
            mtp(r(i)).mol(chosenr(i)).Coordinate=mtpold(r(i)).mol(NExist(r(i))).Coordinate
            NExist(r(i))=NExist(r(i))-1
        end forall
        !Compute their energy
        del=-EnergyChosenCurrent3D(r(1),chosenr(1))-EnergyChosenCurrent3D(r(2),chosenr(2))&
            !By input rule, the smaller molecule has the smaller type number
            +intertpold(r(1),r(2)).molab(chosenr(1),chosenr(2)).Total
        !Insert product at the cavity left by reactant A
        NExist(p)=NExist(p)+1
        chosenp=NExist(p)
        !The position of the centre of mass of the larger reactant
        DisplacementVector=mtpold(r(1)).mol(chosenr(1)).Coordinate(:,1)-MoleculeInput(r(1)).RefConfig(:,1)
        !Form the unit quaternion of the product orientation
            do while(1)
                call random_number(QuaternionSupporter1)
                QuaternionSupporter1=QuaternionSupporter1*2d0-1d0
                if(dot_product(QuaternionSupporter1,QuaternionSupporter1)<1d0) then
                    exit
                end if
            end do
            do while(1)
                call random_number(QuaternionSupporter2)
                QuaternionSupporter2=QuaternionSupporter2*2d0-1d0
                if(dot_product(QuaternionSupporter2,QuaternionSupporter2)<1d0) then
                    exit
                end if
            end do
            RotationQuaternion(1:2)=QuaternionSupporter1
            RotationQuaternion(3:4)=QuaternionSupporter2&
                *Sqrt((1d0-dot_product(QuaternionSupporter1,QuaternionSupporter1))&
                /dot_product(QuaternionSupporter2,QuaternionSupporter2))
        mtp(p).mol(chosenp).Coordinate=struct(DisplacementVector,RotationQuaternion,&
            MoleculeInput(p).RefConfig,MoleculeInput(p).NAtoms)
        del=del+EnergyChosenUpdated3D(p,chosenp)+ReactionInput(rtp).deltaAig&
            -temperature*Log(Real(NExistOld(r(1))*NExistOld(r(2))/NExist(p)))
    !Backward
    else
        !Check the availability of the product
        if(NExist(p)<1) then
            accept=0
            return
        end if
        !Choose a product
        call random_number(rnd)
        chosenp=ceiling(rnd*NExist(p))
        !Delete it by flushing the chosen one with the last exist one
        !Here only flush the coordinate and the existance number
        !Single-body and many-body interactions will be flushed once accepted
        mtp(p).mol(chosenp).Coordinate=mtpold(p).mol(NExist(p)).Coordinate
        NExist(p)=NExist(p)-1
        !Compute its energy
        del=-EnergyChosenCurrent3D(p,chosenp)
        !Insert reactant A at the cavity left by product
        forall(i=1:2)
            chosenr(i)=NExist(r(i))+1
        end forall
        !The position of the centre of mass of the product
        DisplacementVector=mtpold(p).mol(chosenp).Coordinate(:,1)-MoleculeInput(p).RefConfig(:,1)
        !Form the unit quaternion of the reactant orientation
            do while(1)
                call random_number(QuaternionSupporter1)
                QuaternionSupporter1=QuaternionSupporter1*2d0-1d0
                if(dot_product(QuaternionSupporter1,QuaternionSupporter1)<1d0) then
                    exit
                end if
            end do
            do while(1)
                call random_number(QuaternionSupporter2)
                QuaternionSupporter2=QuaternionSupporter2*2d0-1d0
                if(dot_product(QuaternionSupporter2,QuaternionSupporter2)<1d0) then
                    exit
                end if
            end do
            RotationQuaternion(1:2)=QuaternionSupporter1
            RotationQuaternion(3:4)=QuaternionSupporter2&
                *Sqrt((1d0-dot_product(QuaternionSupporter1,QuaternionSupporter1))&
                /dot_product(QuaternionSupporter2,QuaternionSupporter2))
        mtp(r(1)).mol(chosenr(1)).Coordinate=struct(DisplacementVector,RotationQuaternion,&
            MoleculeInput(r(1)).RefConfig,MoleculeInput(r(1)).NAtoms)
        NExist(r(1))=NExist(r(1))+1
        del=del+EnergyChosenUpdated3D(r(1),chosenr(1))
        !Insert B randomly and uniformly
        call InsertMolecule3D(DisplacementVector,RotationQuaternion)
        mtp(r(2)).mol(chosenr(2)).Coordinate=struct(DisplacementVector,RotationQuaternion,&
            MoleculeInput(r(2)).RefConfig,MoleculeInput(r(2)).NAtoms)
        NExist(r(2))=NExist(r(2))+1
        del=del+EnergyChosenUpdated3D(r(2),chosenr(2))-ReactionInput(rtp).deltaAig&
            -temperature*Log(Real(NExistOld(p)/NExist(r(1))/NExist(r(2))))
    end if
    !Accept or reject this trial
    accept=del<ImpossibleDel
    if(accept.and.del>0d0) then
        call random_number(rnd)
        accept=rnd<exp(-del/temperature)
    end if
    if(accept) then
        if(direction<0.5d0) then
            NTotal=NTotal-1
            !Complete flushing the chosen reactant molecules with the last exist reactant molecules
            do k=1,2
                !Currently no position dependent single-body interaction for 3D, old(r_k)=new(r_k) only
                mtpold(r(k)).mol(chosenr(k)).Coordinate=mtp(r(k)).mol(chosenr(k)).Coordinate
                !Flushing many-body interaction, old(r_k)=new(r_k) by the way
                do j=1,r(k)-1
                    forall(i=1:NExistOld(j))
                        intertp(r(k),j).molab(chosenr(k),i)=intertpold(r(k),j).molab(NExistOld(r(k)),i)
                        intertpold(r(k),j).molab(chosenr(k),i)=intertp(r(k),j).molab(chosenr(k),i)
                    end forall
                end do
                forall(i=1:chosenr(k)-1)
                    intertp(r(k),r(k)).molab(chosenr(k),i)=intertpold(r(k),r(k)).molab(NExistOld(r(k)),i)
                    intertpold(r(k),r(k)).molab(chosenr(k),i)=intertp(r(k),r(k)).molab(chosenr(k),i)
                end forall
                !Same type molab(NExistOld(r),chosenr) is not needed as mol(NExistOld(r)) is deleted
                forall(i=chosenr(k)+1:NExist(r(k)))
                    intertp(r(k),r(k)).molab(i,chosenr(k))=intertpold(r(k),r(k)).molab(NExistOld(r(k)),i)
                    intertpold(r(k),r(k)).molab(i,chosenr(k))=intertp(r(k),r(k)).molab(i,chosenr(k))
                end forall
                do j=r(k)+1,MoleculeKinds
                    forall(i=1:NExistOld(j))
                        intertp(j,r(k)).molab(i,chosenr(k))=intertpold(j,r(k)).molab(i,NExistOld(r(k)))
                        intertpold(j,r(k)).molab(i,chosenr(k))=intertp(j,r(k)).molab(i,chosenr(k))
                    end forall
                end do
                !Complete old(r)=new(r)
                NExistOld(r(k))=NExist(r(k))
            end do
            !old(p)=new(p)
            NExistOld(p)=NExist(p)
            !Currently no position dependent single-body interaction for 3D, coordinate only
            mtpold(p).mol(chosenp).Coordinate=mtp(p).mol(chosenp).Coordinate
            !Many-body
            do j=1,p-1
                forall(i=1:NExist(j))
                    intertpold(p,j).molab(chosenp,i)=intertp(p,j).molab(chosenp,i)
                end forall
            end do
            !chosenp=NExist(p) so no need for the latter loop
            forall(i=1:chosenp-1)
                intertpold(p,p).molab(chosenp,i)=intertp(p,p).molab(chosenp,i)
            end forall
            do j=p+1,MoleculeKinds
                forall(i=1:NExist(j))
                    intertpold(j,p).molab(i,chosenp)=intertp(j,p).molab(i,chosenp)
                end forall
            end do
        else
            NTotal=NTotal+1
            !Complete flushing the chosen product molecule with the last exist product molecule
            !Currently no position dependent single-body interaction for 3D, old(p)=new(p) only
            mtpold(p).mol(chosenp).Coordinate=mtp(p).mol(chosenp).Coordinate
            !Flushing many-body interaction, old(p)=new(p) by the way
            do j=1,p-1
                forall(i=1:NExistOld(j))
                    intertp(p,j).molab(chosenp,i)=intertpold(p,j).molab(NExistOld(p),i)
                    intertpold(p,j).molab(chosenp,i)=intertp(p,j).molab(chosenp,i)
                end forall
            end do
            forall(i=1:chosenp-1)
                intertp(p,p).molab(chosenp,i)=intertpold(p,p).molab(NExistOld(p),i)
                intertpold(p,p).molab(chosenp,i)=intertp(p,p).molab(chosenp,i)
            end forall
            forall(i=chosenp+1:NExistOld(p))
                intertp(p,p).molab(i,chosenp)=intertpold(p,p).molab(NExistOld(p),i)
                intertpold(p,p).molab(i,chosenp)=intertp(p,p).molab(i,chosenp)
            end forall
            do j=p+1,MoleculeKinds
                forall(i=1:NExistOld(j))
                    intertp(j,p).molab(i,chosenp)=intertpold(j,p).molab(i,NExistOld(p))
                    intertpold(j,p).molab(i,chosenp)=intertp(j,p).molab(i,chosenp)
                end forall
            end do
            !Complete old(p)=new(p)
            NExistOld(p)=NExist(p)
            !old(r)=new(r)
            do k=1,2
                NExistOld(r(k))=NExist(r(k))
                !Currently no position dependent single-body interaction for 3D, coordinate only
                mtpold(r(k)).mol(chosenr(k)).Coordinate=mtp(r(k)).mol(chosenr(k)).Coordinate
                !Many-body
                do j=1,r(k)-1
                    forall(i=1:NExist(j))
                        intertpold(r(k),j).molab(chosenr(k),i)=intertp(r(k),j).molab(chosenr(k),i)
                    end forall
                end do
                forall(i=1:chosenr(k)-1)
                    intertpold(r(k),r(k)).molab(chosenr(k),i)=intertp(r(k),r(k)).molab(chosenr(k),i)
                end forall
                forall(i=chosenr(k)+1:NExist(r(k)))
                    intertpold(r(k),r(k)).molab(i,chosenr(k))=intertp(r(k),r(k)).molab(i,chosenr(k))
                end forall
                do j=r(k)+1,MoleculeKinds
                    forall(i=1:NExist(j))
                        intertpold(j,r(k)).molab(i,chosenr(k))=intertp(j,r(k)).molab(i,chosenr(k))
                    end forall
                end do
            end do
        end if
    else
        NExist(r(1))=NExistOld(r(1))
        NExist(r(2))=NExistOld(r(2))
        NExist(p)=NExistOld(p)
        !Coordinate only, no need for single-body interaction
        if(direction<0.5d0) then
            forall(i=1:2)
                mtp(r(i)).mol(chosenr(i)).Coordinate=mtpold(r(i)).mol(chosenr(i)).Coordinate
            end forall
        else
            mtp(p).mol(chosenp).Coordinate=mtpold(p).mol(chosenp).Coordinate
        end if
        !No need for many-body interaction
    end if
end subroutine ReactABC3D

end module Trial