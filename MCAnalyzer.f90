!Available analyzers: density (for 1D and 2D), pressure tensor (for 1D and 2D), chemical potential
!Bin order: the bucket at the centre of the box is the first, the bucket at the edge is the last
!
!Comments on pressure tensor:
!Only intermolecular interaction contributes to pressure
!Coulombic systems must adopt Harasima mechanical route to evaluate pressure tensor:
!    1D: to be published, it is hard to account for electro static contribution from the periodic images
!        along the bucket crossing part of the contour, because alpha_k (Gubbins 2018 JCP) is hard to solve
!    2D: Sonne 2005 JCP, the bucket crossing part of the contour has singular electro static contribution from the self-images
!We do not care about the normal pressure as it is much smaller than the tangential pressures
!and pushes mainly on the wall rather than the molecules, so the normal direction is chosen to be the bucket crossing direction
!Although now normal pressure cannot be evaluated directly, it can still be determined from hydrostatic equilibrium:
!    For 1D, p_r = p_normal and p_phi = d( r * dp_normal ) / dr (Yun thesis chapter 3 section 3.1)
!    For 2D, p_z = p_normal = bulk pressure (Sonne 2005 JCP)
module MCAnalyzer
    use General
    use MCBasis
    use ConfigurationalEnergy
    implicit none

!Parameter
    !Number of buckets to be used for counting
    integer::bucketnumbers=100

!Derived type
    !Store single-body virial, and which bucket does the molecule or the atom fit in
    !Example: type(MoleculeTypeVirial),allocatable,dimension(:)::mtpv
    !         mtpv(i).mol(j).atm(k).virial(l) is l-th direction virial of k-th atom in j-th type i molecule
    !         This atom fits in mtpv(i).mol(j).atm(k).bucket-th bucket. Note that for 2D +z and -z fit in a same bucket so dividing 2 is required
    !         The centre of mass of j-th type i molecule fits in mtpv(i).mol(j).bucket-th bucket
    type AtomVirial
        integer::bucket
        !No single-body virial is needed in this program
        !Wall LJ only contributes to p_normal
        !The position independent part of Ewald summation has no virial as it has no force
        !2D shape-dependent part of Ewald summation only contributes to p_normal
        !real*8,dimension(3)::virial=0d0
    end type AtomVirial
    type MoleculeVirial
        integer::bucket
        type(AtomVirial),allocatable,dimension(:)::atm!short for AToM
    end type MoleculeVirial
    type MoleculeTypeVirial
        type(MoleculeVirial),allocatable,dimension(:)::mol!short for MOLecule
    end type MoleculeTypeVirial

    !Store many-body virial
    !Example: type(InterTypeVirial),allocatable,dimension(:,:)::intertpv
    !         intertpv(i,j).molab(k,l).atmab(m,n).virial(o) is o-th direction virial
    !         between m-th atom in k-th type i molecule and n-th atom in l-th type j molecule
    !         To avoid double counting, i>=j, k>l when i=j
    !The many-body virial includes LJ and many-body part of Ewald summation
    type InterAtomVirial
        !For 1D, virial(1) is z direction as usual, while:
        !    virial(2) is the phi direction of the n-th atom in l-th type j molecule
        !    virial(3) is the phi direction of the m-th atom in k-th type i molecule
        !For 2D, we do not compute p_z, so virial(3) is omitted
        real*8,dimension(3)::virial=0d0
    end type InterAtomVirial
    type InterMolecularVirial
        type(InterAtomVirial),allocatable,dimension(:,:)::atmab!short for AToM A and atom B
    end type InterMolecularVirial
    type InterTypeVirial
        type(InterMolecularVirial),allocatable,dimension(:,:)::molab!short for MOLecule A and molecule B
    end type InterTypeVirial

!Analyzation global variable
    !For 1D: radius of bucket, volume of central bucket
    !For 2D: length of bucket, total volume of +z bucket and -z bucket
    real*8::bucketunit,v_bucket
    real*8,allocatable,dimension(:)::bucket
    !Store virial, the existing molecules first
    type(MoleculeTypeVirial),allocatable,dimension(:)::mtpv!short for Molecule TyPe Virial
    type(InterTypeVirial),allocatable,dimension(:,:)::intertpv!short for INTER molecule TyPe Virial

contains
!The initializer for MCAnalyzer module
subroutine InitializeMCAnalyzer(ContinueAnalyzation)
    logical,intent(in)::ContinueAnalyzation
    integer::i,j,k,l
    !3D system does not need virial
    if(PeriodicDims<3) then
        !Initialize bucket
        if(PeriodicDims==1) then
            bucketunit=radius/bucketnumbers
            v_bucket=pi*bucketunit*bucketunit*Boxsize(1)
        else
            bucketunit=HalfLength(3)/bucketnumbers
            v_bucket=BoxSize(1)*BoxSize(2)*bucketunit*2d0
        end if
        allocate(bucket(bucketnumbers))
        forall(i=1:bucketnumbers)
            bucket(i)=bucketunit*(i-0.5d0)
        end forall
        bucket=bucket/AInAU
        !Initialize virial storage
        allocate(mtpv(MoleculeKinds))
        allocate(intertpv(MoleculeKinds,MoleculeKinds))
        do i=1,MoleculeKinds
            !Single-body
            allocate(mtpv(i).mol(MaxPossibleN(i)))
            do j=1,MaxPossibleN(i)
                allocate(mtpv(i).mol(j).atm(MoleculeInput(i).NAtoms))
            end do
            !Many-body
            !Virial between same type molecules
            allocate(intertpv(i,i).molab(MaxPossibleN(i),MaxPossibleN(i)))
            do k=1,MaxPossibleN(i)-1
                do l=k+1,MaxPossibleN(i)
                    allocate(intertpv(i,i).molab(l,k).atmab(MoleculeInput(i).NAtoms,MoleculeInput(i).NAtoms))
                end do
            end do
            !Virial between i-th and larger type molecules
            do j=i+1,MoleculeKinds
                allocate(intertpv(j,i).molab(MaxPossibleN(j),MaxPossibleN(i)))
                do k=1,MaxPossibleN(i)
                    do l=1,MaxPossibleN(j)
                        allocate(intertpv(j,i).molab(l,k).atmab(MoleculeInput(j).NAtoms,MoleculeInput(i).NAtoms))
                    end do
                end do
            end do
        end do
    end if
    !If this is not the 1st analyzation, read analyzation check point files
    if(ContinueAnalyzation) then
        if(PeriodicDims<3) then
            call ReadVirialCheckPoint()
        end if
    !Analyzation requires a higher accuracy level than optimization, so cannot inherit Optimize interactions
    else
        select case(PeriodicDims)
            case(1)
                AverageEnergy=EnergyTotalAverageUpdated1D()
                call UpdateTotalVirial1D()
            case(2)
                AverageEnergy=EnergyTotalAverageUpdated2D()
                call UpdateTotalVirial2D()
            case(3)
                AverageEnergy=EnergyTotalAverageUpdated3D()
        end select
        AverageGroundEnergy=AverageEnergy
        mtpold=mtp
        intertpold=intertp
        call WriteCheckPoint()
    end if
end subroutine InitializeMCAnalyzer

!Virial check point files writer, for 1D or 2D systems only
subroutine WriteVirialCheckPoint()
    integer::i,j,ii,jj,iatm,iiatm
    !Write the current single-body virial
    open(unit=99,file='SingleBodyVirial.CheckPoint',status='replace')
        do j=1,MoleculeKinds
            do i=1,NExist(j)
                write(99,*)mtpv(j).mol(i).bucket
                do iatm=1,MoleculeInput(j).NAtoms
                    write(99,*)mtpv(j).mol(i).atm(iatm).bucket
                end do
            end do
        end do
    close(99)
    !Write the current many-body virial
    open(unit=99,file='ManyBodyVirial.CheckPoint',status='replace')
        do j=1,MoleculeKinds
            do i=1,NExist(j)
                do ii=i+1,NExist(j)
                    do iatm=1,MoleculeInput(j).NAtoms
                        do iiatm=1,MoleculeInput(j).NAtoms
                            write(99,*)intertpv(j,j).molab(ii,i).atmab(iiatm,iatm).virial
                        end do
                    end do
                end do
                do jj=j+1,MoleculeKinds
                    do ii=1,NExist(jj)
                        do iatm=1,MoleculeInput(j).NAtoms
                            do iiatm=1,MoleculeInput(j).NAtoms
                                write(99,*)intertpv(jj,j).molab(ii,i).atmab(iiatm,iatm).virial
                            end do
                        end do
                    end do
                end do
            end do
        end do
    close(99)
end subroutine WriteVirialCheckPoint

!Virial check point files reader
subroutine ReadVirialCheckPoint()
    integer::i,j,ii,jj,iatm,iiatm
    !Read single-body virial
    open(unit=99,file='SingleBodyVirial.CheckPoint',status='replace')
        do j=1,MoleculeKinds
            do i=1,NExist(j)
                write(99,*)mtpv(j).mol(i).bucket
                do iatm=1,MoleculeInput(j).NAtoms
                    write(99,*)mtpv(j).mol(i).atm(iatm).bucket
                end do
            end do
        end do
    close(99)
    !Read many-body virial
    open(unit=99,file='ManyBodyVirial.CheckPoint',status='replace')
        do j=1,MoleculeKinds
            do i=1,NExist(j)
                do ii=i+1,NExist(j)
                    do iatm=1,MoleculeInput(j).NAtoms
                        do iiatm=1,MoleculeInput(j).NAtoms
                            read(99,*)intertpv(j,j).molab(ii,i).atmab(iiatm,iatm).virial
                        end do
                    end do
                end do
                do jj=j+1,MoleculeKinds
                    do ii=1,NExist(jj)
                        do iatm=1,MoleculeInput(j).NAtoms
                            do iiatm=1,MoleculeInput(j).NAtoms
                                read(99,*)intertpv(jj,j).molab(ii,i).atmab(iiatm,iatm).virial
                            end do
                        end do
                    end do
                end do
            end do
        end do
    close(99)
end subroutine ReadVirialCheckPoint

!------------- Observable computer -------------
    !Count each bucket contains how many molecules
    function MoleculeNumberInEachBin()
        integer,dimension(bucketnumbers,MoleculeKinds)::MoleculeNumberInEachBin
        integer::i,j
        !Initialization: set the counter to zero
        MoleculeNumberInEachBin=0
        !Loop over molecule types
        do i=1,MoleculeKinds
            !Loop over existing type i molecules
            do j=1,NExist(i)
                MoleculeNumberInEachBin(mtpv(i).mol(j).bucket,i)=MoleculeNumberInEachBin(mtpv(i).mol(j).bucket,i)+1
            end do
        end do
    end function MoleculeNumberInEachBin
    
    !Compute density from molecule number, then output
    !TotalDensity is meant for kinetic (ideal gas) pressure calculation later
    !For 1D
    subroutine MoleculeNumber2Density1D(MoleculeNumber,TotalDensity)
        real*8,dimension(bucketnumbers,MoleculeKinds)::MoleculeNumber
        real*8,dimension(bucketnumbers),intent(out)::TotalDensity
        integer::i,j
        forall(i=1:bucketnumbers)
            MoleculeNumber(i,:)=MoleculeNumber(i,:)/(2*i-1)/v_bucket
        end forall
        TotalDensity=0d0
        do i=1,MoleculeKinds
            TotalDensity=TotalDensity+MoleculeNumber(:,i)
        end do
        !Convert to human unit
        MoleculeNumber=MoleculeNumber*(AInAU*AInAU*AInAU*1d27)/NAvogadro
        open(unit=99,file='Density.txt',status='replace')
            do i=1,MoleculeKinds
                write(99,'(A8,I3)')'Molecule',i
                write(99,'(A43)')'DistanceFromCentre/A'//char(9)//'Concentration/mol*L^-1'
                do j=1,bucketnumbers
                    write(99,'(F10.5,A1,F20.10)')bucket(j),char(9),MoleculeNumber(j,i)
                end do
                write(99,*)
            end do
        close(99)
    end subroutine MoleculeNumber2Density1D
    !For 2D
    subroutine MoleculeNumber2Density2D(MoleculeNumber,TotalDensity)
        real*8,dimension(bucketnumbers,MoleculeKinds)::MoleculeNumber
        real*8,dimension(bucketnumbers),intent(out)::TotalDensity
        integer::i,j
        MoleculeNumber=MoleculeNumber/v_bucket
        TotalDensity=0d0
        do i=1,MoleculeKinds
            TotalDensity=TotalDensity+MoleculeNumber(:,i)
        end do
        !Convert to human unit: mol/L
        MoleculeNumber=MoleculeNumber*(AInAU*AInAU*AInAU*1d27)/NAvogadro
        open(unit=99,file='Density.txt',status='replace')
            do i=1,MoleculeKinds
                write(99,'(A8,I3)')'Molecule',i
                write(99,'(A43)')'DistanceFromCentre/A'//char(9)//'Concentration/mol*L^-1'
                do j=1,bucketnumbers
                    write(99,'(F10.5,A1,F20.10)')bucket(j),char(9),MoleculeNumber(j,i)
                end do
                write(99,*)
            end do
        close(99)
    end subroutine MoleculeNumber2Density2D

    !Count [Virial_z,Virial_phi]
    function Virial1DInEachBin()
        real*8,dimension(2,bucketnumbers)::Virial1DInEachBin
        integer::i,j,k,l,m,n
        !Initialization: set the counter to zero
        Virial1DInEachBin=0d0
        !Loop over molecule types
        do j=1,MoleculeKinds
            !Single-body: No single-body virial is needed in this program
            !Many-body
            !Same type
            do l=1,NExist(j)-1
                do k=l+1,NExist(j)
                    do n=1,MoleculeInput(j).NAtoms
                        do m=1,MoleculeInput(j).NAtoms
                            Virial1DInEachBin(:,mtpv(j).mol(l).atm(n).bucket)=Virial1DInEachBin(:,mtpv(j).mol(l).atm(n).bucket)&
                                +intertpv(j,j).molab(k,l).atmab(m,n).virial(1:2)
                            Virial1DInEachBin(1,mtpv(j).mol(k).atm(m).bucket)=Virial1DInEachBin(1,mtpv(j).mol(k).atm(m).bucket)&
                                +intertpv(j,j).molab(k,l).atmab(m,n).virial(1)
                            Virial1DInEachBin(2,mtpv(j).mol(k).atm(m).bucket)=Virial1DInEachBin(2,mtpv(j).mol(k).atm(m).bucket)&
                                +intertpv(j,j).molab(k,l).atmab(m,n).virial(3)
                        end do
                    end do
                end do
            end do
            !Different type
            do i=j+1,MoleculeKinds
                do l=1,NExist(j)
                    do k=1,NExist(i)
                        do n=1,MoleculeInput(j).NAtoms
                            do m=1,MoleculeInput(i).NAtoms
                                Virial1DInEachBin(:,mtpv(j).mol(l).atm(n).bucket)=Virial1DInEachBin(:,mtpv(j).mol(l).atm(n).bucket)&
                                    +intertpv(i,j).molab(k,l).atmab(m,n).virial(1:2)
                                Virial1DInEachBin(1,mtpv(i).mol(k).atm(m).bucket)=Virial1DInEachBin(1,mtpv(i).mol(k).atm(m).bucket)&
                                    +intertpv(i,j).molab(k,l).atmab(m,n).virial(1)
                                Virial1DInEachBin(2,mtpv(i).mol(k).atm(m).bucket)=Virial1DInEachBin(2,mtpv(i).mol(k).atm(m).bucket)&
                                    +intertpv(i,j).molab(k,l).atmab(m,n).virial(3)
                            end do
                        end do
                    end do
                end do
            end do
        end do
    end function Virial1DInEachBin
    !Count Virial_x,+z + Virial_y,+z + Virial_x,-z + Virial_y,-z
    function Virial2DInEachBin()
        real*8,dimension(bucketnumbers)::Virial2DInEachBin
        integer::i,j,k,l,m,n
        !Initialization: set the counter to zero
        Virial2DInEachBin=0d0
        !Loop over molecule types
        do i=1,MoleculeKinds
            !Single-body: No single-body virial is needed in this program
            !Many-body
            !Same type
            do k=1,NExist(i)-1
                do l=k+1,NExist(i)
                    do m=1,MoleculeInput(i).NAtoms
                        do n=1,MoleculeInput(i).NAtoms
                            Virial2DInEachBin(mtpv(i).mol(k).atm(m).bucket)=Virial2DInEachBin(mtpv(i).mol(k).atm(m).bucket)&
                                +intertpv(i,i).molab(l,k).atmab(n,m).virial(1)+intertpv(i,i).molab(l,k).atmab(n,m).virial(2)
                            Virial2DInEachBin(mtpv(i).mol(l).atm(n).bucket)=Virial2DInEachBin(mtpv(i).mol(l).atm(n).bucket)&
                                +intertpv(i,i).molab(l,k).atmab(n,m).virial(1)+intertpv(i,i).molab(l,k).atmab(n,m).virial(2)
                        end do
                    end do
                end do
            end do
            !Different type
            do j=i+1,MoleculeKinds
                do k=1,NExist(i)
                    do l=1,NExist(j)
                        do m=1,MoleculeInput(i).NAtoms
                            do n=1,MoleculeInput(j).NAtoms
                                Virial2DInEachBin(mtpv(i).mol(k).atm(m).bucket)=Virial2DInEachBin(mtpv(i).mol(k).atm(m).bucket)&
                                    +intertpv(j,i).molab(l,k).atmab(n,m).virial(1)+intertpv(j,i).molab(l,k).atmab(n,m).virial(2)
                                Virial2DInEachBin(mtpv(j).mol(l).atm(n).bucket)=Virial2DInEachBin(mtpv(j).mol(l).atm(n).bucket)&
                                    +intertpv(j,i).molab(l,k).atmab(n,m).virial(1)+intertpv(j,i).molab(l,k).atmab(n,m).virial(2)
                            end do
                        end do
                    end do
                end do
            end do
        end do
    end function Virial2DInEachBin

    !Compute pressure tensor from virial and density, then output
    !For 1D
    subroutine Virial2Pressure1D(Virial,TotalDensity)
        real*8,dimension(2,bucketnumbers)::Virial
        real*8,dimension(bucketnumbers),intent(in)::TotalDensity
        integer::i
        forall(i=1:bucketnumbers)
            Virial(:,i)=Virial(:,i)/2d0/(2*i-1)/v_bucket
        end forall
        Virial(1,:)=TotalDensity*temperature+Virial(1,:)
        Virial(2,:)=TotalDensity*temperature+Virial(2,:)
        !Convert to human unit
        Virial=Virial/barInAU
        open(unit=99,file='PressureTensor.txt',status='replace')
            write(99,'(A40)')'DistanceFromCentre/A'//char(9)//'p_z/bar'//char(9)//'p_theta/bar'
            do i=1,bucketnumbers
                write(99,'(F10.5,A1,F20.10,A1,F20.10)')bucket(i),char(9),Virial(1,i),char(9),Virial(2,i)
            end do
        close(99)
    end subroutine Virial2Pressure1D
    !For 2D
    subroutine Virial2Pressure2D(Virial,TotalDensity)
        real*8,dimension(bucketnumbers)::Virial
        real*8,dimension(bucketnumbers),intent(in)::TotalDensity
        integer::i
        Virial=TotalDensity*temperature+Virial/4d0/v_bucket
        !Convert to human unit
        Virial=Virial/barInAU
        open(unit=99,file='PressureTensor.txt',status='replace')
            write(99,'(A37)')'DistanceFromCentre/A'//char(9)//'p_tangential/bar'
            do i=1,bucketnumbers
                write(99,'(F10.5,A1,F20.10)')bucket(i),char(9),Virial(i)
            end do
        close(99)
    end subroutine Virial2Pressure2D

    !Widom particle insertion for computing chemical potential：
    !miu = - k_b T ln( <Widom> * / rho / lamda^3 ), where rho = number density
    !For 1D
    function Widom1D()
        real*8,dimension(MoleculeKinds)::Widom1D
        integer::i,j,n
        do i=1,MoleculeKinds
            !Insert a virtual molecule
            n=NExist(i)+1
            call InsertMolecule1D(DisplacementVector,RotationQuaternion,.true.)
            mtp(i).mol(n).Coordinate=struct(DisplacementVector,RotationQuaternion,&
                        MoleculeInput(i).RefConfig,MoleculeInput(i).NAtoms)
            Widom1D(i)=exp(-EnergyChosenUpdated1D(i,n)/temperature)
        end do
    end function Widom1D
    !For 2D
    function Widom2D()
        real*8,dimension(MoleculeKinds)::Widom2D
        integer::i
        do i=1,MoleculeKinds
            !Insert a virtual molecule
            call InsertMolecule2D(DisplacementVector,RotationQuaternion,.true.)
            mtp(i).mol(NExist(i)+1).Coordinate=struct(DisplacementVector,RotationQuaternion,&
                MoleculeInput(i).RefConfig,MoleculeInput(i).NAtoms)
            Widom2D(i)=exp(-EnergyChosenUpdated1D(i,NExist(i)+1)/temperature)
        end do
    end function Widom2D
    !For 3D
    function Widom3D()
        real*8,dimension(MoleculeKinds)::Widom3D
        integer::i
        do i=1,MoleculeKinds
            !Insert a virtual molecule
            call InsertMolecule3D(DisplacementVector,RotationQuaternion)
            mtp(i).mol(NExist(i)+1).Coordinate=struct(DisplacementVector,RotationQuaternion,&
                MoleculeInput(i).RefConfig,MoleculeInput(i).NAtoms)
            Widom3D(i)=exp(-EnergyChosenUpdated1D(i,NExist(i)+1)/temperature)
        end do
    end function Widom3D
!--------------------- End ---------------------

!--------------- Virial updater ----------------
    !Update the chosen molecule's virial
    !For 1D
    subroutine UpdateChosenVirial1D(tp,chosen)
        integer,intent(in)::tp,chosen
        integer::i,j
        !Single-body
        DisplacementVector(1:2)=mtp(tp).mol(chosen).Coordinate(2:3,1)-MoleculeInput(tp).RefConfig(2:3,1)
        mtpv(tp).mol(chosen).bucket=Ceiling(Norm2(DisplacementVector(1:2))/bucketunit)
        forall(i=1:MoleculeInput(tp).NAtoms)
            mtpv(tp).mol(chosen).atm(i).bucket=Ceiling(Norm2(mtp(tp).mol(chosen).Coordinate(2:3,i))/bucketunit)
        end forall
        !Many body
        !Virial between smaller type molecules and chosen molecule
        do j=1,tp-1
            do i=1,NExist(j)
                call UpdateNonbondedVirial1D(tp,j,chosen,i)
            end do
        end do
        !Virial between same type smaller order molecules and chosen molecule
        do i=1,chosen-1
            call UpdateNonbondedVirial1D(tp,tp,chosen,i)
        end do
        !Virial between same type larger order molecules and chosen molecule
        do i=chosen+1,NExist(tp)
            call UpdateNonbondedVirial1D(tp,tp,i,chosen)
        end do
        !Virial between larger type molecules and chosen molecule
        do j=tp+1,MoleculeKinds
            do i=1,NExist(j)
                call UpdateNonbondedVirial1D(j,tp,i,chosen)
            end do
        end do
    end subroutine UpdateChosenVirial1D
    !For 2D
    subroutine UpdateChosenVirial2D(tp,chosen)
        integer,intent(in)::tp,chosen
        integer::i,j
        !Single-body
        mtpv(tp).mol(chosen).bucket=Ceiling(Abs(mtp(tp).mol(chosen).Coordinate(3,1)-MoleculeInput(tp).RefConfig(3,1))/bucketunit)
        forall(i=1:MoleculeInput(tp).NAtoms)
            mtpv(tp).mol(chosen).atm(i).bucket=Ceiling(Abs(mtp(tp).mol(chosen).Coordinate(3,i))/bucketunit)
        end forall
        !Many body
        !Virial between smaller type molecules and chosen molecule
        do j=1,tp-1
            do i=1,NExist(j)
                call UpdateNonbondedVirial2D(tp,j,chosen,i)
            end do
        end do
        !Virial between same type smaller order molecules and chosen molecule
        do i=1,chosen-1
            call UpdateNonbondedVirial2D(tp,tp,chosen,i)
        end do
        !Virial between same type larger order molecules and chosen molecule
        do i=chosen+1,NExist(tp)
            call UpdateNonbondedVirial2D(tp,tp,i,chosen)
        end do
        !Virial between larger type molecules and chosen molecule
        do j=tp+1,MoleculeKinds
            do i=1,NExist(j)
                call UpdateNonbondedVirial2D(j,tp,i,chosen)
            end do
        end do
    end subroutine UpdateChosenVirial2D

    !Update all molecules' virial
    !For 1D
    subroutine UpdateTotalVirial1D()
        integer::i,j,ii,jj
        do j=1,MoleculeKinds
            do i=1,NExist(j)
                !Single-body
                DisplacementVector(1:2)=mtp(j).mol(i).Coordinate(2:3,1)-MoleculeInput(j).RefConfig(2:3,1)
                mtpv(j).mol(i).bucket=Ceiling(Norm2(DisplacementVector(1:2))/bucketunit)
                forall(ii=1:MoleculeInput(j).NAtoms)
                    mtpv(j).mol(i).atm(ii).bucket=Ceiling(Norm2(mtp(j).mol(i).Coordinate(2:3,ii))/bucketunit)
                end forall
                !Many-body
                do ii=i+1,NExist(j)
                    call UpdateNonbondedVirial1D(j,j,ii,i)
                end do
                do jj=j+1,MoleculeKinds
                    do ii=1,NExist(jj)
                        call UpdateNonbondedVirial1D(jj,j,ii,i)
                    end do
                end do
            end do
        end do
    end subroutine UpdateTotalVirial1D
    !For 2D
    subroutine UpdateTotalVirial2D()
        integer::i,j,ii,jj
        do j=1,MoleculeKinds
            do i=1,NExist(j)
                !Single-body
                mtpv(j).mol(i).bucket=Ceiling(Abs(mtp(j).mol(i).Coordinate(3,1)-MoleculeInput(j).RefConfig(3,1))/bucketunit)
                forall(ii=1:MoleculeInput(j).NAtoms)
                    mtpv(j).mol(i).atm(ii).bucket=Ceiling(Abs(mtp(j).mol(i).Coordinate(3,ii))/bucketunit)
                end forall
                !Many-body
                do ii=i+1,NExist(j)
                    call UpdateNonbondedVirial2D(j,j,ii,i)
                end do
                do jj=j+1,MoleculeKinds
                    do ii=1,NExist(jj)
                        call UpdateNonbondedVirial2D(jj,j,ii,i)
                    end do
                end do
            end do
        end do
    end subroutine UpdateTotalVirial2D

    !Update virial between molecule 1 and molecule 2
    !For 1D
    subroutine UpdateNonbondedVirial1D(tp1,tp2,mol1,mol2)
        integer,intent(in)::tp1,tp2,mol1,mol2
        integer::i,j,iclass,jclass
        !Class: 1 epsilon only and 2 non-zero epsilon
        do iclass=1,MoleculeInput(tp1).NEpsilonOnly
            i=MoleculeInput(tp1).EpsilonOnly(iclass)
            do jclass=1,MoleculeInput(tp2).NNon0Epsilon
                j=MoleculeInput(tp2).Non0Epsilon(jclass)
                intertpv(tp1,tp2).molab(mol1,mol2).atmab(i,j).virial=&
                    VanDeWaalsVirial1D(tp1,tp2,i,j,mtp(tp1).mol(mol1).Coordinate(:,i),mtp(tp2).mol(mol2).Coordinate(:,j))
            end do
        end do
        !Class: 1 charge only and 2 non-zero charge
        do iclass=1,MoleculeInput(tp1).NChargeOnly
            i=MoleculeInput(tp1).ChargeOnly(iclass)
            do jclass=1,MoleculeInput(tp2).NNon0Charge
                j=MoleculeInput(tp2).Non0Charge(jclass)
                intertpv(tp1,tp2).molab(mol1,mol2).atmab(i,j).virial=&
                    CoulombVirial1D(tp1,tp2,i,j,mtp(tp1).mol(mol1).Coordinate(:,i),mtp(tp2).mol(mol2).Coordinate(:,j))
            end do
        end do
        !Class: 1 both and...
        do iclass=1,MoleculeInput(tp1).NBoth
            i=MoleculeInput(tp1).Both(iclass)
            !2 epsilon only
            do jclass=1,MoleculeInput(tp2).NEpsilonOnly
                j=MoleculeInput(tp2).EpsilonOnly(jclass)
                intertpv(tp1,tp2).molab(mol1,mol2).atmab(i,j).virial=&
                    VanDeWaalsVirial1D(tp1,tp2,i,j,mtp(tp1).mol(mol1).Coordinate(:,i),mtp(tp2).mol(mol2).Coordinate(:,j))
            end do
            !2 charge only
            do jclass=1,MoleculeInput(tp2).NChargeOnly
                j=MoleculeInput(tp2).ChargeOnly(jclass)
                intertpv(tp1,tp2).molab(mol1,mol2).atmab(i,j).virial=&
                    CoulombVirial1D(tp1,tp2,i,j,mtp(tp1).mol(mol1).Coordinate(:,i),mtp(tp2).mol(mol2).Coordinate(:,j))
            end do
            !2 both
            do jclass=1,MoleculeInput(tp2).NBoth
                j=MoleculeInput(tp2).Both(jclass)
                intertpv(tp1,tp2).molab(mol1,mol2).atmab(i,j).virial=&
                    NonbondedVirial1D(tp1,tp2,i,j,mtp(tp1).mol(mol1).Coordinate(:,i),mtp(tp2).mol(mol2).Coordinate(:,j))
            end do
        end do
    end subroutine UpdateNonbondedVirial1D
    !For 2D
    subroutine UpdateNonbondedVirial2D(tp1,tp2,mol1,mol2)
        integer,intent(in)::tp1,tp2,mol1,mol2
        integer::i,j,iclass,jclass
        !Class: 1 epsilon only and 2 non-zero epsilon
        do iclass=1,MoleculeInput(tp1).NEpsilonOnly
            i=MoleculeInput(tp1).EpsilonOnly(iclass)
            do jclass=1,MoleculeInput(tp2).NNon0Epsilon
                j=MoleculeInput(tp2).Non0Epsilon(jclass)
                intertpv(tp1,tp2).molab(mol1,mol2).atmab(i,j).virial(1:2)=&
                    VanDeWaalsVirial2D(tp1,tp2,i,j,mtp(tp1).mol(mol1).Coordinate(:,i),mtp(tp2).mol(mol2).Coordinate(:,j))
            end do
        end do
        !Class: 1 charge only and 2 non-zero charge
        do iclass=1,MoleculeInput(tp1).NChargeOnly
            i=MoleculeInput(tp1).ChargeOnly(iclass)
            do jclass=1,MoleculeInput(tp2).NNon0Charge
                j=MoleculeInput(tp2).Non0Charge(jclass)
                intertpv(tp1,tp2).molab(mol1,mol2).atmab(i,j).virial(1:2)=&
                    CoulombVirial2D(tp1,tp2,i,j,mtp(tp1).mol(mol1).Coordinate(:,i),mtp(tp2).mol(mol2).Coordinate(:,j))
            end do
        end do
        !Class: 1 both and...
        do iclass=1,MoleculeInput(tp1).NBoth
            i=MoleculeInput(tp1).Both(iclass)
            !2 epsilon only
            do jclass=1,MoleculeInput(tp2).NEpsilonOnly
                j=MoleculeInput(tp2).EpsilonOnly(jclass)
                intertpv(tp1,tp2).molab(mol1,mol2).atmab(i,j).virial(1:2)=&
                    VanDeWaalsVirial2D(tp1,tp2,i,j,mtp(tp1).mol(mol1).Coordinate(:,i),mtp(tp2).mol(mol2).Coordinate(:,j))
            end do
            !2 charge only
            do jclass=1,MoleculeInput(tp2).NChargeOnly
                j=MoleculeInput(tp2).ChargeOnly(jclass)
                intertpv(tp1,tp2).molab(mol1,mol2).atmab(i,j).virial(1:2)=&
                    CoulombVirial2D(tp1,tp2,i,j,mtp(tp1).mol(mol1).Coordinate(:,i),mtp(tp2).mol(mol2).Coordinate(:,j))
            end do
            !2 both
            do jclass=1,MoleculeInput(tp2).NBoth
                j=MoleculeInput(tp2).Both(jclass)
                intertpv(tp1,tp2).molab(mol1,mol2).atmab(i,j).virial(1:2)=&
                    NonbondedVirial2D(tp1,tp2,i,j,mtp(tp1).mol(mol1).Coordinate(:,i),mtp(tp2).mol(mol2).Coordinate(:,j))
            end do
        end do
    end subroutine UpdateNonbondedVirial2D

    !The Van de Waals virial between force centres 1 and 2
    !For 1D
    function VanDeWaalsVirial1D(tp1,tp2,atm1,atm2,co1,co2)
        !VanDeWaalsVirial1D(2) is the phi direction virial of force centre 2
        !VanDeWaalsVirial1D(3) is the phi direction virial of force centre 1
        real*8,dimension(3)::VanDeWaalsVirial1D
        integer,intent(in)::tp1,tp2,atm1,atm2
        real*8,dimension(3),intent(in)::co1,co2
        real*8::r2,r6,r12
        real*8,dimension(3)::f
        !From 1 point to 2
        DisplacementVector=co2-co1
        !Periodic boundary condition
        if(DisplacementVector(1)>HalfLength(1)) then
            DisplacementVector(1)=DisplacementVector(1)-BoxSize(1)
        else if(DisplacementVector(1)<-HalfLength(1)) then
            DisplacementVector(1)=DisplacementVector(1)+BoxSize(1)
        end if
        r2=dot_product(DisplacementVector,DisplacementVector)
        r6=intertpLJP(tp1,tp2).atmab(atm1,atm2).sigmasq/r2
        r6=r6*r6*r6
        r12=r6*r6
        f=intertpLJP(tp1,tp2).atmab(atm1,atm2).epsilonm4*6d0*(2d0*r12-r6)*DisplacementVector/r2
        VanDeWaalsVirial1D(1)=DisplacementVector(1)*f(1)
        !Temporarily r2 = Virial^12_inplane
        r2=dot_product(DisplacementVector(2:3),f(2:3))
        !Temporarily r12 = rho_1
        r12=Norm2(co1(2:3))
        !Temporarily r6 = rho_2
        r6=Norm2(co2(2:3))
        !Temporarily f = f(2:3) * ( rho_2 - rho_1 )
        f(2:3)=f(2:3)*(r6-r12)
        VanDeWaalsVirial1D(2)=r2-dot_product(f(2:3),co1(2:3)/r12)
        VanDeWaalsVirial1D(3)=r2-dot_product(f(2:3),co2(2:3)/r6)
    end function VanDeWaalsVirial1D
    !For 2D
    function VanDeWaalsVirial2D(tp1,tp2,atm1,atm2,co1,co2)
        real*8,dimension(2)::VanDeWaalsVirial2D
        integer,intent(in)::tp1,tp2,atm1,atm2
        real*8,dimension(3),intent(in)::co1,co2
        integer::k
        real*8::r2,r6,r12
        real*8,dimension(2)::runit
        !From 1 point to 2
        DisplacementVector=co2-co1
        !Periodic boundary condition
        do k=1,2
            if(DisplacementVector(k)>HalfLength(k)) then
                DisplacementVector(k)=DisplacementVector(k)-BoxSize(k)
            else if(DisplacementVector(k)<-HalfLength(k)) then
                DisplacementVector(k)=DisplacementVector(k)+BoxSize(k)
            end if
        end do
        r2=dot_product(DisplacementVector,DisplacementVector)
        r6=intertpLJP(tp1,tp2).atmab(atm1,atm2).sigmasq/r2
        r6=r6*r6*r6
        r12=r6*r6
        VanDeWaalsVirial2D=DisplacementVector(1:2)*DisplacementVector(1:2)/r2*&
            intertpLJP(tp1,tp2).atmab(atm1,atm2).epsilonm4*6d0*(2d0*r12-r6)
    end function VanDeWaalsVirial2D
!Need update
    !The Coulomb virial between force centres 1 and 2
    !For 1D
    function CoulombVirial1D(tp1,tp2,atm1,atm2,co1,co2)
        !CoulombVirial1D(2) is the phi direction virial of force centre 2
        !CoulombVirial1D(3) is the phi direction virial of force centre 1
        real*8,dimension(3)::CoulombVirial1D
        integer,intent(in)::tp1,tp2,atm1,atm2
        real*8,dimension(3),intent(in)::co1,co2
        real*8::r,r2,chisqrhosq
        real*8,dimension(2)::rhounit
        real*8,dimension(3)::runit
        !From 1 point to 2
        DisplacementVector=co2-co1
        !Periodic boundary condition
        if(DisplacementVector(1)>HalfLength(1)) then
            DisplacementVector(1)=DisplacementVector(1)-BoxSize(1)
        else if(DisplacementVector(1)<-HalfLength(1)) then
            DisplacementVector(1)=DisplacementVector(1)+BoxSize(1)
        end if
        !Temporarily chisqrhosq = rho^2_12
        chisqrhosq=dot_product(DisplacementVector(2:3),DisplacementVector(2:3))
        rhounit=DisplacementVector(2:3)/Sqrt(chisqrhosq)
        r2=chisqrhosq+DisplacementVector(1)*DisplacementVector(1)
        !Now chisqrhosq is what it is meant to be
        chisqrhosq=chisqrhosq*chisqrhosq*chisq
        r=Sqrt(r2)
        runit=DisplacementVector/r
        !Temporarily runit = Coulomb force_12 - 1D_Gz0_12
        runit=MoleculeInput(tp1).charge(atm1)*MoleculeInput(tp2).charge(atm2)&
            *((Nonbonded_chim2dsqrtpi*Exp(-chisq*r2)/r+Erfc(chi*r)/r2)*runit&
            +pim4dV*ReciprocalForceSummation1D(DisplacementVector))
        !Temporarily runit = Coulomb force_12
        runit(2:3)=runit(2:3)-MoleculeInput(tp1).charge(atm1)*MoleculeInput(tp2).charge(atm2)&
            *(Nonbonded_chisqm2dLz*(1-Exp(-chisqrhosq))/chisqrhosq*rhounit)
        CoulombVirial1D(1)=DisplacementVector(1)*runit(1)
        !Temporarily r = Virial^12_inplane
        r=dot_product(DisplacementVector(2:3),runit(2:3))
        !Temporarily r2 = rho_1
        r2=Norm2(co1(2:3))
        !Temporarily chisqrhosq = rho_2
        chisqrhosq=Norm2(co2(2:3))
        !Temporarily runit = f_12(2:3) * ( rho_2 - rho_1 )
        runit(2:3)=runit(2:3)*(chisqrhosq-r2)
        CoulombVirial1D(2)=r-dot_product(runit(2:3),co1(2:3)/r2)
        CoulombVirial1D(3)=r-dot_product(runit(2:3),co2(2:3)/chisqrhosq)
    end function CoulombVirial1D
    !For 2D
    function CoulombVirial2D(tp1,tp2,atm1,atm2,co1,co2)
        real*8,dimension(2)::CoulombVirial2D
        integer,intent(in)::tp1,tp2,atm1,atm2
        real*8,dimension(3),intent(in)::co1,co2
        integer::k
        real*8::r,r2
        real*8,dimension(2)::runit
        !From 1 point to 2
        DisplacementVector=co2-co1
        !Periodic boundary condition
        do k=1,2
            if(DisplacementVector(k)>HalfLength(k)) then
                DisplacementVector(k)=DisplacementVector(k)-BoxSize(k)
            else if(DisplacementVector(k)<-HalfLength(k)) then
                DisplacementVector(k)=DisplacementVector(k)+BoxSize(k)
            end if
        end do
        r2=dot_product(DisplacementVector,DisplacementVector)
        r=Sqrt(r2)
        runit=DisplacementVector(1:2)/r
        CoulombVirial2D=DisplacementVector(1:2)*(&
            MoleculeInput(tp1).charge(atm1)*MoleculeInput(tp2).charge(atm2)&
            *((Nonbonded_chim2dsqrtpi*Exp(-chisq*r2)/r+Erfc(chi*r)/r2)*runit&
            +pim4dV*ReciprocalForceSummation(DisplacementVector)))
    end function CoulombVirial2D

    !The total nonbonded virial between force centres 1 and 2
    !For 1D
    function NonbondedVirial1D(tp1,tp2,atm1,atm2,co1,co2)
        !NonbondedVirial1D(2) is the phi direction virial of force centre 2
        !NonbondedVirial1D(3) is the phi direction virial of force centre 1
        real*8,dimension(3)::NonbondedVirial1D
        integer,intent(in)::tp1,tp2,atm1,atm2
        real*8,dimension(3),intent(in)::co1,co2
        real*8::r,r2,r6,r12,chisqrhosq
        real*8,dimension(2)::rhounit
        real*8,dimension(3)::runit
        !From 1 point to 2
        DisplacementVector=co2-co1
        !Periodic boundary condition
        if(DisplacementVector(1)>HalfLength(1)) then
            DisplacementVector(1)=DisplacementVector(1)-BoxSize(1)
        else if(DisplacementVector(1)<-HalfLength(1)) then
            DisplacementVector(1)=DisplacementVector(1)+BoxSize(1)
        end if
        !Temporarily chisqrhosq = rho^2_12
        chisqrhosq=dot_product(DisplacementVector(2:3),DisplacementVector(2:3))
        rhounit=DisplacementVector(2:3)/Sqrt(chisqrhosq)
        r2=chisqrhosq+DisplacementVector(1)*DisplacementVector(1)
        !Now chisqrhosq is what it is meant to be
        chisqrhosq=chisqrhosq*chisqrhosq*chisq
        r=Sqrt(r2)
        runit=DisplacementVector/r
        !Temporarily r6 = sigmasq / r2
        r6=intertpLJP(tp1,tp2).atmab(atm1,atm2).sigmasq/r2
        !Now r6 is what it is meant to be
        r6=r6*r6*r6
        r12=r6*r6
        !Temporarily runit = Van de Waals force_12 + Coulomb force_12 - 1D_Gz0_12
        runit=intertpLJP(tp1,tp2).atmab(atm1,atm2).epsilonm4*6d0/r*(2d0*r12-r6)*runit&
            +MoleculeInput(tp1).charge(atm1)*MoleculeInput(tp2).charge(atm2)&
            *((Nonbonded_chim2dsqrtpi*Exp(-chisq*r2)/r+Erfc(chi*r)/r2)*runit&
            +pim4dV*ReciprocalForceSummation1D(DisplacementVector))
        !Temporarily runit = f_12 = Van de Waals force_12 + Coulomb force_12
        runit(2:3)=runit(2:3)-MoleculeInput(tp1).charge(atm1)*MoleculeInput(tp2).charge(atm2)&
            *(Nonbonded_chisqm2dLz*(1-Exp(-chisqrhosq))/chisqrhosq*rhounit)
        NonbondedVirial1D(1)=DisplacementVector(1)*runit(1)
        !Temporarily r = Virial^12_inplane
        r=dot_product(DisplacementVector(2:3),runit(2:3))
        !Temporarily r12 = rho_1
        r12=Norm2(co1(2:3))
        !Temporarily r6 = rho_2
        r6=Norm2(co2(2:3))
        !Temporarily runit = f_12(2:3) * ( rho_2 - rho_1 )
        runit(2:3)=runit(2:3)*(r6-r12)
        NonbondedVirial1D(2)=r-dot_product(runit(2:3),co1(2:3)/r12)
        NonbondedVirial1D(3)=r-dot_product(runit(2:3),co2(2:3)/r6)
    end function NonbondedVirial1D
    !For 2D
    function NonbondedVirial2D(tp1,tp2,atm1,atm2,co1,co2)
        real*8,dimension(2)::NonbondedVirial2D
        integer,intent(in)::tp1,tp2,atm1,atm2
        real*8,dimension(3),intent(in)::co1,co2
        integer::k
        real*8::r,r2,r6,r12
        real*8,dimension(2)::runit
        !From 1 point to 2
        DisplacementVector=co2-co1
        !Periodic boundary condition
        do k=1,2
            if(DisplacementVector(k)>HalfLength(k)) then
                DisplacementVector(k)=DisplacementVector(k)-BoxSize(k)
            else if(DisplacementVector(k)<-HalfLength(k)) then
                DisplacementVector(k)=DisplacementVector(k)+BoxSize(k)
            end if
        end do
        r2=dot_product(DisplacementVector,DisplacementVector)
        r=Sqrt(r2)
        r6=intertpLJP(tp1,tp2).atmab(atm1,atm2).sigmasq/r2
        r6=r6*r6*r6
        r12=r6*r6
        runit=DisplacementVector(1:2)/r
        NonbondedVirial2D=DisplacementVector(1:2)*(&
            !Van de Waals force
            intertpLJP(tp1,tp2).atmab(atm1,atm2).epsilonm4*6d0/r*(2d0*r12-r6)*runit&
            !Coulomb force
            +MoleculeInput(tp1).charge(atm1)*MoleculeInput(tp2).charge(atm2)&
            *((Nonbonded_chim2dsqrtpi*Exp(-chisq*r2)/r+Erfc(chi*r)/r2)*runit&
            +pim4dV*ReciprocalForceSummation(DisplacementVector)))
    end function NonbondedVirial2D
!End need update
!--------------------- End ---------------------

end module MCAnalyzer