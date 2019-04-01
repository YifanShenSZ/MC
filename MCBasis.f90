!Monte Carlo basic routines
!Input and data storage
!
!Implementation detail:
!This program does not guarantee that every atom is inside the box,
!only the centres of mass of all molecules are kept inside the box
!Unit 99 is used in python "with" fasion, i.e., should be closed immediately after writing
module MCBasis
    use General
    implicit none

!Parameter
    !Accuracy control
    real*8::reltol=1d-4,abstol=1d-37
    integer::MaxIteration=100!I've tested and never seen 10 exceeded under reltol=1d-6
    !Wall parameters, default is carbon: C-C bond length 1.44A
    !    layer distance 3.35A, atom density 0.11382A^-3
    !    LJ epsilon 28K, sigma 3.4A
    character*2::WallElement='C'
    real*8::WallBondLength=2.721205632846602d0,&
        LayerDistance=6.330582548636193d0,WallAtomDensity=0.016866383610866523d0,&
        epsilonwall=8.87912d-5,sigmawall=6.42507d0

!Derived type
    !See examples of my special .xyz and .react for the meaning of each variable
    !Store the details of an input molecule, all same type molecules share these details
    type MoleculeDetailsShared
        character*2,allocatable,dimension(:)::ElementSymbol
        !How many atoms have non-zero epsilon/charge; non-zero epsilon but zero charge (vice versa); both non-zero epsilon and charge
        integer::NAtoms,NNon0Epsilon,NNon0Charge,NEpsilonOnly,NChargeOnly,NBoth
        !The indices of the atoms with non-zero epsilon/charge; non-zero epsilon but zero charge (vice versa); both non-zero epsilon and charge
        integer,allocatable,dimension(:)::Non0Epsilon,Non0Charge,EpsilonOnly,ChargeOnly,Both
        real*8::mass
        real*8,allocatable,dimension(:)::epsilon,sigma,charge
        real*8,allocatable,dimension(:,:)::RefConfig!Short for REFerence CONFIGuration
    end type MoleculeDetailsShared
    !Store the details of an input reaction
    type ReactionDetailsShared
        integer::NReactant,NProduct
        integer,allocatable,dimension(:)::ReactantSequence,ProductSequence
        real*8::deltaAig
    end type ReactionDetailsShared

    !Store single-body information: structure, single-body interaction
    !Example: type(moleculetype),allocatable,dimension(:)::mtp
    !         mtp(i).mol(j).Coordinate(k,l) is the k-th dimension displacement
    !         of the l-th atom in the j-th type i molecule
    !         mtp(i).mol(j).Total is the total single-body interaction energy
    type MoleculeDetails
        !ESSBDependent2D is the single-body shape dependent term in 2D Ewald summation
        !WallLJ is the outer field interaction arising from the LJ field of the wall
        !In this program, 1D only has WallLJ, 3D has no single-body interaction
        real*8::Total=0d0,ESSBDependent2D=0d0,WallLJ=0d0
        real*8,allocatable,dimension(:,:)::Coordinate
    end type MoleculeDetails
    type MoleculeType
        type(MoleculeDetails),allocatable,dimension(:)::mol!short for MOLecule
    end type MoleculeType
    !Store many-body information: pairwise interaction
    !Example: type(InterTypeInteraction),allocatable,dimension(:,:)::intertp
    !         intertp(i,j).molab(k,l).Total is the total pairwise interaction energy
    !         between k-th type i molecule and l-th type j molecule
    !         To avoid double counting, i>=j, k>l when i=j
    type InterMolecularInteraction
        !ESMB is the pairwise term in Ewald summation
        !LJ is the intermolecular Lennard-Jones interaction
        real*8::Total=0d0,ESMB=0d0,LJ=0d0
    end type InterMolecularInteraction
    type InterTypeInteraction
        type(InterMolecularInteraction),allocatable,dimension(:,:)::molab!short for MOLecule A and molecule B
    end type InterTypeInteraction

!Input variable
    !See input example for the meaning of each variable
    !Main input
    character*32::jobtype
    integer::TotalSteps,PeriodicDims
    real*8,dimension(3)::BoxSize
    integer::MoleculeKinds
    character*32,allocatable,dimension(:)::MoleculeFiles
    integer,allocatable,dimension(:)::MaxPossibleN,NExist
    character*32::ensemble
    real*8::temperature,pressure
    real*8,allocatable,dimension(:)::ChemicalPotential
    integer::ReactionKinds
    character*32,allocatable,dimension(:)::ReactionFiles
    !Molecule input
    type(MoleculeDetailsShared),allocatable,dimension(:)::MoleculeInput
    !Reaction input
    type(ReactionDetailsShared),allocatable,dimension(:)::ReactionInput

!Global variable
    real*8,allocatable,dimension(:)::oneMoleculeKinds
    !Public work space
        !When generating a unit quaternion, 2 length 2 vectors are required
        real*8,dimension(2)::QuaternionSupporter1,QuaternionSupporter2
        real*8,dimension(3)::DisplacementVector
        real*8,dimension(4)::RotationQuaternion
    !Supports Trial step
        !TrialCount(i,j) counts how many times translational trial is performed for the i-th direction of type j molecule
        !TrialSuccessCount counts only how many times this translational trial succeed
        integer,allocatable,dimension(:,:)::TrialCount,TrialSuccessCount
        !A translational trial for the i-th direction of type j molecule is within [-TrialStep(i,j),TrialStep(i,j)]
        real*8,allocatable,dimension(:,:)::TrialStep
    !Store structure and interaction, the existing molecules first. See "Derived type" for details
        type(moleculetype),allocatable,dimension(:)::mtp,mtpold!short for Molecule TyPe
        type(InterTypeInteraction),allocatable,dimension(:,:)::intertp,intertpold!short for INTER molecule TyPe interaction
    !Molecule and atom number counter
        integer::NTotal!The total number of existing molecules
        !The number of existing molecules before the trial move, the number of existing molecules at last check point
        integer,allocatable,dimension(:)::NExistOld,NExistLastTime
    !Size correlated variables
        !Radius and square of radius (for 1D only), volume of the simulation box
        real*8::radius,radiussq,volume
        !Half of the length of the simulation box, max translational trial move step
        real*8,dimension(3)::HalfLength,MaxTranStep
    !Translational partition function for each kind of molecule
        real*8,allocatable,dimension(:)::TranslationalPartitionFunction
    !Cumulative distribution function for choosing a type of molecule to run insertion trial
        real*8,allocatable,dimension(:)::ExChoose
    !Whether output the warning that MaxPossibleN reached
        logical,allocatable,dimension(:)::MaxPossibleNReachedWarning
    !Threshold to reject a trial without rolling dice, current potential energy per molecule,
    !the lowest potential energy per molecule, the lowest potential energy per molecule at last check point
        real*8::ImpossibleDel,AverageEnergy,AverageGroundEnergy,AverageGroundEnergyOld

contains
!The initializer for MCBasis module
subroutine InitializeMCBasis(ContinueAnalyzation)
    logical,intent(out)::ContinueAnalyzation
    integer::i,j,index,k
    real*8::temp
    call BetterRandomSeed()
    !A randomly generated structure usually has a high energy
    !It needs more (thus faster) trials rather than accuracy
    if(jobtype=='NewStructure') then
        reltol=1d-2
    !Generally, the magnitude of intermolecular interaction is tens kJ/mol per molecule,
    !This package usually deals with hundreds of molecules, and my ideal accuracy
    !for my result is the absolute error of total energy <= 0.1 kJ/mol
    else if(jobtype=='Analyze') then
        reltol=1d-6
    end if
    !Variables needed to be initialized before reading check point files
    allocate(TrialStep(3,MoleculeKinds))
    allocate(mtp(MoleculeKinds))
    allocate(mtpold(MoleculeKinds))
    allocate(intertp(MoleculeKinds,MoleculeKinds))
    allocate(intertpold(MoleculeKinds,MoleculeKinds))
    do i=1,MoleculeKinds
        allocate(mtp(i).mol(MaxPossibleN(i)))
        allocate(mtpold(i).mol(MaxPossibleN(i)))
        do j=1,MaxPossibleN(i)
            allocate(mtp(i).mol(j).Coordinate(3,MoleculeInput(i).NAtoms))
            allocate(mtpold(i).mol(j).Coordinate(3,MoleculeInput(i).NAtoms))
        end do
        do j=i,MoleculeKinds
            allocate(intertp(j,i).molab(MaxPossibleN(j),MaxPossibleN(i)))
            allocate(intertpold(j,i).molab(MaxPossibleN(j),MaxPossibleN(i)))
        end do
    end do
    call ReadCheckPoint(ContinueAnalyzation)
    !Initialize global variables by definition order
    allocate(oneMoleculeKinds(MoleculeKinds))
    oneMoleculeKinds=1d0
    !DisplacementVector, RotationQuaternion do not need initialization
    allocate(TrialCount(3,MoleculeKinds))
        TrialCount=0
    allocate(TrialSuccessCount(3,MoleculeKinds))
        TrialSuccessCount=0
    !TrialStep has been initialized in subroutine ReadCheckPoint
    !mtp.mol.WallJ, mtpold, intertp, intertpold will be initialized in main program for NewStructure,
    !otherwise have been initialized in subroutine ReadCheckPoint
    NTotal=0
    do i=1,MoleculeKinds
        NTotal=NTotal+NExist(i)
    end do
    allocate(NExistOld(MoleculeKinds))
    NExistOld=NExist
    allocate(NExistLastTime(MoleculeKinds))
    NExistLastTime=NExistOld
    !radius, radiussq, volume, HalfLength have been initialized in subroutine ReadCheckPoint
    MaxTranStep=HalfLength
    allocate(TranslationalPartitionFunction(MoleculeKinds))
    forall(i=1:MoleculeKinds)
        TranslationalPartitionFunction(i)=(MoleculeInput(i).mass*temperature/pim2)**1.5d0*volume
    end forall
    !Here we let the probability of choosing a type of molecule to run insertion trial
    !be propotional to the input existance number of each type of molecules
    allocate(ExChoose(MoleculeKinds))
    temp=0d0
    do i=1,MoleculeKinds
        ExChoose(i)=temp+NExist(i)
        temp=ExChoose(i)
    end do
    ExChoose=ExChoose/NTotal
    !Default value is True, will be set to False once a warning is outputed
    allocate(MaxPossibleNReachedWarning(MoleculeKinds))
    MaxPossibleNReachedWarning=1
    !The minimum random number Fortran can produce is 1d0/2147483647d0
    !So if Del>ImpossibleDel, this trial move will be rejected for sure
    ImpossibleDel=temperature*21.487562596892646d0
    !AverageEnergy,AverageGroundEnergy,AverageGroundEnergyOld will be initialized in main program for NewStructure,
    !otherwise have been initialized in subroutine ReadCheckPoint
end subroutine InitializeMCBasis

!---------- Monte Carlo specific basic subroutines and function ----------
    !Generate the displacement vector of the centre of mass, and the rotational quaternion about the centre of mass
    !If don't want to insert uniformly, will be inserted at least sigmawall from the wall
    !For 1D
    subroutine InsertMolecule1D(tran,rt,uniformly)
        real*8,dimension(3),intent(inout)::tran
        real*8,dimension(4),intent(inout)::rt
        logical,intent(in)::uniformly
        !Uniformly sample a point inside the unit circle
        do while(1)
            call random_number(QuaternionSupporter1)
            QuaternionSupporter1=QuaternionSupporter1*2d0-1d0
            if(dot_product(QuaternionSupporter1,QuaternionSupporter1)<1d0) then
                exit
            end if
        end do
        !Convert the sampling from unit circle to actual tube
        if(uniformly) then
            tran(2:3)=QuaternionSupporter1*radius
        else
            tran(2:3)=QuaternionSupporter1*(radius-sigmawall)
        end if
        call random_number(tran(1))
        tran(1)=(tran(1)*2d0-1d0)*HalfLength(1)
        !Form the unit quaternion
            !Uniformly sample a point inside the unit circle
            do while(1)
                call random_number(QuaternionSupporter1)
                QuaternionSupporter1=QuaternionSupporter1*2d0-1d0
                if(dot_product(QuaternionSupporter1,QuaternionSupporter1)<1d0) then
                    exit
                end if
            end do
            !Uniformly sample another point inside the unit circle
            do while(1)
                call random_number(QuaternionSupporter2)
                QuaternionSupporter2=QuaternionSupporter2*2d0-1d0
                if(dot_product(QuaternionSupporter2,QuaternionSupporter2)<1d0) then
                    exit
                end if
            end do
            !Convert to the unit quaternion
            rt(1:2)=QuaternionSupporter1
            rt(3:4)=QuaternionSupporter2&
                *Sqrt((1d0-dot_product(QuaternionSupporter1,QuaternionSupporter1))&
                /dot_product(QuaternionSupporter2,QuaternionSupporter2))
    end subroutine InsertMolecule1D
    !For 2D
    subroutine InsertMolecule2D(tran,rt,uniformly)
        real*8,dimension(3),intent(inout)::tran
        real*8,dimension(4),intent(inout)::rt
        logical,intent(in)::uniformly
        call RANDOM_NUMBER(tran)
        if(uniformly) then
            tran=(tran*2d0-1d0)*HalfLength
        else
            tran(1:2)=(tran(1:2)*2d0-1d0)*HalfLength(1:2)
            tran(3)=tran(3)*(HalfLength(3)-sigmawall)
        end if
        !Form the unit quaternion
            !Uniformly sample a point inside the unit circle
            do while(1)
                call random_number(QuaternionSupporter1)
                QuaternionSupporter1=QuaternionSupporter1*2d0-1d0
                if(dot_product(QuaternionSupporter1,QuaternionSupporter1)<1d0) then
                    exit
                end if
            end do
            !Uniformly sample another point inside the unit circle
            do while(1)
                call random_number(QuaternionSupporter2)
                QuaternionSupporter2=QuaternionSupporter2*2d0-1d0
                if(dot_product(QuaternionSupporter2,QuaternionSupporter2)<1d0) then
                    exit
                end if
            end do
            !Convert to the unit quaternion
            rt(1:2)=QuaternionSupporter1
            rt(3:4)=QuaternionSupporter2&
                *Sqrt((1d0-dot_product(QuaternionSupporter1,QuaternionSupporter1))&
                /dot_product(QuaternionSupporter2,QuaternionSupporter2))
    end subroutine InsertMolecule2D
    !For 3D
    subroutine InsertMolecule3D(tran,rt)
        real*8,dimension(3),intent(inout)::tran
        real*8,dimension(4),intent(inout)::rt
        call RANDOM_NUMBER(tran)
        tran=(tran*2d0-1d0)*HalfLength
        !Form the unit quaternion
            !Uniformly sample a point inside the unit circle
            do while(1)
                call random_number(QuaternionSupporter1)
                QuaternionSupporter1=QuaternionSupporter1*2d0-1d0
                if(dot_product(QuaternionSupporter1,QuaternionSupporter1)<1d0) then
                    exit
                end if
            end do
            !Uniformly sample another point inside the unit circle
            do while(1)
                call random_number(QuaternionSupporter2)
                QuaternionSupporter2=QuaternionSupporter2*2d0-1d0
                if(dot_product(QuaternionSupporter2,QuaternionSupporter2)<1d0) then
                    exit
                end if
            end do
            !Convert to the unit quaternion
            rt(1:2)=QuaternionSupporter1
            rt(3:4)=QuaternionSupporter2&
                *Sqrt((1d0-dot_product(QuaternionSupporter1,QuaternionSupporter1))&
                /dot_product(QuaternionSupporter2,QuaternionSupporter2))
    end subroutine InsertMolecule3D
    
    !Struct the present coordinate of an atom from translating and rotating the reference coordinate
    function Struct(tran,rt,co0,n)
        integer,intent(in)::n
        real*8,dimension(3,n)::struct
        real*8,dimension(3),intent(in)::tran
        real*8,dimension(4),intent(in)::rt
        real*8,dimension(3,n),intent(in)::co0
        integer::i
        do i=1,n
            struct(:,i)=tran+rotate(rt,co0(:,i))
        end do
    end function Struct
!---------------------------------- End ----------------------------------

!-------------------- Non-analyzation input / output ---------------------
    !Read all input files: input, .xyz, .react (optional)
    subroutine ReadInput(AdvancedInput)
        logical,intent(out)::AdvancedInput
        logical::NeedRTrial
        integer::i,j
        integer,allocatable,dimension(:)::IndicesTemp,IndicesTemp1,IndicesTemp2
        !Set default value
            AdvancedInput=.false.
            NeedRTrial=.false.
        !Read main input, convert to atomic unit,
        !initialize variables required in other inputs, write some job comment
        open(unit=99,file='input',status='old')
            read(99,*)
            read(99,*)jobtype
                write(*,*)'Job type: '//jobtype
            read(99,*)
            read(99,*)TotalSteps
            read(99,*)
            read(99,*)PeriodicDims
            read(99,*)
            read(99,*)
            read(99,'(F20.15,F20.15,F20.15)')BoxSize
                !Convert to atomic unit
                BoxSize=BoxSize*AInAU
            read(99,*)
            read(99,*)MoleculeKinds
            read(99,*)
            allocate(MoleculeFiles(MoleculeKinds))
            do i=1,MoleculeKinds
                read(99,*)MoleculeFiles(i)
            end do
            read(99,*)
            allocate(MaxPossibleN(MoleculeKinds))
            do i=1,MoleculeKinds
                read(99,*)MaxPossibleN(i)
            end do
            read(99,*)
            allocate(NExist(MoleculeKinds))
            do i=1,MoleculeKinds
                read(99,*)NExist(i)
            end do
            read(99,*)
            read(99,*)ensemble
                if(jobtype=='NewStructure') then
                    write(*,'(1x,A55)')'Optimize a randomly generated structure in NVT ensemble'
                else
                    if(ensemble=='RVT'.or.ensemble=='RPT') then
                        NeedRTrial=.true.!Reaction requested, read reaction input
                    end if
                    write(*,*)'Ensemble: '//ensemble
                end if
            read(99,*)
            read(99,*)temperature
                !Convert to atomic unit
                temperature=temperature*KInAU
            read(99,*)
            read(99,*)pressure
                !Convert to atomic unit
                pressure=pressure*barInAU
            read(99,*)
            allocate(ChemicalPotential(MoleculeKinds))
            do i=1,MoleculeKinds
                read(99,*)ChemicalPotential(i)
            end do
                !Convert to atomic unit
                ChemicalPotential=ChemicalPotential*kJmolInAU
            read(99,*)
            read(99,*)AdvancedInput
            if(AdvancedInput) write(*,*)'Advanced input requested, parameters are set to user specification'
            if(NeedRTrial) then
                !If reaction is requested, read reaction input
                read(99,*)
                read(99,*)ReactionKinds
                read(99,*)
                read(99,*)
                allocate(ReactionFiles(ReactionKinds))
                do i=1,ReactionKinds
                    read(99,*)ReactionFiles(i)
                end do
            end if
        close(99)
        !Read molecule input
        allocate(MoleculeInput(MoleculeKinds))
        do i=1,MoleculeKinds
            open(unit=99,file=MoleculeFiles(i),status='old')
                read(99,*)
                read(99,*)
                read(99,*)
                read(99,*)MoleculeInput(i).NAtoms
                read(99,*)
                allocate(MoleculeInput(i).ElementSymbol(MoleculeInput(i).NAtoms))
                !The j-th column is the displacement vector of the j-th atom
                allocate(MoleculeInput(i).RefConfig(3,MoleculeInput(i).NAtoms))
                do j=1,MoleculeInput(i).NAtoms
                    read(99,'(A2,F20.15,F20.15,F20.15)')MoleculeInput(i).ElementSymbol(j),&
                        MoleculeInput(i).RefConfig(1,j),MoleculeInput(i).RefConfig(2,j),MoleculeInput(i).RefConfig(3,j)
                end do
                    !Convert to atomic unit
                    MoleculeInput(i).RefConfig=MoleculeInput(i).RefConfig*AInAU
                read(99,*)
                read(99,*)MoleculeInput(i).mass
                    !Convert to atomic unit
                    MoleculeInput(i).mass=MoleculeInput(i).mass*AMUInAU
                read(99,*)
                !Work space for temporarily saving the indices of the atoms with non-zero epsilon/charge
                allocate(IndicesTemp(MoleculeInput(i).NAtoms))
                !Read epsilon and count how many atoms have non-zero epsilon
                allocate(MoleculeInput(i).epsilon(MoleculeInput(i).NAtoms))
                MoleculeInput(i).NNon0Epsilon=0
                do j=1,MoleculeInput(i).NAtoms
                    read(99,*)MoleculeInput(i).epsilon(j)
                    if(MoleculeInput(i).epsilon(j)/=0d0) then
                        MoleculeInput(i).NNon0Epsilon=MoleculeInput(i).NNon0Epsilon+1
                        !Temporarily save the indices of the atoms with non-zero epsilon
                        IndicesTemp(MoleculeInput(i).NNon0Epsilon)=j
                    end if
                end do
                    !Convert to atomic unit
                    MoleculeInput(i).epsilon=MoleculeInput(i).epsilon*kJmolInAU
                !Permanently save the indices of the atoms with non-zero epsilon
                allocate(MoleculeInput(i).Non0Epsilon(MoleculeInput(i).NNon0Epsilon))
                MoleculeInput(i).Non0Epsilon=IndicesTemp(1:MoleculeInput(i).NNon0Epsilon)
                read(99,*)
                !Read sigma
                allocate(MoleculeInput(i).sigma(MoleculeInput(i).NAtoms))
                do j=1,MoleculeInput(i).NAtoms
                    read(99,*)MoleculeInput(i).sigma(j)
                end do
                    !Convert to atomic unit
                    MoleculeInput(i).sigma=MoleculeInput(i).sigma*AInAU
                read(99,*)
                !Read charge and count how many atoms have non-zero charge
                allocate(MoleculeInput(i).charge(MoleculeInput(i).NAtoms))
                MoleculeInput(i).NNon0Charge=0
                do j=1,MoleculeInput(i).NAtoms
                    read(99,*)MoleculeInput(i).charge(j)
                    if(MoleculeInput(i).charge(j)/=0d0) then
                        MoleculeInput(i).NNon0Charge=MoleculeInput(i).NNon0Charge+1
                        !Temporarily save the indices of the atoms with non-zero charge
                        IndicesTemp(MoleculeInput(i).NNon0Charge)=j
                    end if
                end do
                !Permanently save the indices of the atoms with non-zero epsilon
                allocate(MoleculeInput(i).Non0Charge(MoleculeInput(i).NNon0Charge))
                MoleculeInput(i).Non0Charge=IndicesTemp(1:MoleculeInput(i).NNon0Charge)
                !Count how many atoms have non-zero epsilon but zero charge (vice versa); both non-zero epsilon and charge
                allocate(IndicesTemp1(MoleculeInput(i).NAtoms))
                allocate(IndicesTemp2(MoleculeInput(i).NAtoms))
                MoleculeInput(i).NEpsilonOnly=0
                MoleculeInput(i).NChargeOnly=0
                MoleculeInput(i).NBoth=0
                do j=1,MoleculeInput(i).NAtoms
                    if(MoleculeInput(i).epsilon(j)/=0.and.MoleculeInput(i).charge(j)==0) then
                        MoleculeInput(i).NEpsilonOnly=MoleculeInput(i).NEpsilonOnly+1
                        !Temporarily save the indices of the atoms with non-zero epsilon but zero charge
                        IndicesTemp(MoleculeInput(i).NEpsilonOnly)=j
                    else if(MoleculeInput(i).epsilon(j)==0.and.MoleculeInput(i).charge(j)/=0) then
                        MoleculeInput(i).NChargeOnly=MoleculeInput(i).NChargeOnly+1
                        !Temporarily save the indices of the atoms with non-zero charge but zero epsilon
                        IndicesTemp1(MoleculeInput(i).NChargeOnly)=j
                    else if(MoleculeInput(i).epsilon(j)/=0.and.MoleculeInput(i).charge(j)/=0) then
                        MoleculeInput(i).NBoth=MoleculeInput(i).NBoth+1
                        !Temporarily save the indices of the atoms with both non-zero epsilon and charge
                        IndicesTemp2(MoleculeInput(i).NBoth)=j
                    end if
                end do
                !Permanently save the indices of the atoms with non-zero epsilon but zero charge (vice versa); both non-zero epsilon and charge
                allocate(MoleculeInput(i).EpsilonOnly(MoleculeInput(i).NEpsilonOnly))
                MoleculeInput(i).EpsilonOnly=IndicesTemp(1:MoleculeInput(i).NEpsilonOnly)
                allocate(MoleculeInput(i).ChargeOnly(MoleculeInput(i).NChargeOnly))
                MoleculeInput(i).ChargeOnly=IndicesTemp1(1:MoleculeInput(i).NChargeOnly)
                allocate(MoleculeInput(i).Both(MoleculeInput(i).NBoth))
                MoleculeInput(i).Both=IndicesTemp2(1:MoleculeInput(i).NBoth)
                !Clean work space for next loop
                deallocate(IndicesTemp)
                deallocate(IndicesTemp1)
                deallocate(IndicesTemp2)
            close(99)
        end do
        !If reaction is requested, read reaction input
        if(NeedRTrial) then
            allocate(ReactionInput(ReactionKinds))
            do i=1,ReactionKinds
                open(unit=99,file=ReactionFiles(i),status='old')
                    read(99,*)
                    read(99,*)ReactionInput(i).NReactant
                    read(99,*)
                    allocate(ReactionInput(i).ReactantSequence(ReactionInput(i).NReactant))
                    do j=1,ReactionInput(i).NReactant
                        read(99,*)ReactionInput(i).ReactantSequence(j)
                    end do
                    read(99,*)
                    read(99,*)ReactionInput(i).NProduct
                    read(99,*)
                    allocate(ReactionInput(i).ProductSequence(ReactionInput(i).NProduct))
                    do j=1,ReactionInput(i).NProduct
                        read(99,*)ReactionInput(i).ProductSequence(j)
                    end do
                    read(99,*)
                    read(99,*)ReactionInput(i).deltaAig
                close(99)
            end do
        end if
    end subroutine ReadInput
    
    !If not a NewStructure job, read the check point files
    !Otherwise, initialize same variables read or initialized during reading the check point files
    subroutine ReadCheckPoint(ContinueAnalyzation)
        logical,intent(out)::ContinueAnalyzation
        character*2::char2temp
        character*17::char17temp
        integer::i,j,ii,jj
        real*8::temp
        !If not a NewStructure job, read the check point files
        if(jobtype/='NewStructure') then
            !Read main check point file, convert to atomic unit,
            !initialize variables required in other inputs, write some job comment
            open(unit=99,file='MonteCarlo.CheckPoint',status='old')
                read(99,*)char17temp
                    if(char17temp=='Enjoy your result') then
                        ContinueAnalyzation=.true.
                    end if
                read(99,*)
                !Replace input box size with what in check point file
                read(99,'(F20.15,F20.15,F20.15)')BoxSize
                    !Convert to atomic unit
                    BoxSize=BoxSize*AInAU
                read(99,*)
                do i=1,MoleculeKinds
                    read(99,*)TrialStep(:,i)
                end do
                    !Convert to atomic unit
                    TrialStep=TrialStep*AInAU
                read(99,*)
                read(99,*)NExist
                read(99,*)
                read(99,*)AverageEnergy
                    !Convert to atomic unit
                    AverageEnergy=AverageEnergy*kJmolInAU
                read(99,*)
                read(99,*)AverageGroundEnergyOld
                    !Convert to atomic unit
                    AverageGroundEnergyOld=AverageGroundEnergyOld*kJmolInAU
                    AverageGroundEnergy=AverageGroundEnergyOld
            close(99)
            !Read single-body interaction
            open(unit=99,file='SingleBody.CheckPoint',status='old')
                do i=1,MoleculeKinds
                    do j=1,NExist(i)
                        read(99,*)mtp(i).mol(j).ESSBDependent2D
                        read(99,*)mtp(i).mol(j).WallLJ
                        !Sum up as total single-body interaction
                        mtp(i).mol(j).Total=mtp(i).mol(j).ESSBDependent2D+mtp(i).mol(j).WallLJ
                    end do
                end do
            close(99)
            !Read many-body interaction
            open(unit=99,file='ManyBody.CheckPoint',status='old')
                do j=1,MoleculeKinds
                    do i=1,NExist(j)
                        do ii=i+1,NExist(j)
                            read(99,*)intertp(j,j).molab(ii,i).ESMB
                            read(99,*)intertp(j,j).molab(ii,i).LJ
                            !Sum up as total many-body interaction
                            intertp(j,j).molab(ii,i).Total=intertp(j,j).molab(ii,i).ESMB+intertp(j,j).molab(ii,i).LJ
                        end do
                        do jj=j+1,MoleculeKinds
                            do ii=1,NExist(jj)
                                read(99,*)intertp(jj,j).molab(ii,i).ESMB
                                read(99,*)intertp(jj,j).molab(ii,i).LJ
                                !Sum up as total many-body interaction
                                intertp(jj,j).molab(ii,i).Total=intertp(jj,j).molab(ii,i).ESMB+intertp(jj,j).molab(ii,i).LJ
                            end do
                        end do
                    end do
                end do
            close(99)
            !Read system structure
            open(unit=99,file='Structure.xyz',status='old')
                read(99,*)
                read(99,*)
                do i=1,MoleculeKinds
                    do j=1,NExist(i)
                        do ii=1,MoleculeInput(i).NAtoms
                            read(99,'(A2,F20.15,F20.15,F20.15)')char2temp,&
                                mtp(i).mol(j).Coordinate(1,ii),mtp(i).mol(j).Coordinate(2,ii),mtp(i).mol(j).Coordinate(3,ii)
                        end do
                        !Convert to atomic unit
                        mtp(i).mol(j).Coordinate=mtp(i).mol(j).Coordinate*AInAU
                    end do
                end do
            close(99)
            mtpold=mtp
            intertpold=intertp
            !Same things are done in else
            if(PeriodicDims==1) then
                radius=BoxSize(2)
                radiussq=radius*radius
                HalfLength(1)=BoxSize(1)/2d0
                !Used as the maxmium translational trial move step
                HalfLength(2)=radius/Sqrt(2d0)
                HalfLength(3)=HalfLength(2)
                !Volume of the box
                volume=pi*radiussq*BoxSize(1)
            else if(PeriodicDims==2) then
                HalfLength=BoxSize/2d0
                !Volume of the box
                volume=BoxSize(1)*BoxSize(2)*BoxSize(3)
            else
                HalfLength=BoxSize/2d0
                !Volume of the box
                volume=BoxSize(1)*BoxSize(2)*BoxSize(3)
            end if
        !Otherwise, initialize same variables read or initialized during reading the check point files
        else
            !Rescale z size to fit the period of the wall
            if(PeriodicDims==1) then
                !The height of the hexagon
                temp=Sqrt(3d0)*WallBondLength
                !The box should contain integer number of hexagons in z direction
                BoxSize(1)=ceiling(BoxSize(1)/temp)*temp
                write(*,'(A51,F20.15,A1)')'z size is rescaled to fit the period of the wall: ',BoxSize(1)/AInAU,'A'
            !Rescale xy size to fit the period of the wall
            else if(PeriodicDims==2) then
                !Hexagonal symmetry requires a unit cell with 3 atoms align as ethylene
                !1 at each corner (H place) and 2 at center (C place)
                !The width of the unit cell (perpendicular to C=C)
                temp=Sqrt(3d0)*WallBondLength
                !The box should contain integer number of hexagons in x direction
                BoxSize(1)=ceiling(BoxSize(1)/temp)*temp
                !The length of the unit cell (parallel to C=C)
                temp=3d0*WallBondLength
                !The box should contain integer number of hexagons in y direction
                BoxSize(2)=ceiling(BoxSize(2)/temp)*temp
                write(*,'(A51,F20.15,A1)')'x size is rescaled to fit the period of the wall: ',BoxSize(1)/AInAU,'A'
                write(*,'(A51,F20.15,A1)')'y size is rescaled to fit the period of the wall: ',BoxSize(2)/AInAU,'A'
            end if
            !Default initial value of Trial steps
            TrialStep=1d0
            !Cannot initialize AverageEnergy, AverageGroundEnergy, intertp, intertpold now since they require functions in later modules
            !NewStrucutre does not need AverageGroundEnergyOld
            if(PeriodicDims==1) then
                !Before inserting molecules, radius and HalfLength must be ready
                radius=BoxSize(2)
                radiussq=radius*radius
                HalfLength(1)=BoxSize(1)/2d0
                !Used as the maxmium translational trial move step
                HalfLength(2)=radius/Sqrt(2d0)
                HalfLength(3)=HalfLength(2)
                !Volume of the box
                volume=pi*radiussq*BoxSize(1)
                !Insert molecules
                do i=1,MoleculeKinds
                    do j=1,NExist(i)
                        call InsertMolecule1D(DisplacementVector,RotationQuaternion,.false.)
                        mtp(i).mol(j).Coordinate=struct(DisplacementVector,RotationQuaternion,&
                            MoleculeInput(i).RefConfig,MoleculeInput(i).NAtoms)
                    end do
                end do
            else if(PeriodicDims==2) then
                !Before inserting molecules, HalfLength must be ready
                HalfLength=BoxSize/2d0
                !Volume of the box
                volume=BoxSize(1)*BoxSize(2)*BoxSize(3)
                !Insert molecules
                do i=1,MoleculeKinds
                    do j=1,NExist(i)
                        call InsertMolecule2D(DisplacementVector,RotationQuaternion,.false.)
                        mtp(i).mol(j).Coordinate=struct(DisplacementVector,RotationQuaternion,&
                            MoleculeInput(i).RefConfig,MoleculeInput(i).NAtoms)
                    end do
                end do
            else
                !Before inserting molecules, HalfLength must be ready
                HalfLength=BoxSize/2d0
                !Volume of the box
                volume=BoxSize(1)*BoxSize(2)*BoxSize(3)
                !Insert molecules
                do i=1,MoleculeKinds
                    do j=1,NExist(i)
                        call InsertMolecule3D(DisplacementVector,RotationQuaternion)
                        mtp(i).mol(j).Coordinate=struct(DisplacementVector,RotationQuaternion,&
                            MoleculeInput(i).RefConfig,MoleculeInput(i).NAtoms)
                    end do
                end do
            end if
            !Cannot initialize mtp.mol.SingleBodyInteraction, mtpold since they requrie functions in later modules
        end if
    end subroutine ReadCheckPoint

    !Check point files writer
    subroutine WriteCheckPoint()
        logical::flag
        integer::i,j,k,ii,jj
        !Write main check point
        open(unit=99,file='MonteCarlo.CheckPoint',status='replace')
            if(jobtype=='NewStructure') then
                write(99,'(A21)')'Structure initialized'
            else if(jobtype=='Optimize') then
                if(jobtype=='NVT'.or.jobtype=='NPT') then
                    if(AverageGroundEnergy<AverageGroundEnergyOld) then
                        AverageGroundEnergyOld=AverageGroundEnergy
                        write(99,'(A25)')'Need further optimization'
                    else
                        write(99,'(A21)')'Optimization complete'
                    end if
                else
                    if(AverageGroundEnergy<AverageGroundEnergyOld) then
                        AverageGroundEnergyOld=AverageGroundEnergy
                        write(99,'(A25)')'Need further optimization'
                        flag=0
                    else
                        flag=1
                    end if
                    if(flag) then
                        do i=1,MoleculeKinds
                            if(NExist(i)/=NExistLastTime(i)) then
                                write(99,'(A25)')'Need further optimization'
                                NExistLastTime(i:MoleculeKinds)=NExist(i:MoleculeKinds)
                                flag=0
                                exit
                            end if
                        end do
                    end if
                    if(flag) then
                        write(99,'(A21)')'Optimization complete'
                    end if
                end if
            else
                write(99,'(A17)')'Enjoy your result'
            end if
            write(99,'(A17)')'Current box size:'
            write(99,'(F20.15,F20.15,F20.15)')BoxSize/AInAU
            write(99,'(A33)')'Current translational trail step:'
            do i=1,MoleculeKinds
                write(99,*)TrialStep(:,i)/AInAU
            end do
            write(99,'(A40)')'Current number of each kind of molecule:'
            write(99,*)NExist
            write(99,'(A51)')'Current average configurational energy: (In kJ/mol)'
            write(99,*)AverageEnergy/kJmolInAU
            write(99,'(A66)')'The lowest average configurational energy experienced: (In kJ/mol)'
            write(99,*)AverageGroundEnergy/kJmolInAU
        close(99)
        !Write the current single-body interaction
        open(unit=99,file='SingleBody.CheckPoint',status='replace')
            do i=1,MoleculeKinds
                do j=1,NExist(i)
                    write(99,*)mtp(i).mol(j).ESSBDependent2D
                    write(99,*)mtp(i).mol(j).WallLJ
                end do
            end do
        close(99)
        !Write the current many-body interaction
        open(unit=99,file='ManyBody.CheckPoint',status='replace')
            do j=1,MoleculeKinds
                do i=1,NExist(j)
                    do ii=i+1,NExist(j)
                        write(99,*)intertp(j,j).molab(ii,i).ESMB
                        write(99,*)intertp(j,j).molab(ii,i).LJ
                    end do
                    do jj=j+1,MoleculeKinds
                        do ii=1,NExist(jj)
                            write(99,*)intertp(jj,j).molab(ii,i).ESMB
                            write(99,*)intertp(jj,j).molab(ii,i).LJ
                        end do
                    end do
                end do
            end do
        close(99)
        !Write the current structure
        open(unit=99,file='Structure.xyz',status='replace')
            !xyz file requires the total number of atoms
            j=0
            do i=1,MoleculeKinds
                j=j+NExist(i)*MoleculeInput(i).NAtoms
            end do
            write(99,*)j
            write(99,*)
            do i=1,MoleculeKinds
                do j=1,NExist(i)
                    do k=1,MoleculeInput(i).NAtoms
                        write(99,'(A2,F20.15,F20.15,F20.15)')MoleculeInput(i).ElementSymbol(k),&
                            mtp(i).mol(j).Coordinate(1,k)/AInAU,mtp(i).mol(j).Coordinate(2,k)/AInAU,mtp(i).mol(j).Coordinate(3,k)/AInAU
                    end do
                end do
            end do
        close(99)
    end subroutine WriteCheckPoint
!---------------------------------- End ----------------------------------

end module MCBasis