!This is the module to account for nonbonded (intermolecular) energy
!Including Lennard-Jones potential and Coulomb potential
!Lennard-Jones potential is computed directly with cut-off = half length of the box
!Coulomb potential is computed by Ewald summation technique
!Reference: 1D: Brodka 2002 Chem. Phys. Lett.
!           2D: Yeh and Berkowitz 1999 JCP
!           3D: https://www.nist.gov/mml/csd/chemical-informatics-research-group/spce-water-reference-calculations-10%C3%A5-cutoff
!A molecule is specified by its kind number and coordinate
module Nonbonded
    use General
    use MCBasis
    implicit none

!Derived type
    type LennardJonesParameter
        real*8::epsilonm4,sigmasq
    end type LennardJonesParameter
    !Example: type(MoleculeTypeWallLJParameter),allocatable,dimension(:)::mtpWLJP
    !         mtpWLJP(i).atm(j).epsilonm4 is 4 times of the epsilon,
    !         mtpWLJP(i).atm(j).sigmasq is the sqaure of the sigma,
    !         between the j-th atom in type i molecule and the wall
    type MoleculeTypeWallLJParameter
        type(LennardJonesParameter),allocatable,dimension(:)::atm!short for AToM
    end type MoleculeTypeWallLJParameter
    !Example: type(InterTypeLJParameter),allocatable,dimension(:,:)::intertpLJP
    !         intertpLJP(i,j).atmab(k,l).epsilonm4 is 4 times of the epsilon,
    !         intertpLJP(i,j).atmab(k,l).sigmasq is the sqaure of the sigma,
    !         between the k-th atom in type i molecule and the l-th atom in type j molecule
    type InterTypeLJParameter
        type(LennardJonesParameter),allocatable,dimension(:,:)::atmab!short for AToM A and atom B
    end type InterTypeLJParameter

!Lennard-Jones potential variable
    !The coefficient except epsilon and sigma square in 10-4-3 steele potential, for 2D only
        real*8::SteeleCoefficient
    !Wall structure, for 1D only
        real*8,allocatable,dimension(:,:)::wall
    !Store Lennard-Jones parameters. See "Derived type" for details
        type(MoleculeTypeWallLJParameter),allocatable,dimension(:)::mtpWLJP!short for Molecule TyPe Wall LJ Parameter
        type(InterTypeLJParameter),allocatable,dimension(:,:)::intertpLJP!short for INTER molecule TyPe LJ Parameter

!Ewald summation variable
    !1D PBC is Ewald summed by somewhat a vacuum layer in x and y directions making them periodic 
    !2D PBC is Ewald summed by adding a vacuum layer in z direction making it periodic
    ! pim2dV = 2pi / volume, pim4dV = 4pi / volume (or volume_vacuum for 1D and 2D)
    real*8::volume_vacuum,pim2dV,pim4dV
    real*8,dimension(3)::wavenumber
    !chi is the dump factor in Ewald summation, chisq = chi**2, minus4chisq = -4chi**2
    real*8::chi,chisq,minus4chisq
    !Single-body position independent term in Ewald summation
    real*8,allocatable,dimension(:)::ESSBIndependent,ESSBIndependentOld
    !Force calculation
    real*8::Nonbonded_chim2dsqrtpi,Nonbonded_chisqm2dLz!For 1D only

contains
!The initializer for Nonbonded module
subroutine InitializeNonbonded()
    call InitializeLJ()
    call InitializeEwaldSummation()
end subroutine InitializeNonbonded

!------------------------ Van de Waals interaction: Lennard-Jones potential --------------------------
    !The initializer for Lennard-Jones potential
    subroutine InitializeLJ()
        integer::i,j,k,l
        SteeleCoefficient=pi/2d0*WallAtomDensity*LayerDistance
        if(PeriodicDims==1) then
            call GenerateWall1D()
        end if
        if(PeriodicDims<3) then
            allocate(mtpWLJP(MoleculeKinds))
            do i=1,MoleculeKinds
                allocate(mtpWLJP(i).atm(MoleculeInput(i).NAtoms))
                forall(j=1:MoleculeInput(i).NAtoms)
                    !Lorentz-Berthelot mixing rule
                    mtpWLJP(i).atm(j).epsilonm4=4d0*Sqrt(MoleculeInput(i).epsilon(j)*epsilonwall)
                    mtpWLJP(i).atm(j).sigmasq=(MoleculeInput(i).sigma(j)+sigmawall)/2d0
                    mtpWLJP(i).atm(j).sigmasq=mtpWLJP(i).atm(j).sigmasq*mtpWLJP(i).atm(j).sigmasq
                end forall
            end do
        end if
        allocate(intertpLJP(MoleculeKinds,MoleculeKinds))
        do i=1,MoleculeKinds
            do j=1,MoleculeKinds
                allocate(intertpLJP(i,j).atmab(MoleculeInput(i).NAtoms,MoleculeInput(j).NAtoms))
                forall(k=1:MoleculeInput(i).NAtoms,l=1:MoleculeInput(j).NAtoms)
                    !Lorentz-Berthelot mixing rule
                    intertpLJP(i,j).atmab(k,l).epsilonm4=4d0*Sqrt(MoleculeInput(i).epsilon(k)*MoleculeInput(j).epsilon(l))
                    intertpLJP(i,j).atmab(k,l).sigmasq=(MoleculeInput(i).sigma(k)+MoleculeInput(j).sigma(l))/2d0
                    intertpLJP(i,j).atmab(k,l).sigmasq=intertpLJP(i,j).atmab(k,l).sigmasq*intertpLJP(i,j).atmab(k,l).sigmasq
                end forall
            end do
        end do
    end subroutine InitializeLJ

    !Armchair cylindrical wall at r=radius
    subroutine GenerateWall1D()
        integer::atomsinlayer,atomsinheight,i,j
        real*8::temporary
        real*8,allocatable,dimension(:)::temporarytheta,eventheta,oddtheta  
        atomsinlayer=floor(pim2*radius/(3d0*WallBondLength))*2
        !Sometimes there is a round-off error due to double precision
        temporary=HalfLength(1)/(Sqrt(3d0)/2d0*WallBondLength)
        if(ceiling(temporary)-temporary<1d-14) then
            atomsinheight=2*ceiling(temporary)
        else
            atomsinheight=2*floor(temporary)
        end if
        allocate(wall(3,atomsinlayer*atomsinheight))
        allocate(temporarytheta(atomsinlayer/2))
        temporarytheta(1)=0
        temporary=pim2/size(temporarytheta)
        do i=2,size(temporarytheta)
            temporarytheta(i)=temporarytheta(i-1)+temporary
        end do
        allocate(eventheta(2*size(temporarytheta)))
        eventheta(1:size(temporarytheta))=temporarytheta
        eventheta((size(temporarytheta)+1):(2*size(temporarytheta)))=temporarytheta+temporary/3
        allocate(oddtheta(2*size(temporarytheta)))
        oddtheta=eventheta+temporary/2
        temporary=-HalfLength(1)
        do i=0,atomsinheight-1
            if(mod(i,2)==0) then
                do j=1,atomsinlayer
                    wall(2,i*atomsinlayer+j)=radius*Cos(eventheta(j))
                    wall(3,i*atomsinlayer+j)=radius*Sin(eventheta(j))
                    wall(1,i*atomsinlayer+j)=temporary
                end do
            else
                do j=1,atomsinlayer
                    wall(2,i*atomsinlayer+j)=radius*Cos(oddtheta(j))
                    wall(3,i*atomsinlayer+j)=radius*Sin(oddtheta(j))
                    wall(1,i*atomsinlayer+j)=temporary
                end do
            end if
            temporary=temporary+WallBondLength*Sqrt(3d0)/2d0
        end do
        deallocate(temporarytheta)
        deallocate(eventheta)
        deallocate(oddtheta)
    end subroutine GenerateWall1D
    !Slit wall in xy plane at top and bottom of the box
    !Not used for wall Lennard-Jones potential because steele potential is adopted
    !Steele potential accounts for all interactions from infinite plane and depth
    subroutine GenerateWall2D()
        integer::atomsinlayer,atomsinheight,i,j,k
        integer,dimension(2)::blin
        real*8::temporary
        real*8,dimension(2)::bl,bld2,centre
        real*8,allocatable,dimension(:)::temporarytheta,eventheta,oddtheta  
        bl=[Sqrt(3d0),3d0]*WallBondLength
        bld2=bl/2d0
        !Sometimes there is a round-off error due to double precision
        centre=HalfLength(1:2)/bl
        do i=1,2
            if(ceiling(centre(i))-centre(i)<1d-14) then
                blin(i)=2*ceiling(centre(i))
            else
                blin(i)=2*floor(centre(i))
            end if
        end do
        atomsinlayer=4*blin(1)*blin(2)
        allocate(wall(3,2*atomsinlayer))
        wall(3,1:atomsinlayer)=-HalfLength(3)
        wall(3,atomsinlayer+1:2*atomsinlayer)=HalfLength(3)
        centre=bld2-HalfLength(1:2)
        k=0
        do i=1,blin(1)
            do j=1,blin(2)
                k=k+1
                wall(1,k)=centre(1)
                wall(2,k)=centre(2)-WallBondLength/2d0
                k=k+1
                wall(1,k)=centre(1)
                wall(2,k)=centre(2)+WallBondLength/2d0
                k=k+1
                wall(1,k)=centre(1)-bld2(1)
                wall(2,k)=centre(2)-WallBondLength
                k=k+1
                wall(1,k)=centre(1)-bld2(1)
                wall(2,k)=centre(2)+WallBondLength
                centre(2)=centre(2)+bl(2)
            end do
            centre(1)=centre(1)+bl(1)
            centre(2)=bld2(2)-HalfLength(2)
        end do
        wall(1:2,atomsinlayer+1:2*atomsinlayer)=wall(1:2,1:atomsinlayer)
    end subroutine GenerateWall2D

    !Lennard-Jones potential between molecule and the wall
    !For 1D, explicit wall atom
    function WallLJPotential1D(tp,co)
        real*8::WallLJPotential1D
        integer,intent(in)::tp
        real*8,dimension(3,MoleculeInput(tp).NAtoms),intent(in)::co
        integer::i,j,iNon0Epsilon
        real*8::r2,r6,r12
        !Initialization: set the counter to zero
        WallLJPotential1D=0d0
        !Loop over molecular force centres with non-zero epsilon
        do iNon0Epsilon=1,MoleculeInput(tp).NNon0Epsilon
            i=MoleculeInput(tp).Non0Epsilon(iNon0Epsilon)
            !Loop over wall force centres
            do j=1,size(wall,2)
                DisplacementVector=co(:,i)-wall(:,j)
                !Periodic boundary condition
                if(DisplacementVector(1)>HalfLength(1)) then
                    DisplacementVector(1)=DisplacementVector(1)-BoxSize(1)
                else if(DisplacementVector(1)<-HalfLength(1)) then
                    DisplacementVector(1)=DisplacementVector(1)+BoxSize(1)
                end if
                r2=mtpWLJP(tp).atm(i).sigmasq/dot_product(DisplacementVector,DisplacementVector)
                r6=r2*r2*r2
                r12=r6*r6
                WallLJPotential1D=WallLJPotential1D+mtpWLJP(tp).atm(i).epsilonm4*(r12-r6)
            end do
        end do
    end function WallLJPotential1D
    !For 2D, 10-4-3 steele potential
    function WallLJPotential2D(tp,co)
        real*8::WallLJPotential2D
        integer,intent(in)::tp
        real*8,dimension(3,MoleculeInput(tp).NAtoms),intent(in)::co
        integer::i,iNon0Epsilon
        real*8::xtop,xbottom,x2,x4,sigma4
        !Initialization: set the counter to zero
        WallLJPotential2D=0d0
        !Loop over molecular force centres with non-zero epsilon
        do iNon0Epsilon=1,MoleculeInput(tp).NNon0Epsilon
            i=MoleculeInput(tp).Non0Epsilon(iNon0Epsilon)
            !Interaction from the top wall (without the coefficient and 4*epsilon*sigma^2)
                xtop=HalfLength(3)-co(3,i)
                sigma4=mtpWLJP(tp).atm(i).sigmasq*mtpWLJP(tp).atm(i).sigmasq
                x2=mtpWLJP(tp).atm(i).sigmasq/(xtop*xtop)
                x4=x2*x2
                xtop=xtop+0.61d0*LayerDistance
                xtop=xtop*xtop*xtop
                xtop=sigma4/(3d0*LayerDistance*xtop)
                xtop=0.4d0*x2*x4*x4-x4-xtop
            !Interaction from the bottom wall (without the coefficient and 4*epsilon*sigma^2)
                xbottom=co(3,i)+HalfLength(3)
                sigma4=mtpWLJP(tp).atm(i).sigmasq*mtpWLJP(tp).atm(i).sigmasq
                x2=mtpWLJP(tp).atm(i).sigmasq/(xbottom*xbottom)
                x4=x2*x2
                xbottom=xbottom+0.61d0*LayerDistance
                xbottom=xbottom*xbottom*xbottom
                xbottom=sigma4/(3d0*LayerDistance*xbottom)
                xbottom=0.4d0*x2*x4*x4-x4-xbottom
            !Times 4*epsilon*sigma^2
            WallLJPotential2D=WallLJPotential2D+mtpWLJP(tp).atm(i).epsilonm4*mtpWLJP(tp).atm(i).sigmasq*(xtop+xbottom)
        end do
        !Times the coefficient
        WallLJPotential2D=SteeleCoefficient*WallLJPotential2D
    end function WallLJPotential2D

    !Lennard-Jones potential between molecule1 and molecule2
    !For 1D
    function LJInter1D(tp1,tp2,co1,co2)
        real*8::LJInter1D
        integer,intent(in)::tp1,tp2
        real*8,dimension(3,MoleculeInput(tp1).NAtoms),intent(in)::co1
        real*8,dimension(3,MoleculeInput(tp2).NAtoms),intent(in)::co2
        integer::i,j,iNon0Epsilon,jNon0Epsilon
        real*8::epsilon,sig,r2,r6,r12
        !Initialization: set the counter to zero
        LJInter1D=0d0
        !Loop over molecule1 force centres with non-zero epsilon
        do iNon0Epsilon=1,MoleculeInput(tp1).NNon0Epsilon
            i=MoleculeInput(tp1).Non0Epsilon(iNon0Epsilon)
            !Loop over molecule2 force centres with non-zero epsilon
            do jNon0Epsilon=1,MoleculeInput(tp2).NNon0Epsilon
                j=MoleculeInput(tp2).Non0Epsilon(jNon0Epsilon)
                DisplacementVector=co2(:,j)-co1(:,i)
                !Periodic boundary condition
                if(DisplacementVector(1)>HalfLength(1)) then
                    DisplacementVector(1)=DisplacementVector(1)-BoxSize(1)
                else if(DisplacementVector(1)<-HalfLength(1)) then
                    DisplacementVector(1)=DisplacementVector(1)+BoxSize(1)
                end if
                r2=intertpLJP(tp1,tp2).atmab(i,j).sigmasq/dot_product(DisplacementVector,DisplacementVector)
                r6=r2*r2*r2
                r12=r6*r6
                LJInter1D=LJInter1D+intertpLJP(tp1,tp2).atmab(i,j).epsilonm4*(r12-r6)
            end do
        end do
    end function LJInter1D
    !For 2D
    function LJInter2D(tp1,tp2,co1,co2)
        real*8::LJInter2D
        integer,intent(in)::tp1,tp2
        real*8,dimension(3,MoleculeInput(tp1).NAtoms),intent(in)::co1
        real*8,dimension(3,MoleculeInput(tp2).NAtoms),intent(in)::co2
        integer::i,j,k,iNon0Epsilon,jNon0Epsilon
        real*8::epsilon,sig,r2,r6,r12
        !Initialization: set the counter to zero
        LJInter2D=0d0
        !Loop over molecule1 force centres with non-zero epsilon
        do iNon0Epsilon=1,MoleculeInput(tp1).NNon0Epsilon
            i=MoleculeInput(tp1).Non0Epsilon(iNon0Epsilon)
            !Loop over molecule2 force centres with non-zero epsilon
            do jNon0Epsilon=1,MoleculeInput(tp2).NNon0Epsilon
                j=MoleculeInput(tp2).Non0Epsilon(jNon0Epsilon)
                DisplacementVector=co2(:,j)-co1(:,i)
                !Periodic boundary condition
                do k=1,2
                    if(DisplacementVector(k)>HalfLength(k)) then
                        DisplacementVector(k)=DisplacementVector(k)-BoxSize(k)
                    else if(DisplacementVector(k)<-HalfLength(k)) then
                        DisplacementVector(k)=DisplacementVector(k)+BoxSize(k)
                    end if
                end do
                r2=intertpLJP(tp1,tp2).atmab(i,j).sigmasq/dot_product(DisplacementVector,DisplacementVector)
                r6=r2*r2*r2
                r12=r6*r6
                LJInter2D=LJInter2D+intertpLJP(tp1,tp2).atmab(i,j).epsilonm4*(r12-r6)
            end do
        end do
    end function LJInter2D
    !For 3D
    function LJInter3D(tp1,tp2,co1,co2)
        real*8::LJInter3D
        integer,intent(in)::tp1,tp2
        real*8,dimension(3,MoleculeInput(tp1).NAtoms),intent(in)::co1
        real*8,dimension(3,MoleculeInput(tp2).NAtoms),intent(in)::co2
        integer::i,j,k,iNon0Epsilon,jNon0Epsilon
        real*8::epsilon,sig,r2,r6,r12
        !Initialization: set the counter to zero
        LJInter3D=0d0
        !Loop over molecule1 force centres with non-zero epsilon
        do iNon0Epsilon=1,MoleculeInput(tp1).NNon0Epsilon
            i=MoleculeInput(tp1).Non0Epsilon(iNon0Epsilon)
            !Loop over molecule2 force centres with non-zero epsilon
            do jNon0Epsilon=1,MoleculeInput(tp2).NNon0Epsilon
                j=MoleculeInput(tp2).Non0Epsilon(jNon0Epsilon)
                DisplacementVector=co2(:,j)-co1(:,i)
                !Periodic boundary condition
                do k=1,3
                    if(DisplacementVector(k)>HalfLength(k)) then
                        DisplacementVector(k)=DisplacementVector(k)-BoxSize(k)
                    else if(DisplacementVector(k)<-HalfLength(k)) then
                        DisplacementVector(k)=DisplacementVector(k)+BoxSize(k)
                    end if
                end do
                r2=intertpLJP(tp1,tp2).atmab(i,j).sigmasq/dot_product(DisplacementVector,DisplacementVector)
                r6=r2*r2*r2
                r12=r6*r6
                LJInter3D=LJInter3D+intertpLJP(tp1,tp2).atmab(i,j).epsilonm4*(r12-r6)
            end do
        end do
    end function LJInter3D
!--------------------------------------------- End ---------------------------------------------------
    
!------------------------ Coulomb interaction: Ewald summation technique -----------------------------
    !The initializer for Ewald summation
    subroutine InitializeEwaldSummation()
        !ESSBIndependent is computed in subroutine EwaldSummationSingleBodyPositionIndependent
        allocate(ESSBIndependent(MoleculeKinds))
        select case(PeriodicDims)
            case(1)
                !Brodka: accuracy 1d-6 when Lx=3*Lx0, where Lx is the x component of the pseudo lattice vector,
                !Lx0 is the max x component absolute value of displacement vectors between to point charges
                !Similarly Ly 
                HalfLength(2)=3d0*radius
                HalfLength(3)=HalfLength(2)
                volume_vacuum=8d0*HalfLength(1)*HalfLength(2)*HalfLength(3)
                wavenumber=pi/HalfLength
                pim2dV=pim2/volume_vacuum
                pim4dV=pim4/volume_vacuum
                call EwaldSummationSingleBodyPositionIndependent1D()
                !chi, chisq, minus4chisq are computed in subroutine EwaldSummationDumpFactor
                call EwaldSummationDumpFactor1D()
                !For electro static force calculation
                Nonbonded_chim2dsqrtpi=chi*2d0/sqrtpi
                Nonbonded_chisqm2dLz=chisq/HalfLength(1)
            case(2)
                !2D PBC is Ewald summed by adding a vacuum layer in z direction then making z direction periodic
                !Spohr 1997 JCP: accuracy 1d-3 when the height of the vacuum layer = maxval(BoxSize(1:2))
                !So I add another BoxSize(3) similar to 1D case
                volume_vacuum=volume*2d0+BoxSize(1)*BoxSize(2)*maxval(BoxSize(1:2))
                wavenumber=pi/HalfLength
                wavenumber(3)=pi/(2d0*HalfLength(3)+maxval(HalfLength(1:2)))
                pim2dV=pim2/volume_vacuum
                pim4dV=pim4/volume_vacuum
                call EwaldSummationSingleBodyPositionIndependent2D()
                !chi, chisq, minus4chisq are computed in subroutine EwaldSummationDumpFactor
                call EwaldSummationDumpFactor2D()
                !For electro static force calculation
                Nonbonded_chim2dsqrtpi=chi*2d0/sqrtpi
            case(3)
                wavenumber=pi/HalfLength
                pim2dV=pim2/volume
                pim4dV=pim4/volume
                call EwaldSummationSingleBodyPositionIndependent3D()
                !chi, chisq, minus4chisq are computed in subroutine EwaldSummationDumpFactor
                call EwaldSummationDumpFactor3D()
                !Only 3D needs to store old ESSBIndependent, as ESSBIndependent changes only due to volume change
                allocate(ESSBIndependentOld(MoleculeKinds))
                ESSBIndependentOld=ESSBIndependent
        end select
    end subroutine InitializeEwaldSummation
    
    !Ewald summation dump factor and relevant terms. Call at start or when volume changes
    !For 1D
    subroutine EwaldSummationDumpFactor1D()
        !Choose the minimum chi satisfying accuracy. Smaller chi leads to faster convergence of reciprocal summation
        chi=inverse_erfc(reltol)/BoxSize(1)
        chisq=chi*chi
        minus4chisq=-4d0*chisq
    end subroutine EwaldSummationDumpFactor1D
    !For 2D
    subroutine EwaldSummationDumpFactor2D()
        !Choose the minimum chi satisfying accuracy. Smaller chi leads to faster convergence of reciprocal summation
        chi=inverse_erfc(reltol)/minval(BoxSize(1:2))
        chisq=chi*chi
        minus4chisq=-4d0*chisq
    end subroutine EwaldSummationDumpFactor2D
    !For 3D
    subroutine EwaldSummationDumpFactor3D()
        !Choose the minimum chi satisfying accuracy. Smaller chi leads to faster convergence of reciprocal summation
        chi=inverse_erfc(reltol)/minval(BoxSize)
        chisq=chi*chi
        minus4chisq=-4d0*chisq
    end subroutine EwaldSummationDumpFactor3D
    
    !The single-body position independent term. On exit, the dump factor has changed!
    !For 1D
    subroutine EwaldSummationSingleBodyPositionIndependent1D()
        integer::iKind,i,j,iNon0Charge,jNon0Charge
        real*8::reltolold,sccoefficient,r
        !This is computed once for all for 1D and 2D, once per volume change for 3D. We can afford the best accuracy
        reltolold=reltol
        reltol=1d-15
        call EwaldSummationDumpFactor1D()
        !Initialization: set the counter to zero
        ESSBIndependent=0d0
        !This is the coefficient for the single-charge term
        sccoefficient=pim2dV*ReciprocalSummationSingleCharge1D()-chi/sqrtpi
        !Loop over molecule kinds
        do iKind=1,MoleculeKinds
            !Single-charge
            ESSBIndependent(iKind)=sccoefficient*dot_product(MoleculeInput(iKind).charge,MoleculeInput(iKind).charge)
            !Intramolecular bicharge
            do iNon0Charge=1,MoleculeInput(iKind).NNon0Charge-1
                i=MoleculeInput(iKind).Non0Charge(iNon0Charge)
                do jNon0Charge=iNon0Charge+1,MoleculeInput(iKind).NNon0Charge
                    j=MoleculeInput(iKind).Non0Charge(jNon0Charge)
                    DisplacementVector=MoleculeInput(iKind).RefConfig(:,i)-MoleculeInput(iKind).RefConfig(:,j)
                    r=Norm2(DisplacementVector)
                    ESSBIndependent(iKind)=ESSBIndependent(iKind)+MoleculeInput(iKind).charge(i)*MoleculeInput(iKind).charge(j)&
                        *(pim4dV*ReciprocalSummation1D(DisplacementVector)-Erf(chi*r)/r&
                        -ES1D_Gz0(chisq*dot_product(DisplacementVector(2:3),DisplacementVector(2:3)))/BoxSize(1))
                end do
            end do
        end do
        !Recover accuracy
        reltol=reltolold
    end subroutine EwaldSummationSingleBodyPositionIndependent1D
    !For 2D
    subroutine EwaldSummationSingleBodyPositionIndependent2D()
        integer::iKind,i,j,iNon0Charge,jNon0Charge
        real*8::reltolold,sccoefficient,r
        !This is computed once for all for 1D and 2D, once per volume change for 3D. We can afford the best accuracy
        reltolold=reltol
        reltol=1d-15
        call EwaldSummationDumpFactor2D()
        !Initialization: set the counter to zero
        ESSBIndependent=0d0
        !This is the coefficient for the single-charge term
        sccoefficient=pim2dV*ReciprocalSummationSingleCharge()-chi/sqrtpi
        !Loop over molecule kinds
        do iKind=1,MoleculeKinds
            !Single-charge
            ESSBIndependent(iKind)=sccoefficient*dot_product(MoleculeInput(iKind).charge,MoleculeInput(iKind).charge)
            !Intramolecular bicharge
            do iNon0Charge=1,MoleculeInput(iKind).NNon0Charge-1
                i=MoleculeInput(iKind).Non0Charge(iNon0Charge)
                do jNon0Charge=iNon0Charge+1,MoleculeInput(iKind).NNon0Charge
                    j=MoleculeInput(iKind).Non0Charge(jNon0Charge)
                    DisplacementVector=MoleculeInput(iKind).RefConfig(:,i)-MoleculeInput(iKind).RefConfig(:,j)
                    r=Norm2(DisplacementVector)
                    ESSBIndependent(iKind)=ESSBIndependent(iKind)+MoleculeInput(iKind).charge(i)*MoleculeInput(iKind).charge(j)&
                        *(pim4dV*ReciprocalSummation(DisplacementVector)-Erf(chi*r)/r)
                end do
            end do
        end do
        !Recover accuracy
        reltol=reltolold
    end subroutine EwaldSummationSingleBodyPositionIndependent2D
    !For 3D
    subroutine EwaldSummationSingleBodyPositionIndependent3D()
        integer::iKind,i,j,iNon0Charge,jNon0Charge
        real*8::reltolold,sccoefficient,r
        !This is computed once for all for 1D and 2D, once per volume change for 3D. We can afford the best accuracy
        reltolold=reltol
        reltol=1d-15
        call EwaldSummationDumpFactor3D()
        !Initialization: set the counter to zero
        ESSBIndependent=0d0
        !This is the coefficient for the single-charge term
        sccoefficient=pim2dV*ReciprocalSummationSingleCharge()-chi/sqrtpi
        !Loop over molecule kinds
        do iKind=1,MoleculeKinds
            !Single-charge
            ESSBIndependent(iKind)=sccoefficient*dot_product(MoleculeInput(iKind).charge,MoleculeInput(iKind).charge)
            !Intramolecular bicharge
            do iNon0Charge=1,MoleculeInput(iKind).NNon0Charge-1
                i=MoleculeInput(iKind).Non0Charge(iNon0Charge)
                do jNon0Charge=iNon0Charge+1,MoleculeInput(iKind).NNon0Charge
                    j=MoleculeInput(iKind).Non0Charge(jNon0Charge)
                    DisplacementVector=MoleculeInput(iKind).RefConfig(:,i)-MoleculeInput(iKind).RefConfig(:,j)
                    r=Norm2(DisplacementVector)
                    ESSBIndependent(iKind)=ESSBIndependent(iKind)+MoleculeInput(iKind).charge(i)*MoleculeInput(iKind).charge(j)&
                        *(pim4dV*ReciprocalSummation(DisplacementVector)-Erf(chi*r)/r)
                end do
            end do
        end do
        !Recover accuracy
        reltol=reltolold
    end subroutine EwaldSummationSingleBodyPositionIndependent3D
    
    !For 2D only, the single-body position dependent term
    function EwaldSummationSingleBodyPositionDependent(tp,z)
        real*8::EwaldSummationSingleBodyPositionDependent
        integer,intent(in)::tp
        real*8,dimension(MoleculeInput(tp).NAtoms)::z
        integer::i,j,iNon0Charge,jNon0Charge
        real*8::temp,bicharge
        !Initialization: set the counters to zero
        EwaldSummationSingleBodyPositionDependent=0d0
        bicharge=0d0
        do iNon0Charge=1,MoleculeInput(tp).NNon0Charge
            i=MoleculeInput(tp).Non0Charge(iNon0Charge)
            !Single-charge
            temp=MoleculeInput(tp).charge(i)*z(i)
            EwaldSummationSingleBodyPositionDependent=EwaldSummationSingleBodyPositionDependent+temp*temp
            !Intramolecular bicharge
            do jNon0Charge=iNon0Charge+1,MoleculeInput(tp).NNon0Charge
                j=MoleculeInput(tp).Non0Charge(jNon0Charge)
                bicharge=bicharge+temp*MoleculeInput(tp).charge(j)*z(j)
            end do
        end do
        EwaldSummationSingleBodyPositionDependent=pim2dV*(EwaldSummationSingleBodyPositionDependent+2d0*bicharge)
    end function EwaldSummationSingleBodyPositionDependent
    
    !The many-body term of Ewald summation between molecule1 and molecule2
    !For 1D
    function EwaldSummationManyBody1D(tp1,tp2,co1,co2)
        real*8::EwaldSummationManyBody1D
        integer,intent(in)::tp1,tp2
        real*8,dimension(3,MoleculeInput(tp1).NAtoms),intent(in)::co1
        real*8,dimension(3,MoleculeInput(tp2).NAtoms),intent(in)::co2
        integer::i,j,iNon0Charge,jNon0Charge
        real*8::r
        !Initialization: set the counter to zero
        EwaldSummationManyBody1D=0d0
        !Loop over molecule1 force centres with non-zero charge
        do iNon0Charge=1,MoleculeInput(tp1).NNon0Charge
            i=MoleculeInput(tp1).Non0Charge(iNon0Charge)
            !Loop over molecule2 force centres with non-zero charge
            do jNon0Charge=1,MoleculeInput(tp2).NNon0Charge
                j=MoleculeInput(tp2).Non0Charge(jNon0Charge)
                DisplacementVector=co2(:,j)-co1(:,i)
                !Periodic boundary condition
                if(DisplacementVector(1)>HalfLength(1)) then
                    DisplacementVector(1)=DisplacementVector(1)-BoxSize(1)
                else if(DisplacementVector(1)<-HalfLength(1)) then
                    DisplacementVector(1)=DisplacementVector(1)+BoxSize(1)
                end if
                r=Norm2(DisplacementVector)
                EwaldSummationManyBody1D=EwaldSummationManyBody1D+MoleculeInput(tp1).charge(i)*MoleculeInput(tp2).charge(j)&
                    *(Erfc(chi*r)/r+pim4dV*ReciprocalSummation1D(DisplacementVector)&
                    -ES1D_Gz0(chisq*dot_product(DisplacementVector(2:3),DisplacementVector(2:3)))/BoxSize(1))
            end do
        end do
    end function EwaldSummationManyBody1D
    !For 2D
    function EwaldSummationManyBody2D(tp1,tp2,co1,co2)
        real*8::EwaldSummationManyBody2D
        integer,intent(in)::tp1,tp2
        real*8,dimension(3,MoleculeInput(tp1).NAtoms),intent(in)::co1
        real*8,dimension(3,MoleculeInput(tp2).NAtoms),intent(in)::co2
        integer::i,j,k,iNon0Charge,jNon0Charge
        real*8::r
        !Initialization: set the counter to zero
        EwaldSummationManyBody2D=0d0
        !Loop over molecule1 force centres with non-zero charge
        do iNon0Charge=1,MoleculeInput(tp1).NNon0Charge
            i=MoleculeInput(tp1).Non0Charge(iNon0Charge)
            !Loop over molecule2 force centres with non-zero charge
            do jNon0Charge=1,MoleculeInput(tp2).NNon0Charge
                j=MoleculeInput(tp2).Non0Charge(jNon0Charge)
                DisplacementVector=co2(:,j)-co1(:,i)
                !Periodic boundary condition
                do k=1,2
                    if(DisplacementVector(k)>HalfLength(k)) then
                        DisplacementVector(k)=DisplacementVector(k)-BoxSize(k)
                    else if(DisplacementVector(k)<-HalfLength(k)) then
                        DisplacementVector(k)=DisplacementVector(k)+BoxSize(k)
                    end if
                end do
                r=Norm2(DisplacementVector)
                EwaldSummationManyBody2D=EwaldSummationManyBody2D+MoleculeInput(tp1).charge(i)*MoleculeInput(tp2).charge(j)&
                    *(Erfc(chi*r)/r+pim4dV*ReciprocalSummation(DisplacementVector))
            end do
        end do
    end function EwaldSummationManyBody2D
    !For 3D
    function EwaldSummationManyBody3D(tp1,tp2,co1,co2)
        real*8::EwaldSummationManyBody3D
        integer,intent(in)::tp1,tp2
        real*8,dimension(3,MoleculeInput(tp1).NAtoms),intent(in)::co1
        real*8,dimension(3,MoleculeInput(tp2).NAtoms),intent(in)::co2
        integer::i,j,k,iNon0Charge,jNon0Charge
        real*8::r
        !Initialization: set the counter to zero
        EwaldSummationManyBody3D=0d0
        !Loop over molecule1 force centres with non-zero charge
        do iNon0Charge=1,MoleculeInput(tp1).NNon0Charge
            i=MoleculeInput(tp1).Non0Charge(iNon0Charge)
            !Loop over molecule2 force centres with non-zero charge
            do jNon0Charge=1,MoleculeInput(tp2).NNon0Charge
                j=MoleculeInput(tp2).Non0Charge(jNon0Charge)
                DisplacementVector=co2(:,j)-co1(:,i)
                !Periodic boundary condition
                do k=1,3
                    if(DisplacementVector(k)>HalfLength(k)) then
                        DisplacementVector(k)=DisplacementVector(k)-BoxSize(k)
                    else if(DisplacementVector(k)<-HalfLength(k)) then
                        DisplacementVector(k)=DisplacementVector(k)+BoxSize(k)
                    end if
                end do
                r=Norm2(DisplacementVector)
                EwaldSummationManyBody3D=EwaldSummationManyBody3D+MoleculeInput(tp1).charge(i)*MoleculeInput(tp2).charge(j)&
                    *(Erfc(chi*r)/r+pim4dV*ReciprocalSummation(DisplacementVector))
            end do
        end do
    end function EwaldSummationManyBody3D
    
    !The reciprocal space summation in single-charge term has higher symmetry
    !For 1D
    function ReciprocalSummationSingleCharge1D()
        !For 1D, kz=0 term is treated specially
        !This summation depends only on the 2-norm of the reciprocal vector, so according to the mirror symmetry,
        !we sum only: 1 of 2 half axes, 2 of 8 octant dividing planes, 1 of 8 octants
        !Namely: +z half axis, positive z/=0 octant dividing planes, positive octant
        real*8::ReciprocalSummationSingleCharge1D
        integer::k1norm,k1,k2
        real*8::axes,planes,octant,axeschange,planeschange,octantchange,kk
        !Initialization: set the counters to zero
        axes=0d0
        planes=0d0
        octant=0d0
        !Loop over 1-norm of the reciprocal lattice vector
        !Note that in 1D, the notation is z-x-y
        do k1norm=1,MaxIteration
            !Reset the counters for the changes of axes, planes, octant
            axeschange=0d0
            planeschange=0d0
            octantchange=0d0
            !+z half axis
                kk=k1norm*wavenumber(1)
                kk=kk*kk
                axeschange=axeschange+Exp(kk/minus4chisq)/kk
            do k1=1,k1norm-1
                !Positive z/=0 octant dividing planes
                    !x=0
                    DisplacementVector=[k1,0,k1norm-k1]*wavenumber
                    kk=dot_product(DisplacementVector,DisplacementVector)
                    planeschange=planeschange+Exp(kk/minus4chisq)/kk
                    !y=0
                    DisplacementVector=[k1,k1norm-k1,0]*wavenumber
                    kk=dot_product(DisplacementVector,DisplacementVector)
                    planeschange=planeschange+Exp(kk/minus4chisq)/kk
                !Positive octant
                do k2=1,k1norm-1-k1
                    DisplacementVector=[k1,k2,k1norm-k1-k2]*wavenumber
                    kk=dot_product(DisplacementVector,DisplacementVector)
                    octantchange=octantchange+Exp(kk/minus4chisq)/kk
                end do
            end do
            !Add the changes of axes, planes, octant
            axes=axes+axeschange
            planes=planes+planeschange
            octant=octant+octantchange
            !Convergence check
            axeschange=abs(axeschange)
            planeschange=abs(planeschange)
            octantchange=abs(octantchange)
            if( (axeschange<abstol.or.axeschange/abs(axes)<reltol).and.&
                (planeschange<abstol.or.planeschange/abs(planes)<reltol).and.&
                (octantchange<abstol.or.octantchange/abs(octant)<reltol) ) then
                exit
            end if
        end do
        !Times the symmetry
        ReciprocalSummationSingleCharge1D=2d0*axes+4d0*planes+8d0*octant
    end function ReciprocalSummationSingleCharge1D
    !For 2D and 3D
    function ReciprocalSummationSingleCharge()
        !This summation depends only on the 2-norm of the reciprocal vector, so according to the mirror symmetry,
        !we sum only: 3 of 6 half axes, 3 of 12 octant dividing planes, 1 of 8 octants
        !Namely: positive half axes, positive octant dividing planes, positive octant
        real*8::ReciprocalSummationSingleCharge
        integer::k1norm,k1,k2
        real*8::axes,planes,octant,axeschange,planeschange,octantchange,kk
        !Initialization: set the counters to zero
        axes=0d0
        planes=0d0
        octant=0d0
        !Loop over 1-norm of the reciprocal lattice vector
        do k1norm=1,MaxIteration
            !Reset the counters for the changes of axes, planes, octant
            axeschange=0d0
            planeschange=0d0
            octantchange=0d0
            !Positive half axes
                !+x
                kk=k1norm*wavenumber(1)
                kk=kk*kk
                axeschange=axeschange+Exp(kk/minus4chisq)/kk
                !+y
                kk=k1norm*wavenumber(2)
                kk=kk*kk
                axeschange=axeschange+Exp(kk/minus4chisq)/kk
                !+z
                kk=k1norm*wavenumber(3)
                kk=kk*kk
                axeschange=axeschange+Exp(kk/minus4chisq)/kk
            do k1=1,k1norm-1
                !Positive octant dividing planes
                    !x=0
                    DisplacementVector=[0,k1,k1norm-k1]*wavenumber
                    kk=dot_product(DisplacementVector,DisplacementVector)
                    planeschange=planeschange+Exp(kk/minus4chisq)/kk
                    !y=0
                    DisplacementVector=[k1,0,k1norm-k1]*wavenumber
                    kk=dot_product(DisplacementVector,DisplacementVector)
                    planeschange=planeschange+Exp(kk/minus4chisq)/kk
                    !z=0
                    DisplacementVector=[k1,k1norm-k1,0]*wavenumber
                    kk=dot_product(DisplacementVector,DisplacementVector)
                    planeschange=planeschange+Exp(kk/minus4chisq)/kk
                !Positive octant
                do k2=1,k1norm-1-k1
                    DisplacementVector=[k1,k2,k1norm-k1-k2]*wavenumber
                    kk=dot_product(DisplacementVector,DisplacementVector)
                    octantchange=octantchange+Exp(kk/minus4chisq)/kk
                end do
            end do
            !Add the changes of axes, planes, octant
            axes=axes+axeschange
            planes=planes+planeschange
            octant=octant+octantchange
            !Convergence check
            axeschange=abs(axeschange)
            planeschange=abs(planeschange)
            octantchange=abs(octantchange)
            if( (axeschange<abstol.or.axeschange/abs(axes)<reltol).and.&
                (planeschange<abstol.or.planeschange/abs(planes)<reltol).and.&
                (octantchange<abstol.or.octantchange/abs(octant)<reltol) ) then
                exit
            end if
        end do
        !Times the symmetry
        ReciprocalSummationSingleCharge=2d0*axes+4d0*planes+8d0*octant
    end function ReciprocalSummationSingleCharge
    
    !The usual bicharge reciprocal space summation
    !For 1D
    function ReciprocalSummation1D(r)
        !For 1D, kz=0 term is treated specially
        !This summation has point symmetry across the gamma point, so we sum only half of the space
        !Namely: +z half axis, four z>0 octant dividing planes, four z>0 octants
        real*8::ReciprocalSummation1D
        real*8,dimension(3),intent(in)::r
        integer::k1norm,k1,k2
        real*8::change,kk,temp
        real*8,dimension(3)::k
        !Initialization: set the counter to zero
        ReciprocalSummation1D=0d0
        !Loop over 1-norm of the reciprocal lattice vector
        !Note that in 1D, the notation is z-x-y
        do k1norm=1,MaxIteration
            !Reset the counter for the change of ReciprocalSummation
            change=0d0
            !+z half axis
                temp=k1norm*wavenumber(1)
                kk=temp*temp
                change=change+Exp(kk/minus4chisq)/kk*Cos(temp*r(1))
            do k1=1,k1norm-1
                !z>0 octant dividing planes
                    !+x+z
                    k=[k1,k1norm-k1,0]*wavenumber
                    kk=dot_product(k,k)
                    temp=dot_product(k,r)
                    change=change+Exp(kk/minus4chisq)/kk*Cos(temp)
                    !-x+z
                    temp=temp-2d0*k(2)*r(2)
                    k(2)=-k(2)
                    change=change+Exp(kk/minus4chisq)/kk*Cos(temp)
                    !+y+z
                    k=[k1,0,k1norm-k1]*wavenumber
                    kk=dot_product(k,k)
                    temp=dot_product(k,r)
                    change=change+Exp(kk/minus4chisq)/kk*Cos(temp)
                    !-y+z
                    temp=temp-2d0*k(3)*r(3)
                    k(3)=-k(3)
                    change=change+Exp(kk/minus4chisq)/kk*Cos(temp)
                !z>0 octants
                do k2=1,k1norm-1-k1
                    !+x+y+z
                    k=[k1,k2,k1norm-k1-k2]*wavenumber
                    kk=dot_product(k,k)
                    temp=dot_product(k,r)
                    change=change+Exp(kk/minus4chisq)/kk*Cos(temp)
                    !-x+y+z
                    temp=temp-2d0*k(2)*r(2)
                    k(2)=-k(2)
                    change=change+Exp(kk/minus4chisq)/kk*Cos(temp)
                    !-x-y+z
                    temp=temp-2d0*k(3)*r(3)
                    k(3)=-k(3)
                    change=change+Exp(kk/minus4chisq)/kk*Cos(temp)
                    !+x-y+z
                    temp=temp-2d0*k(2)*r(2)
                    k(2)=-k(2)
                    change=change+Exp(kk/minus4chisq)/kk*Cos(temp)
                end do
            end do
            !Add the change
            ReciprocalSummation1D=ReciprocalSummation1D+change
            !Convergence check
            change=abs(change)
            if(change<abstol.or.change/abs(ReciprocalSummation1D)<reltol) then
                exit
            end if
        end do
        !Times the symmetry
        ReciprocalSummation1D=ReciprocalSummation1D*2d0
    end function ReciprocalSummation1D
    !For 2D and 3D
    function ReciprocalSummation(r)
        !This summation has point symmetry across the gamma point, so we sum only half of the space
        !Namely: positive half axes, four x>0 and two z>0 octant dividing planes, four x>0 octants
        real*8::ReciprocalSummation
        real*8,dimension(3),intent(in)::r
        integer::k1norm,k1,k2
        real*8::change,kk,temp
        real*8,dimension(3)::k
        !Initialization: set the counter to zero
        ReciprocalSummation=0d0
        !Loop over 1-norm of the reciprocal lattice vector
        do k1norm=1,MaxIteration
            !Reset the counter for the change of ReciprocalSummation
            change=0d0
            !Positive half axes
                !+x
                temp=k1norm*wavenumber(1)
                kk=temp*temp
                change=change+Exp(kk/minus4chisq)/kk*Cos(temp*r(1))
                !+y
                temp=k1norm*wavenumber(2)
                kk=temp*temp
                change=change+Exp(kk/minus4chisq)/kk*Cos(temp*r(2))
                !+z
                temp=k1norm*wavenumber(3)
                kk=temp*temp
                change=change+Exp(kk/minus4chisq)/kk*Cos(temp*r(3))
            do k1=1,k1norm-1
                !four x>0 and two z>0 octant dividing planes
                    !+x+y
                    k=[k1,k1norm-k1,0]*wavenumber
                    kk=dot_product(k,k)
                    temp=dot_product(k,r)
                    change=change+Exp(kk/minus4chisq)/kk*Cos(temp)
                    !+x-y
                    temp=temp-2d0*k(2)*r(2)
                    k(2)=-k(2)
                    change=change+Exp(kk/minus4chisq)/kk*Cos(temp)
                    !+x+z
                    k=[k1,0,k1norm-k1]*wavenumber
                    kk=dot_product(k,k)
                    temp=dot_product(k,r)
                    change=change+Exp(kk/minus4chisq)/kk*Cos(temp)
                    !+x-z
                    temp=temp-2d0*k(3)*r(3)
                    k(3)=-k(3)
                    change=change+Exp(kk/minus4chisq)/kk*Cos(temp)
                    !+y+z
                    k=[0,k1,k1norm-k1]*wavenumber
                    kk=dot_product(k,k)
                    temp=dot_product(k,r)
                    change=change+Exp(kk/minus4chisq)/kk*Cos(temp)
                    !-y+z
                    temp=temp-2d0*k(2)*r(2)
                    k(2)=-k(2)
                    change=change+Exp(kk/minus4chisq)/kk*Cos(temp)
                !x>0 octants
                do k2=1,k1norm-1-k1
                    !+x+y+z
                    k=[k1,k2,k1norm-k1-k2]*wavenumber
                    kk=dot_product(k,k)
                    temp=dot_product(k,r)
                    change=change+Exp(kk/minus4chisq)/kk*Cos(temp)
                    !+x-y+z
                    temp=temp-2d0*k(2)*r(2)
                    k(2)=-k(2)
                    change=change+Exp(kk/minus4chisq)/kk*Cos(temp)
                    !+x-y-z
                    temp=temp-2d0*k(3)*r(3)
                    k(3)=-k(3)
                    change=change+Exp(kk/minus4chisq)/kk*Cos(temp)
                    !+x+y-z
                    temp=temp-2d0*k(2)*r(2)
                    k(2)=-k(2)
                    change=change+Exp(kk/minus4chisq)/kk*Cos(temp)
                end do
            end do
            !Add the change
            ReciprocalSummation=ReciprocalSummation+change
            !Convergence check
            change=abs(change)
            if(change<abstol.or.change/abs(ReciprocalSummation)<reltol) then
                exit
            end if
        end do
        !Times the symmetry
        ReciprocalSummation=ReciprocalSummation*2d0
    end function ReciprocalSummation
    
    !For 1D only, the kz=0 term of reciprocal space summation
    function ES1D_Gz0(x)
        integer::i
        real*8::ES1D_Gz0
        real*8,intent(in)::x!=chi**2*r**2
        !10 significant digits
        real*8::temp,temp2,temp3,temp6
        !When reltol=1d-10, chi**2>18.78741462316355d0, so in descending order
        if(x>18.78741462316355d0) then
            ES1D_Gz0=EulerGamma+log(x)
        else if(x>13.29d0) then
            temp=x-15.99d0
            temp2=temp*temp
            temp3=temp*temp2
            temp6=temp3*temp3
            ES1D_Gz0=3.3491791984585926d0+temp*0.06253907982075221d0-temp2*(0.0019555649204058842d0+temp2*3.823889781974916d-6)&
                +temp3*(8.153164493904012d-5+temp2*1.9125482223137583d-7)&
                +temp6*(-9.95754526846618d-9+temp*5.323631012903137d-10-temp2*2.8955510452797626d-11&
                +temp3*1.5900561610283544d-12-temp*temp3*8.753766072858468d-14)
        else if(x>9.29d0) then
            temp=x-11.26d0
            temp2=temp*temp
            temp3=temp*temp2
            temp6=temp3*temp3
            ES1D_Gz0=2.998473344163554d0+temp*0.08880880303082772d0-temp2*(0.003942980690850256d0+temp2*1.548920196156643d-5)&
                +temp3*(0.00023325998949853596d0+temp2*1.0909454673108656d-6)&
                +temp6*(-7.915055845502563d-8+temp*5.798241768265512d-9-temp2*4.222086908296498d-10&
                +temp3*3.01783879825018d-11-temp*temp3*2.096958494268975d-12)
        else if(x>6.24d0) then
            temp=x-7.74d0
            temp2=temp*temp
            temp3=temp*temp2
            temp6=temp3*temp3
            ES1D_Gz0=2.623667685122219d0+temp*0.12914275561045493d0-temp2*(0.008314449873086318+temp2*6.614403525922297d-5)&
                +temp3*(7.067770869351106d-4+temp2*6.368169476536038d-6)&
                +temp6*(-6.075635427151326d-7+temp*5.6129848279670643d-8-temp2*4.951311575317511d-9&
                +temp3*4.1372418821110844d-10-temp*temp3*3.2617274753582414d-11)
        else if(x>3.82d0) then
            temp=x-5.01d0
            temp2=temp*temp
            temp3=temp*temp2
            temp6=temp3*temp3
            ES1D_Gz0=2.1897864802195994d0+temp*0.19826928077719455d0-temp2*(0.01912159455797797d0+temp2*0.00029220534183235477d0)&
                +temp3*(0.002322537422676534d0+temp2*3.5563555405274084d-5)&
                +temp6*(-4.066098341347632d-6+temp*4.3146411654336176d-7-temp2*4.2331758731517724d-8&
                +temp3*3.841319120017832d-9-temp*temp3*3.231267642291151d-10)
        else if(x>1.85d0) then
            temp=x-2.81d0
            temp2=temp*temp
            temp3=temp*temp2
            temp6=temp3*temp3
            ES1D_Gz0=1.6270397285101232d0+temp*0.33444662192442226d0-temp2*(0.048797442977232874d0+temp2*0.0012441706978081175)&
                +temp3*(0.008006214680697388d0+temp2*0.00017566843000357574d0)&
                +temp6*(-2.233890154694154d-5+temp*2.5630586828795826d-6-temp2*2.667249013198523d-7&
                +temp3*2.533095435935187d-8-temp*temp3*2.2088886053649257d-9)
        else if(x>0.723d0) then
            temp=x-1.25d0
            temp2=temp*temp
            temp3=temp*temp2
            temp6=temp3*temp3
            ES1D_Gz0=0.9467725887416528d0+temp*0.570796162511848d0-temp2*(0.11371654626066321d0+temp2*0.003918751159257355d0)&
                +temp3*(0.022448185090995054d0+temp2*0.000597968762856782d0)&
                +temp6*(-8.030717872653526d-5+temp*9.590827815628838d-6-temp2*1.0289604856217327d-6&
                +temp3*1.0008090251673041d-7-temp*temp3*8.895816642345844d-9)
        else
            temp2=x*x
            temp3=x*temp2
            temp6=temp3*temp3
            ES1D_Gz0=x-temp2*(0.25d0+temp2/96d0)+temp3*(0.05555555555555555d0+temp2/600d0)&
                +temp6*(-2.314814814814815d-4+x/35280d0-temp2/322560d0+temp3/3265920d0-x*temp3/36288d3)
        end if
    end function ES1D_Gz0

    !The Coulomb force between force centres 1 and 2
    !For 1D
    function EwaldSummationManyBodyForce1D(tp1,tp2,atm1,atm2,co1,co2)
        real*8,dimension(3)::EwaldSummationManyBodyForce1D
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
        EwaldSummationManyBodyForce1D=MoleculeInput(tp1).charge(atm1)*MoleculeInput(tp2).charge(atm2)&
            *((Nonbonded_chim2dsqrtpi*Exp(-chisq*r2)/r+Erfc(chi*r)/r2)*runit&!Real space
            +pim4dV*ReciprocalForceSummation1D(DisplacementVector)&!Reciprocal space except Gz=0
            -Nonbonded_chisqm2dLz*(1-Exp(-chisqrhosq))/chisqrhosq*[rhounit(1),rhounit(2),0d0])!1D_Gz0_12
    end function EwaldSummationManyBodyForce1D
    !For 2D and 3D
    function EwaldSummationManyBodyForce(tp1,tp2,atm1,atm2,co1,co2)
        real*8,dimension(2)::EwaldSummationManyBodyForce
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
        EwaldSummationManyBodyForce=MoleculeInput(tp1).charge(atm1)*MoleculeInput(tp2).charge(atm2)&
            *((Nonbonded_chim2dsqrtpi*Exp(-chisq*r2)/r+Erfc(chi*r)/r2)*runit&
            +pim4dV*ReciprocalForceSummation(DisplacementVector))
    end function EwaldSummationManyBodyForce

    !Reciprocal space summation for force
    !For 1D
    function ReciprocalForceSummation1D(r)
        !For 1D, kz=0 term is treated specially
        !This summation has point symmetry across the gamma point, so we sum only half of the space
        !Namely: +z half axis, four z>0 octant dividing planes, four z>0 octants
        real*8,dimension(3)::ReciprocalForceSummation1D
        real*8,dimension(3),intent(in)::r
        integer::k1norm,k1,k2
        real*8::kk,temp
        real*8,dimension(3)::change,k
        !Initialization: set the counter to zero
        ReciprocalForceSummation1D=0d0
        !Loop over 1-norm of the reciprocal lattice vector
        !Note that in 1D, the notation is z-x-y
        do k1norm=1,MaxIteration
            !Reset the counter for the change of ReciprocalSummation
            change=0d0
            !+z half axis
                temp=k1norm*wavenumber(1)
                kk=temp*temp
                change(1)=change(1)+Exp(kk/minus4chisq)/kk*Sin(temp*r(1))*temp
            do k1=1,k1norm-1
                !z>0 octant dividing planes
                    !+x+z
                    k=[k1,k1norm-k1,0]*wavenumber
                    kk=dot_product(k,k)
                    temp=dot_product(k,r)
                    change=change+Exp(kk/minus4chisq)/kk*Sin(temp)*k
                    !-x+z
                    temp=temp-2d0*k(2)*r(2)
                    k(2)=-k(2)
                    change=change+Exp(kk/minus4chisq)/kk*Sin(temp)*k
                    !+y+z
                    k=[k1,0,k1norm-k1]*wavenumber
                    kk=dot_product(k,k)
                    temp=dot_product(k,r)
                    change=change+Exp(kk/minus4chisq)/kk*Sin(temp)*k
                    !-y+z
                    temp=temp-2d0*k(3)*r(3)
                    k(3)=-k(3)
                    change=change+Exp(kk/minus4chisq)/kk*Sin(temp)*k
                !z>0 octants
                do k2=1,k1norm-1-k1
                    !+x+y+z
                    k=[k1,k2,k1norm-k1-k2]*wavenumber
                    kk=dot_product(k,k)
                    temp=dot_product(k,r)
                    change=change+Exp(kk/minus4chisq)/kk*Sin(temp)*k
                    !-x+y+z
                    temp=temp-2d0*k(2)*r(2)
                    k(2)=-k(2)
                    change=change+Exp(kk/minus4chisq)/kk*Sin(temp)*k
                    !-x-y+z
                    temp=temp-2d0*k(3)*r(3)
                    k(3)=-k(3)
                    change=change+Exp(kk/minus4chisq)/kk*Sin(temp)*k
                    !+x-y+z
                    temp=temp-2d0*k(2)*r(2)
                    k(2)=-k(2)
                    change=change+Exp(kk/minus4chisq)/kk*Sin(temp)*k
                end do
            end do
            !Add the change
            ReciprocalForceSummation1D=ReciprocalForceSummation1D+change
            !Convergence check
            temp=dot_product(change,change)
            if(temp<abstol.or.temp/dot_product(ReciprocalForceSummation1D,ReciprocalForceSummation1D)<reltol) then
                exit
            end if
        end do
        !Times the symmetry
        ReciprocalForceSummation1D=ReciprocalForceSummation1D*2d0
    end function ReciprocalForceSummation1D
    !For 2D and 3D, x and y directions only
    function ReciprocalForceSummation(r)
        !This summation has point symmetry across the gamma point, so we sum only half of the space
        !Namely: positive half axes, four x>0 and two z>0 octant dividing planes, four x>0 octants
        !We do not care about p_z so +z half axis is ignored
        real*8,dimension(2)::ReciprocalForceSummation
        real*8,dimension(3),intent(in)::r
        integer::k1norm,k1,k2
        real*8::kk,temp
        real*8,dimension(2)::change
        real*8,dimension(3)::k
        !Initialization: set the counter to zero
        ReciprocalForceSummation=0d0
        !Loop over 1-norm of the reciprocal lattice vector
        do k1norm=1,MaxIteration
            !Reset the counter for the change of ReciprocalForceSummation
            change=0d0
            !Positive half axes
                !+x
                temp=k1norm*wavenumber(1)
                kk=temp*temp
                change(1)=change(1)+Exp(kk/minus4chisq)/kk*Sin(temp*r(1))*temp
                !+y
                temp=k1norm*wavenumber(2)
                kk=temp*temp
                change(2)=change(2)+Exp(kk/minus4chisq)/kk*Sin(temp*r(2))*temp
            do k1=1,k1norm-1
                !four x>0 and two z>0 octant dividing planes
                    !+x+y
                    k=[k1,k1norm-k1,0]*wavenumber
                    kk=dot_product(k,k)
                    temp=dot_product(k,r)
                    change=change+Exp(kk/minus4chisq)/kk*Sin(temp)*k(1:2)
                    !+x-y
                    temp=temp-2d0*k(2)*r(2)
                    k(2)=-k(2)
                    change=change+Exp(kk/minus4chisq)/kk*Sin(temp)*k(1:2)
                    !+x+z
                    k=[k1,0,k1norm-k1]*wavenumber
                    kk=dot_product(k,k)
                    temp=dot_product(k,r)
                    change=change+Exp(kk/minus4chisq)/kk*Sin(temp)*k(1:2)
                    !+x-z
                    temp=temp-2d0*k(3)*r(3)
                    k(3)=-k(3)
                    change=change+Exp(kk/minus4chisq)/kk*Sin(temp)*k(1:2)
                    !+y+z
                    k=[0,k1,k1norm-k1]*wavenumber
                    kk=dot_product(k,k)
                    temp=dot_product(k,r)
                    change=change+Exp(kk/minus4chisq)/kk*Sin(temp)*k(1:2)
                    !-y+z
                    temp=temp-2d0*k(2)*r(2)
                    k(2)=-k(2)
                    change=change+Exp(kk/minus4chisq)/kk*Sin(temp)*k(1:2)
                !x>0 octants
                do k2=1,k1norm-1-k1
                    !+x+y+z
                    k=[k1,k2,k1norm-k1-k2]*wavenumber
                    kk=dot_product(k,k)
                    temp=dot_product(k,r)
                    change=change+Exp(kk/minus4chisq)/kk*Sin(temp)*k(1:2)
                    !+x-y+z
                    temp=temp-2d0*k(2)*r(2)
                    k(2)=-k(2)
                    change=change+Exp(kk/minus4chisq)/kk*Sin(temp)*k(1:2)
                    !+x-y-z
                    temp=temp-2d0*k(3)*r(3)
                    k(3)=-k(3)
                    change=change+Exp(kk/minus4chisq)/kk*Sin(temp)*k(1:2)
                    !+x+y-z
                    temp=temp-2d0*k(2)*r(2)
                    k(2)=-k(2)
                    change=change+Exp(kk/minus4chisq)/kk*Sin(temp)*k(1:2)
                end do
            end do
            !Add the change
            ReciprocalForceSummation=ReciprocalForceSummation+change
            !Convergence check
            temp=dot_product(change,change)
            if(temp<abstol.or.temp/dot_product(ReciprocalForceSummation,ReciprocalForceSummation)<reltol) then
                exit
            end if
        end do
        !Times the symmetry
        ReciprocalForceSummation=ReciprocalForceSummation*2d0
    end function ReciprocalForceSummation
!---------------------------------------------- End --------------------------------------------------

!---------- Van de Waals interaction and Coulomb interaction can be computed in a same loop ----------
    !Compute all nonbonded many-body interaction between molecule1 and molecule2,
    !including: the many-body term of Ewald summation (returned in ESMB)
    !           the Lennard-Jones potential (returned in LJ)
    !The return value of the function is the total nonbonded many-body energy
    !For 1D
    function NonbondedManyBody1D(tp1,tp2,co1,co2,ESMB,LJ)
        real*8::NonbondedManyBody1D
        integer,intent(in)::tp1,tp2
        real*8,dimension(3,MoleculeInput(tp1).NAtoms),intent(in)::co1
        real*8,dimension(3,MoleculeInput(tp2).NAtoms),intent(in)::co2
        real*8,intent(out)::ESMB,LJ
        integer::i,j,iclass,jclass
        real*8::r,r2,r6,r12
        !Initialization: set the counters to zero
        ESMB=0d0
        LJ=0d0
        !Class: 1 epsilon only and 2 non-zero epsilon 
        do iclass=1,MoleculeInput(tp1).NEpsilonOnly
            i=MoleculeInput(tp1).EpsilonOnly(iclass)
            do jclass=1,MoleculeInput(tp2).NNon0Epsilon
                j=MoleculeInput(tp2).Non0Epsilon(jclass)
                DisplacementVector=co2(:,j)-co1(:,i)
                !Periodic boundary condition
                if(DisplacementVector(1)>HalfLength(1)) then
                    DisplacementVector(1)=DisplacementVector(1)-BoxSize(1)
                else if(DisplacementVector(1)<-HalfLength(1)) then
                    DisplacementVector(1)=DisplacementVector(1)+BoxSize(1)
                end if
                r2=intertpLJP(tp1,tp2).atmab(i,j).sigmasq/dot_product(DisplacementVector,DisplacementVector)
                r6=r2*r2*r2
                r12=r6*r6
                LJ=LJ+intertpLJP(tp1,tp2).atmab(i,j).epsilonm4*(r12-r6)
            end do
        end do
        !Class: 1 charge only and 2 non-zero charge
        do iclass=1,MoleculeInput(tp1).NChargeOnly
            i=MoleculeInput(tp1).ChargeOnly(iclass)
            do jclass=1,MoleculeInput(tp2).NNon0Charge
                j=MoleculeInput(tp2).Non0Charge(jclass)
                DisplacementVector=co2(:,j)-co1(:,i)
                !Periodic boundary condition
                if(DisplacementVector(1)>HalfLength(1)) then
                    DisplacementVector(1)=DisplacementVector(1)-BoxSize(1)
                else if(DisplacementVector(1)<-HalfLength(1)) then
                    DisplacementVector(1)=DisplacementVector(1)+BoxSize(1)
                end if
                r=Norm2(DisplacementVector)
                ESMB=ESMB+MoleculeInput(tp1).charge(i)*MoleculeInput(tp2).charge(j)&
                    *(Erfc(chi*r)/r+pim4dV*ReciprocalSummation1D(DisplacementVector)&
                    -ES1D_Gz0(chisq*dot_product(DisplacementVector(2:3),DisplacementVector(2:3)))/BoxSize(1))
            end do
        end do
        !Class: 1 both and...
        do iclass=1,MoleculeInput(tp1).NBoth
            i=MoleculeInput(tp1).Both(iclass)
            !2 epsilon only
            do jclass=1,MoleculeInput(tp2).NEpsilonOnly
                j=MoleculeInput(tp2).EpsilonOnly(jclass)
                DisplacementVector=co2(:,j)-co1(:,i)
                !Periodic boundary condition
                if(DisplacementVector(1)>HalfLength(1)) then
                    DisplacementVector(1)=DisplacementVector(1)-BoxSize(1)
                else if(DisplacementVector(1)<-HalfLength(1)) then
                    DisplacementVector(1)=DisplacementVector(1)+BoxSize(1)
                end if
                r2=intertpLJP(tp1,tp2).atmab(i,j).sigmasq/dot_product(DisplacementVector,DisplacementVector)
                r6=r2*r2*r2
                r12=r6*r6
                LJ=LJ+intertpLJP(tp1,tp2).atmab(i,j).epsilonm4*(r12-r6)
            end do
            !2 charge only
            do jclass=1,MoleculeInput(tp2).NChargeOnly
                j=MoleculeInput(tp2).ChargeOnly(jclass)
                DisplacementVector=co2(:,j)-co1(:,i)
                !Periodic boundary condition
                if(DisplacementVector(1)>HalfLength(1)) then
                    DisplacementVector(1)=DisplacementVector(1)-BoxSize(1)
                else if(DisplacementVector(1)<-HalfLength(1)) then
                    DisplacementVector(1)=DisplacementVector(1)+BoxSize(1)
                end if
                r=Norm2(DisplacementVector)
                ESMB=ESMB+MoleculeInput(tp1).charge(i)*MoleculeInput(tp2).charge(j)&
                    *(Erfc(chi*r)/r+pim4dV*ReciprocalSummation1D(DisplacementVector)&
                    -ES1D_Gz0(chisq*dot_product(DisplacementVector(2:3),DisplacementVector(2:3)))/BoxSize(1))
            end do
            !2 both
            do jclass=1,MoleculeInput(tp2).NBoth
                j=MoleculeInput(tp2).Both(jclass)
                DisplacementVector=co2(:,j)-co1(:,i)
                !Periodic boundary condition
                if(DisplacementVector(1)>HalfLength(1)) then
                    DisplacementVector(1)=DisplacementVector(1)-BoxSize(1)
                else if(DisplacementVector(1)<-HalfLength(1)) then
                    DisplacementVector(1)=DisplacementVector(1)+BoxSize(1)
                end if
                r2=dot_product(DisplacementVector,DisplacementVector)
                r=Sqrt(r2)
                ESMB=ESMB+MoleculeInput(tp1).charge(i)*MoleculeInput(tp2).charge(j)&
                    *(Erfc(chi*r)/r+pim4dV*ReciprocalSummation1D(DisplacementVector)&
                    -ES1D_Gz0(chisq*dot_product(DisplacementVector(2:3),DisplacementVector(2:3)))/BoxSize(1))
                r6=intertpLJP(tp1,tp2).atmab(i,j).sigmasq/r2
                r6=r6*r6*r6
                r12=r6*r6
                LJ=LJ+intertpLJP(tp1,tp2).atmab(i,j).epsilonm4*(r12-r6)
            end do
        end do
        NonbondedManyBody1D=ESMB+LJ
    end function NonbondedManyBody1D
    !For 2D
    function NonbondedManyBody2D(tp1,tp2,co1,co2,ESMB,LJ)
        real*8::NonbondedManyBody2D
        integer,intent(in)::tp1,tp2
        real*8,dimension(3,MoleculeInput(tp1).NAtoms),intent(in)::co1
        real*8,dimension(3,MoleculeInput(tp2).NAtoms),intent(in)::co2
        real*8,intent(out)::ESMB,LJ
        integer::i,j,iclass,jclass,k
        real*8::r,r2,r6,r12
        !Initialization: set the counters to zero
        ESMB=0d0
        LJ=0d0
        !Class: 1 epsilon only and 2 non-zero epsilon
        do iclass=1,MoleculeInput(tp1).NEpsilonOnly
            i=MoleculeInput(tp1).EpsilonOnly(iclass)
            do jclass=1,MoleculeInput(tp2).NNon0Epsilon
                j=MoleculeInput(tp2).Non0Epsilon(jclass)
                DisplacementVector=co2(:,j)-co1(:,i)
                !Periodic boundary condition
                do k=1,2
                    if(DisplacementVector(k)>HalfLength(k)) then
                        DisplacementVector(k)=DisplacementVector(k)-BoxSize(k)
                    else if(DisplacementVector(k)<-HalfLength(k)) then
                        DisplacementVector(k)=DisplacementVector(k)+BoxSize(k)
                    end if
                end do
                r2=intertpLJP(tp1,tp2).atmab(i,j).sigmasq/dot_product(DisplacementVector,DisplacementVector)
                r6=r2*r2*r2
                r12=r6*r6
                LJ=LJ+intertpLJP(tp1,tp2).atmab(i,j).epsilonm4*(r12-r6)
            end do
        end do
        !Class: 1 charge only and 2 non-zero charge
        do iclass=1,MoleculeInput(tp1).NChargeOnly
            i=MoleculeInput(tp1).ChargeOnly(iclass)
            do jclass=1,MoleculeInput(tp2).NNon0Charge
                j=MoleculeInput(tp2).Non0Charge(jclass)
                DisplacementVector=co2(:,j)-co1(:,i)
                !Periodic boundary condition
                do k=1,2
                    if(DisplacementVector(k)>HalfLength(k)) then
                        DisplacementVector(k)=DisplacementVector(k)-BoxSize(k)
                    else if(DisplacementVector(k)<-HalfLength(k)) then
                        DisplacementVector(k)=DisplacementVector(k)+BoxSize(k)
                    end if
                end do
                r=Norm2(DisplacementVector)
                ESMB=ESMB+MoleculeInput(tp1).charge(i)*MoleculeInput(tp2).charge(j)&
                    *(Erfc(chi*r)/r+pim4dV*ReciprocalSummation(DisplacementVector))
            end do
        end do
        !Class: 1 both and...
        do iclass=1,MoleculeInput(tp1).NBoth
            i=MoleculeInput(tp1).Both(iclass)
            !2 epsilon only
            do jclass=1,MoleculeInput(tp2).NEpsilonOnly
                j=MoleculeInput(tp2).EpsilonOnly(jclass)
                DisplacementVector=co2(:,j)-co1(:,i)
                !Periodic boundary condition
                do k=1,2
                    if(DisplacementVector(k)>HalfLength(k)) then
                        DisplacementVector(k)=DisplacementVector(k)-BoxSize(k)
                    else if(DisplacementVector(k)<-HalfLength(k)) then
                        DisplacementVector(k)=DisplacementVector(k)+BoxSize(k)
                    end if
                end do
                r2=intertpLJP(tp1,tp2).atmab(i,j).sigmasq/dot_product(DisplacementVector,DisplacementVector)
                r6=r2*r2*r2
                r12=r6*r6
                LJ=LJ+intertpLJP(tp1,tp2).atmab(i,j).epsilonm4*(r12-r6)
            end do
            !2 charge only
            do jclass=1,MoleculeInput(tp2).NChargeOnly
                j=MoleculeInput(tp2).ChargeOnly(jclass)
                DisplacementVector=co2(:,j)-co1(:,i)
                !Periodic boundary condition
                do k=1,2
                    if(DisplacementVector(k)>HalfLength(k)) then
                        DisplacementVector(k)=DisplacementVector(k)-BoxSize(k)
                    else if(DisplacementVector(k)<-HalfLength(k)) then
                        DisplacementVector(k)=DisplacementVector(k)+BoxSize(k)
                    end if
                end do
                r=Norm2(DisplacementVector)
                ESMB=ESMB+MoleculeInput(tp1).charge(i)*MoleculeInput(tp2).charge(j)&
                    *(Erfc(chi*r)/r+pim4dV*ReciprocalSummation(DisplacementVector))
            end do
            !2 both
            do jclass=1,MoleculeInput(tp2).NBoth
                j=MoleculeInput(tp2).Both(jclass)
                DisplacementVector=co2(:,j)-co1(:,i)
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
                ESMB=ESMB+MoleculeInput(tp1).charge(i)*MoleculeInput(tp2).charge(j)&
                    *(Erfc(chi*r)/r+pim4dV*ReciprocalSummation(DisplacementVector))
                r6=intertpLJP(tp1,tp2).atmab(i,j).sigmasq/r2
                r6=r6*r6*r6
                r12=r6*r6
                LJ=LJ+intertpLJP(tp1,tp2).atmab(i,j).epsilonm4*(r12-r6)
            end do
        end do
        NonbondedManyBody2D=ESMB+LJ
    end function NonbondedManyBody2D
    !For 3D
    function NonbondedManyBody3D(tp1,tp2,co1,co2,ESMB,LJ)
        real*8::NonbondedManyBody3D
        integer,intent(in)::tp1,tp2
        real*8,dimension(3,MoleculeInput(tp1).NAtoms),intent(in)::co1
        real*8,dimension(3,MoleculeInput(tp2).NAtoms),intent(in)::co2
        real*8,intent(out)::ESMB,LJ
        integer::i,j,iclass,jclass,k
        real*8::r,r2,r6,r12
        !Initialization: set the counters to zero
        ESMB=0d0
        LJ=0d0
        !Class: 1 epsilon only and 2 non-zero epsilon
        do iclass=1,MoleculeInput(tp1).NEpsilonOnly
            i=MoleculeInput(tp1).EpsilonOnly(iclass)
            do jclass=1,MoleculeInput(tp2).NNon0Epsilon
                j=MoleculeInput(tp2).Non0Epsilon(jclass)
                DisplacementVector=co2(:,j)-co1(:,i)
                !Periodic boundary condition
                do k=1,3
                    if(DisplacementVector(k)>HalfLength(k)) then
                        DisplacementVector(k)=DisplacementVector(k)-BoxSize(k)
                    else if(DisplacementVector(k)<-HalfLength(k)) then
                        DisplacementVector(k)=DisplacementVector(k)+BoxSize(k)
                    end if
                end do
                r2=intertpLJP(tp1,tp2).atmab(i,j).sigmasq/dot_product(DisplacementVector,DisplacementVector)
                r6=r2*r2*r2
                r12=r6*r6
                LJ=LJ+intertpLJP(tp1,tp2).atmab(i,j).epsilonm4*(r12-r6)
            end do
        end do
        !Class: 1 charge only and 2 non-zero charge
        do iclass=1,MoleculeInput(tp1).NChargeOnly
            i=MoleculeInput(tp1).ChargeOnly(iclass)
            do jclass=1,MoleculeInput(tp2).NNon0Charge
                j=MoleculeInput(tp2).Non0Charge(jclass)
                DisplacementVector=co2(:,j)-co1(:,i)
                !Periodic boundary condition
                do k=1,3
                    if(DisplacementVector(k)>HalfLength(k)) then
                        DisplacementVector(k)=DisplacementVector(k)-BoxSize(k)
                    else if(DisplacementVector(k)<-HalfLength(k)) then
                        DisplacementVector(k)=DisplacementVector(k)+BoxSize(k)
                    end if
                end do
                r=Norm2(DisplacementVector)
                ESMB=ESMB+MoleculeInput(tp1).charge(i)*MoleculeInput(tp2).charge(j)&
                    *(Erfc(chi*r)/r+pim4dV*ReciprocalSummation(DisplacementVector))
            end do
        end do
        !Class: 1 both and...
        do iclass=1,MoleculeInput(tp1).NBoth
            i=MoleculeInput(tp1).Both(iclass)
            !2 epsilon only
            do jclass=1,MoleculeInput(tp2).NEpsilonOnly
                j=MoleculeInput(tp2).EpsilonOnly(jclass)
                DisplacementVector=co2(:,j)-co1(:,i)
                !Periodic boundary condition
                do k=1,3
                    if(DisplacementVector(k)>HalfLength(k)) then
                        DisplacementVector(k)=DisplacementVector(k)-BoxSize(k)
                    else if(DisplacementVector(k)<-HalfLength(k)) then
                        DisplacementVector(k)=DisplacementVector(k)+BoxSize(k)
                    end if
                end do
                r2=intertpLJP(tp1,tp2).atmab(i,j).sigmasq/dot_product(DisplacementVector,DisplacementVector)
                r6=r2*r2*r2
                r12=r6*r6
                LJ=LJ+intertpLJP(tp1,tp2).atmab(i,j).epsilonm4*(r12-r6)
            end do
            !2 charge only
            do jclass=1,MoleculeInput(tp2).NChargeOnly
                j=MoleculeInput(tp2).ChargeOnly(jclass)
                DisplacementVector=co2(:,j)-co1(:,i)
                !Periodic boundary condition
                do k=1,3
                    if(DisplacementVector(k)>HalfLength(k)) then
                        DisplacementVector(k)=DisplacementVector(k)-BoxSize(k)
                    else if(DisplacementVector(k)<-HalfLength(k)) then
                        DisplacementVector(k)=DisplacementVector(k)+BoxSize(k)
                    end if
                end do
                r=Norm2(DisplacementVector)
                ESMB=ESMB+MoleculeInput(tp1).charge(i)*MoleculeInput(tp2).charge(j)&
                    *(Erfc(chi*r)/r+pim4dV*ReciprocalSummation(DisplacementVector))
            end do
            !2 both
            do jclass=1,MoleculeInput(tp2).NBoth
                j=MoleculeInput(tp2).Both(jclass)
                DisplacementVector=co2(:,j)-co1(:,i)
                !Periodic boundary condition
                do k=1,3
                    if(DisplacementVector(k)>HalfLength(k)) then
                        DisplacementVector(k)=DisplacementVector(k)-BoxSize(k)
                    else if(DisplacementVector(k)<-HalfLength(k)) then
                        DisplacementVector(k)=DisplacementVector(k)+BoxSize(k)
                    end if
                end do
                r2=dot_product(DisplacementVector,DisplacementVector)
                r=Sqrt(r2)
                ESMB=ESMB+MoleculeInput(tp1).charge(i)*MoleculeInput(tp2).charge(j)&
                    *(Erfc(chi*r)/r+pim4dV*ReciprocalSummation(DisplacementVector))
                r6=intertpLJP(tp1,tp2).atmab(i,j).sigmasq/r2
                r6=r6*r6*r6
                r12=r6*r6
                LJ=LJ+intertpLJP(tp1,tp2).atmab(i,j).epsilonm4*(r12-r6)
            end do
        end do
        NonbondedManyBody3D=ESMB+LJ
    end function NonbondedManyBody3D
!---------------------------------------------- End --------------------------------------------------

!---------------- These are just functions for reading. They are not actually called -----------------
    !Lennard-Jones potential
    function LJ(e,s,rr)
        real*8::LJ
        real*8,intent(in)::e,s,rr
        real*8::r6,r12
        r6=s*s/rr
        r6=r6*r6*r6
        r12=r6*r6
        LJ=4d0*e*(r12-r6)
    end function LJ
    
    !The force of Lennard-Jones potential
    function LJForce(e,s,r)
        real*8,dimension(3)::LJForce
        real*8,intent(in)::e,s
        real*8,dimension(3),intent(in)::r
        real*8::r2,r6,r12
        r2=dot_product(r,r)
        r6=s*s/r2
        r6=r6*r6*r6
        r12=r6*r6
        LJForce=24d0*e*(2d0*r12-r6)*r/r2
    end function LJForce
    
    !Lorentz-Berthelot mixing rule for Lennard-Jones potential between different type of force centers
    subroutine lbmix(e1,s1,e2,s2,e,s)
        real*8::e1,s1,e2,s2,e,s
        e=Sqrt(e1*e2)
        s=(s1+s2)/2d0
    end subroutine lbmix
    
    !Coulomb potential between two point charges
    function Coulomb(q1,q2,r)
        real*8::Coulomb
        real*8,intent(in)::q1,q2,r
        Coulomb=q1*q2/r
    end function Coulomb
    
    !Coulomb force between two point charges
    function CoulombForce(q1,q2,r)
        real*8,dimension(3)::CoulombForce
        real*8,intent(in)::q1,q2
        real*8,dimension(3),intent(in)::r
        real*8::rr
        rr=dot_product(r,r)
        CoulombForce=q1*q2*r/(rr*sqrt(rr))
    end function CoulombForce
!---------------------------------------------- End --------------------------------------------------

end module Nonbonded