!General mathematical and physical constant
!General basic routines and math routines
module General
    implicit none

!Constant
    !Mathematical constant
    real*8,parameter::pi=3.141592653589793d0,EulerGamma=0.57721566490153286060651209d0
        !Frequently used form (m is short for multiply, d is short for divide)
        real*8,parameter::pim2=6.283185307179586d0,pim4=12.566370614359172d0,sqrtpi=1.7724538509055159d0
    !Physical constant
    real*8,parameter::NAvogadro=6.02214076d23
    !Unit conversion
    real*8,parameter::KInAU=3.166813539739535d-6,barInAU=3.39882737736419d-9,&
        AInAU=1.8897261339212517d0,kJmolInAU=0.00038087967507991464d0,&
        AMUInAU=1822.888486192d0

contains
!------------ System general basic subroutine -------------
    !Show date hour minute second
    subroutine ShowTime()
        integer,dimension(8)::value
        call date_and_time(values=value)
        write(*,*)value(3),'d',value(5),':',value(6),':',value(7)
    end subroutine ShowTime

    !A better version of random_seed
    subroutine BetterRandomSeed()
        integer ii,nn,value(1:8)
        integer,allocatable :: seed(:)
        double precision flagd
        call random_seed(size=nn)
        allocate(seed(nn))
        call date_and_time(values=value)
        seed = value(8)+37*(/(ii-1,ii=1,nn)/)
        call random_seed(put=seed)
        deallocate(seed)
        do ii=1,value(6)*3600+value(7)*60+value(8)
            call random_number(flagd)
        end do
    end subroutine BetterRandomSeed
!-------------------------- End ---------------------------

!---------- Mathematical general basic functions ----------
    !Get a dimension(i) array full of 1d0
    function ones(i)
        integer,intent(in)::i
        real*8,dimension(i)::ones
        ones=1d0
    end function ones
    
    function inverse_erfc(x)
        logical::flag
        real*8::inverse_erfc,x,temp
        flag=0
        if(x>1d0) then
            x=2d0-x
            flag=1
        end if
        if(x<0.1d0) then
            temp=log(2d0/pi)-2d0*log(x)
            inverse_erfc=sqrt((temp-log(temp))/2d0)
        else if(x<0.25d0) then
            temp=4.995264535887506d0*(x-0.15)
            inverse_erfc=1.0179024648320276d0-temp/2d0+temp**2*0.2544756162080069d0&
                -temp**3*0.21435423798518619d0+temp**4*0.20605638309556853d0
        else 
            temp=(1d0-x)*sqrtpi
            inverse_erfc=temp/2d0+temp**3/24d0+temp**5*7d0/960d0+temp**7*127d0/80640d0&
                +temp**9*4369d0/11612160d0
        end if
        if(flag) then
            x=2d0-x
            inverse_erfc=-inverse_erfc
        end if
    end function inverse_erfc
    
    !Quaternion's multiplication, return the product of a * b
    function quamul(a,b)
        real*8,dimension(4)::quamul
        real*8,dimension(4),intent(in)::a,b
        quamul(1)=a(1)*b(1)-a(2)*b(2)-a(3)*b(3)-a(4)*b(4)
        quamul(2)=a(2)*b(1)+a(1)*b(2)+a(4)*b(3)-a(3)*b(4)
        quamul(3)=a(3)*b(1)+a(1)*b(3)+a(2)*b(4)-a(4)*b(2)
        quamul(4)=a(4)*b(1)+a(1)*b(4)+a(3)*b(2)-a(2)*b(3)
    end function quamul
    
    !Rotation in 3 dimensional space by quaternion rather than matrix
    !q=(cos(theta/2),sin(theta/2)*axis)
    function rotate(q,co)
        real*8,dimension(3)::rotate
        real*8,dimension(4),intent(in)::q
        real*8,dimension(3),intent(in)::co
        real*8,dimension(4)::qstar,qtemp
        qstar(1)=q(1)
        qstar(2:4)=-q(2:4)
        qtemp(1)=0d0
        qtemp(2:4)=co
        qtemp=quamul(quamul(qstar,qtemp),q)
        rotate=qtemp(2:4)
    end function rotate
    
    !Return a unit quaternion
    function UnitQuaternion()
        real*8,dimension(4)::UnitQuaternion
        real*8,dimension(2)::s1,s2
        do while(1)
            call random_number(s1)
            s1=s1*2d0-1d0
            if(dot_product(s1,s1)<1d0) then
                exit
            end if
        end do
        do while(1)
            call random_number(s2)
            s2=s2*2d0-1d0
            if(dot_product(s2,s2)<1d0) then
                exit
            end if
        end do
        UnitQuaternion(1:2)=s1
        UnitQuaternion(3:4)=s2*Sqrt((1d0-dot_product(s1,s1))/dot_product(s2,s2))
    end function UnitQuaternion
!-------------------------- End ---------------------------

end module General