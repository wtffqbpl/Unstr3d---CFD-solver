subroutine bcond(ghocell,adcell,adflag,iface)
    implicit none 
!**********************************************************************
!
! Function: Update boundary conditions
!
!**********************************************************************

! Arguments
    integer :: ghocell,adcell,adflag,iface

! Local variables
    integer :: i

! All kinds of boundary conditions
    selectcase(adflag)
        case(11)   ! Farfield
            call farfield(ghocell,adcell,iface)
        case(12)   ! Symmetry plane
            call symmetry(ghocell,adcell,iface)
        case(13)   ! inviscid wall
            call inviscid(ghocell,adcell,iface)
        case(14)   ! no-slip wall
            call noslipwal(ghocell,adcell,iface)
    end select

    return 
end subroutine bcond


subroutine farfield(ghocell, adcell,iface)
    use globalvar
    implicit none 
!**********************************************************************
!
! Function: Farfield boundary condition
!
!**********************************************************************

! Arguments
    integer :: ghocell,adcell,iface

! Local variables
    real*8  :: rieminf,rieminn,norm(ndimn),enorm,ex,ey,ez,vari(nmet),varb(nmet),&
               pp,vn0,cbond,entropy0,tempa,tempma,temp1,temp2,temp3, &
               temp4,preangle

    norm(:)=cellinfo%facenorm(:,iface,adcell)
    if(3 .eq. ndimn) then
        enorm=dsqrt(norm(1)**2+norm(2)**2+norm(3)**2)
        ex   =norm(1)/enorm
        ey   =norm(2)/enorm
        ey   =norm(3)/enorm
        call ctop(vari(:),cellinfo%cvar(:,adcell))

        vn0   =vari(2)*ex+vari(3)*ey+vari(4)*ez
        tempa =dsqrt(gamma*vari(5)/vari(1))
        tempma=dsqrt(vari(2)**2+vari(3)**2+vari(4)**2)/tempa

        if(tempma .ge. 1.0D0) then   ! Supersonic boundary conditions
            if(vn0 .lt. 0.0D0) then  ! inflow
                varb(1)=rhoinf
                varb(2)=uinf
                varb(3)=vinf
                varb(4)=winf
                varb(5)=pinf
            else                     ! Outflow
                varb(1)=vari(1)
                varb(2)=vari(2)
                varb(3)=vari(3)
                varb(4)=vari(4)
                varb(5)=vari(5)
            end if
        else                         ! Subsonic boundary 
            rieminn=vn0-2.0D0*tempa/gam1
            rieminf=uinf*ex+vinf*ey+winf*ez+2.0D0*dsqrt(gamma*rgas*tinf)/gam1
            temp1  =0.5D0*(rieminn+rieminf)
            temp2  =0.5D0*(rieminn-rieminf)/gam1
            if(temp1 .lt. 0.0D0) then !inflow
                temp4=rhoinf**gamma/pinf
                temp3=(temp2**2+temp4)**(1.0D0/gam1)
                varb(5)=0.5D0*(pinf+vari(5)-temp3*dabs(temp2)*(ex*uinf-vari(2))+ &
                        ey*(vinf-vari(3))+ez(vinf-vari(4)))
                varb(1)=rhoinf+(varb(5)-pinf)/temp2**2
                varb(2)=uinf-ex*(pinf-varb(5))/(temp3*dabs(temp2))
                varb(3)=uinf-ey*(pinf-varb(5))/(temp3*dabs(temp2))
                varb(4)=uinf-ez*(pinf-varb(5))/(temp3*dabs(temp2))
            else    
                temp4=vari(i)**gamma/vari(5)
                temp3=(temp2**2*temp4)**(1.0D0/gam1)
                varb(5)=pinf
                varb(1)=vari(1)+(varb(5)-vari(5))/temp2**2
                varb(2)=varb(2)+ex*(vari(5)-varb(5))/(temp3*dabs(temp2))
                varb(3)=varb(3)+ey*(vari(5)-varb(5))/(temp3*dabs(temp2))
                varb(4)=varb(4)+ez*(vari(5)-varb(5))/(temp3*dabs(temp2))
            end if
        end if
    else if( 2 .eq. ndimn) then
        enorm=dsqrt(norm(1)**2+norm(2)**2)
        ex   =norm(1)/enorm
        ey   =norm(2)/enorm
        call ctop(vari(:),cellinfo%cvar(:,adcell))

        vn0   =vari(2)*ex+vari(3)*ey
        tempa =dsqrt(gamma*vari(4)/vari(1))
        tempma=dsqrt(vari(2)**2+vari(3)**2)/tempa

        if(tempma .ge. 1.0D0) then  ! Supersonic boundary conditions
            if(vn0 .lt. 0.0D0) then ! Inflow
                varb(1)=rhoinf
                varb(2)=uinf
                varb(3)=vinf
                varb(4)=pinf
            else                    ! Outflow
                varb(1)=vari(1)
                varb(2)=vari(2)
                varb(3)=vari(3)
                varb(4)=vari(4)
            end if
        else                        ! Subsonic boundary
            rieminn=vn0-2.0D0*tempa/gam1
            rieminf=uinf*ex+vinf*ey+2.0D0*dsqrt(gamma*rgas*tinf)/gam1
            temp1  =0.50D0*(rieminn+rieminf)
            temp2  =0.25D0*(rieminn-rieminf)*gam1
            if(temp1 .lt. 0.0D0) then ! inflow
                temp4=rhoinf**gamma/pinf
                temp3=(temp2**2*temp4)**(1.0D0/gam1)
                varb(4)=0.50D0*(pinf+vari(4)-temp3*dabs(temp2)*(ex*(uinf-vari(2))+ &
                            ey*(vinf-vari(3)))
                varb(1)=rhoinf+(varb(4)-pinf)/temp2**2
                varb(2)=uinf-ex*(pinf-varb(4))/(temp3*dabs(temp2))
                varb(3)=vinf-ey*(pinf-varb(4))/(temp3*dabs(temp2))
            else
                temp4=vari(1)**gamma/vari(4)
                temp3=(temp2**2*temp4)**(1.0D0/gam1)
                varb(4)=pinf
                varb(1)=vari(1)+(varb(4)-vari(4))/temp2**2
                varb(2)=vari(2)+ex*(vari(4)-varb(4))/(temp3*dabs(temp2))
                varb(3)=vari(3)+ex*(vari(4)-varb(4))/(temp3*dabs(temp2))
            end if
        end if
    end if

    call ptoc(varb(:),cellinfo%cvar(:,ghocell))

    return 
end subroutine farfield



subroutine symmetry(ghocell,adcell,iface)
    use globalvar
    implicit none 
!**********************************************************************
!
! Function: Symmetry boundary conditions
!
!**********************************************************************

! Arguments
    integer :: ghocell,adcell,iface

! Local variables
    real*8  :: norm(ndim),enorm,ex,ey,ez,uu,cvar0(nmet),cvar1(nmet)

    norm(:) = cellinfo%facenorm(:,iface,adcell)
    if(3 .eq. ndimn) then
        enorm =dsqrt(norm(1)**2+norm(2)**2+norm(3)**2)
        ex    =norm(1)/enorm
        ey    =norm(2)/enorm
        ez    =norm(3)/enorm
        call ctop(cvar0(:),cellinfo%cvar(:,adcell))

        uu    =cvar0(2)*ex+cvar0(3)*ey+cvar0(4)*ez
        cvar1(1)=cvar0(1)
        cvar1(2)=cvar0(2)-2.0D0*uu*ex
        cvar1(3)=cvar0(3)-2.0D0*uu*ey
        cvar1(4)=cvar0(4)-2.0D0*uu*ez
        cvar1(5)=cvar0(5)
    else if( 2 .eq. ndimn) then
        enorm =dsqrt(norm(1)**2+norm(2)**2)
        ex    =norm(1)/enorm
        ex    =norm(2)/enorm
        call ctop(cvar0(:),cellinfo%cvar(:,adcell))
        
        uu    =cvar0(2)*ex+cvar0(3)*ey
        cvar1(1)=cvar0(1)
        cvar1(2)=cvar0(2)-2.0D0*uu*ex
        cvar1(3)=cvar0(3)-2.0D0*uu*ey
        cvar1(4)=cvar0(4)
    end if
    call ptoc(cvar1(:),cellinfo%cvar(:,ghocell))

    return 
end subroutine symmetry



subroutine inviscid(ghocell,adcell,iface)
    use globalvar
    implicit none 
!**********************************************************************
!
! Function: Inviscid wall boundary condition
!
!**********************************************************************

! Arguments
    integer :: ghocell,adcell,iface

! Local variables
    real*8  :: norm(ndimn),enorm,ex,ey,ez,uu,cvar0(nmet),cvar1(nmet)

    norm(:) = cellinfo%facenorm(:,iface,adcell)
    if(3 .eq. ndimn) then
        enorm =dsqrt(norm(1)**2+norm(2)**2+norm(3)**2)
        ex    =norm(1)/enorm
        ey    =norm(2)/enorm
        ez    =norm(3)/enorm

        call ctop(cvar0(:),cellinfo%cvar(:,adcell))
        uu    =cvar0(2)*ex+cvar0(3)*ey+cvar0(4)*ez
        cvar1(1)=cvar0(1)
        cvar1(2)=cvar0(2)-2.0D0*uu*ex
        cvar1(3)=cvar0(3)-2.0D0*uu*ey
        cvar1(4)=cvar0(4)-2.0D0*uu*ez
        cvar1(5)=cvar0(5)
    else if( 2 .eq. ndimn) then
        enorm =dsqrt(norm(1)**2+norm(2)**2)
        ex    =norm(1)/enorm
        ex    =norm(2)/enorm

        call ctop(cvar0(:),cellinfo%cvar(:,adcell))
        uu    =cvar0(2)*ex+cvar0(3)*ey
        cvar1(1)=cvar0(1)
        cvar1(2)=cvar0(2)-2.0D0*uu*ex
        cvar1(3)=cvar0(3)-2.0D0*uu*ey
        cvar1(4)=cvar0(4)
    end if
    call ptoc(cvar1(:),cellinfo%cvar(:,ghocell))

    return 
end subroutine inviscid



subroutine noslipwal(ghocell,adcell,iface)
    use globalvar
    implicit none 
!**********************************************************************
!
! Function: No-slip wall boundary condition
!
!**********************************************************************

! Arguments
    integer :: ghocell,adcell,iface

! Local variables
    real*8  :: cvar0(nmet),cvar1(nmet)

    cvar0(:) = cellinfo%cvar(:,adcell)
    if(3 .eq. ndmin) then
        cvar1(1)= cvar0(1)
        cvar1(2)=-cvar0(2)
        cvar1(3)=-cvar0(3)
        cvar1(4)=-cvar0(4)
        cvar1(5)= cvar0(5)
    else if( 2 .eq. ndmin) then
        cvar1(1)= cvar0(1)
        cvar1(2)=-cvar0(2)
        cvar1(3)=-cvar0(3)
        cvar1(4)= cvar0(4)
    end if
    
    cellinfo%cvar(:,ghocell)=cvar1(:)

    return 
end subroutine noslipwal
