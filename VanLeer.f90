subroutine VanLeer(RHS)
    use globalvar
    implicit none 
!**********************************************************************
!
!  Function: Van Leer scheme
!
!**********************************************************************

! Arguments
    real*8  :: RHS(nmet,ncell)

! Local variables
    integer :: i,j,k,adcell,adface
    real*8  :: um(nmet),up(nmet),ul(nmet),ur(nmet),flx(nmet)

    do i=1, ncell           ! Cell-cycle
        do j=1,ncellface    ! face-cycle
            adcell=cellinfo%neighbour(j,i)
            um(:) =cellinfo%cvar(:,i)
            up(:) =cellinfo%cvar(:,adcell)

            call MUSCL(ul,ur,um,up)
            call flux_vl_c(i,j,ul,ur,flx)
    !  Accumulate residual
            do k=1,nmet
                RHS(k,i)=RHS(k,i)+flx(k)
            end do
        end do  ! Face-cycle
    end do      ! Cell-cycle

    return 
end subroutine VanLeer



subroutine flux_vl_c(icell,iface,ul,ur,flux)
    use globalvar
    implicit none 
!**********************************************************************
!
! Function: compute flux accumulation
!
!**********************************************************************

! Arguments
    real*8  :: ul(nmet),ur(nmet),flux(nmet)
    integer :: icell,iface

! Local variables
    integer :: i
    real*8  :: fluxm(nmet),fluxp(nmet),norm(ndimn),enorm,uul,uur,mal,mar, &
               ccl,ccr,machm,machp,maface,ppl,ppr,ex,ey,ez

    norm(:)=cellinfo%faenorm(:,iface,icell)
    if(3 .eq. ndimn) then
        enorm=dsqrt(norm(1)**2+norm(2)**2+norm(3)**2)
        ex   =norm(1)/enorm
        ey   =norm(2)/enorm
        ez   =norm(3)/enorm
        
    ! Velocity in normal vector direction, left & right
        uul  =(ex*ul(2)+ey*ul(3)+ez*ul(4))/ul(1)
        uur  =(ex*ur(2)+ey*ur(3)+ez*ur(4))/ur(1)
    else if(2 .eq. ndimn) then
        enorm=dsqrt(norm(1)**2+norm(2)**2)
        ex   =norm(1)/enorm
        ey   =norm(2)/enorm
        
    ! Velocity in normal vector direction, left & right
        uul  =(ex*ul(2)+ey*ul(3))/ul(1)
        uur  =(ex*ur(2)+ey*ur(3))/ur(1)
    end if

    call compp(ul,ppl)
    call compp(ur,ppr)
    ccl=dsqrt(gamma*ppl/ul(1))      ! Left sound speed
    ccr=dsqrt(gamma*ppr/ur(1))      ! Right sound speed

    mal=uul/ccl     ! Left Mach number
    mar=uur/ccr     ! Right Mach number

    ! positive fluxes
    if(dabs(mal) .lt. 1.0D0) then
        call flux_vl_p(fluxp,ul,ccl,mal,norm)
    else if(mal .ge. 1.0D0) then
        call flux_vl(fluxp,ul,norm)
    else
        fluxp=0.0D0
    end if

    ! Negative fluxes
    if(dabs(mar) .lt. 1.0D0) then
        call flux_vl_m(fluxm,ur,ccr,mar,norm)
    else if(mar .ge. 1.0D0) then
        fluxm=0.0D0
    else
        call flux_vl(fluxm,ur,norm)
    end if


    do i=1,nmet
        flux(i)=fluxp(i)+fluxm(i)
    end do

    return
end subroutine flux_vl_c



subroutine flux_vl_p(fluxp,u,cc,ma,norm)
    use globalvar
    implicit none 
!**********************************************************************
!
! Function: Compute positive flux
!
!**********************************************************************

! Arguments
    real*8  :: fluxp(nmet),u(nmet),norm(ndimn),cc,ma

! Local variables
    real*8  :: enorm,exe,eye,eze,uu,pvari(nmet)
    integer :: i

    call ctop(pvari(:),u(:))
    if(3 .eq. ndimn) then
        enorm=dsqrt(norm(1)**2+norm(2)**2+norm(3)**2)
        exe  =norm(1)/enorm
        eye  =norm(2)/enorm
        eze  =norm(3)/enorm
        uu   =exe*pvari(2)+eye*pvari(3)+eze*pvari(4)
    
        fluxp(1) =0.25D0*enorm*pvari(1)*cc*(ma+1.0D0)**2
        fluxp(2) =fluxp(1)*(exe*(-uu+2.0D0*cc)/gamma+pvari(2))
        fluxp(3) =fluxp(1)*(exe*(-uu+2.0D0*cc)/gamma+pvari(3))
        fluxp(4) =fluxp(1)*(exe*(-uu+2.0D0*cc)/gamma+pvari(4))
        fluxp(5) =fluxp(1)*((-gam1*uu**2+2.0D0*gam1*uu*cc+ &
                    2.0D0*cc**2)/(gamma**2-1.0D0)+0.5D0*(pvari(2)**2+pvari(3)**2+pvari(4)**2))
    else if( 2 .eq. ndimn) then
        enorm=dsqrt(norm(1)**2+norm(2)**2)
        exe  =norm(1)/enorm
        eye  =norm(2)/enorm
        uu   =exe*pvari(2)+eye*pvari(3)
    
        fluxp(1) =0.25D0*enorm*pvari(1)*cc*(ma+1.0D0)**2
        fluxp(2) =fluxp(1)*(exe*(-uu+2.0D0*cc)/gamma+pvari(2))
        fluxp(3) =fluxp(1)*(exe*(-uu+2.0D0*cc)/gamma+pvari(3))
        fluxp(4) =fluxp(1)*((-gam1*uu**2+2.0D0*gam1*uu*cc+ &
                    2.0D0*cc**2)/(gamma**2-1.0D0)+0.5D0*(pvari(2)**2+pvari(3)**2))

    end if

    return 
end subroutine flux_vl_p


subroutine flux_vl_m(fluxm,u,cc,ma,norm)
    use globalvar
    implicit none 
!**********************************************************************
!
! Function: compute positive flux
!
!**********************************************************************

! Arguments
    real*8  :: fluxm(nmet),u(nmet),norm(ndmin),cc,ma

! Local variables
    real*8  :: enorm,exe,eye,eze,uu,pvari(nmet)
    integer :: i

    call ctop(pvari(:),u(:))
    if(3 .eq. ndimn) then
        enorm=dsqrt(norm(1)**2+norm(2)**2+norm(3)**2)
        exe  =norm(1)/enorm
        eye  =norm(2)/enorm
        eze  =norm(3)/enorm
        uu   =exe*pvari(2)+eye*pvari(3)+eze*pvari(4)
    
        fluxp(1) =-0.25D0*enorm*pvari(1)*cc*(ma-1.0D0)**2
        fluxp(2) = fluxp(1)*(exe*(-uu-2.0D0*cc)/gamma+pvari(2))
        fluxp(3) = fluxp(1)*(exe*(-uu-2.0D0*cc)/gamma+pvari(3))
        fluxp(4) = fluxp(1)*(exe*(-uu-2.0D0*cc)/gamma+pvari(4))
        fluxp(5) = fluxp(1)*((-gam1*uu**2-2.0D0*gam1*uu*cc+ &
                    2.0D0*cc**2)/(gamma**2-1.0D0)+0.5D0*(pvari(2)**2+pvari(3)**2+pvari(4)**2))
    else if(2 .eq. ndimn) then
        enorm=dsqrt(norm(1)**2+norm(2)**2)
        exe  =norm(1)/enorm
        eye  =norm(2)/enorm
        uu   =exe*pvari(2)+eye*pvari(3)
    
        fluxp(1) =-0.25D0*enorm*pvari(1)*cc*(ma-1.0D0)**2
        fluxp(2) = fluxp(1)*(exe*(-uu-2.0D0*cc)/gamma+pvari(2))
        fluxp(3) = fluxp(1)*(exe*(-uu-2.0D0*cc)/gamma+pvari(3))
        fluxp(4) = fluxp(1)*((-gam1*uu**2-2.0D0*gam1*uu*cc+ &
                    2.0D0*cc**2)/(gamma**2-1.0D0)+0.5D0*(pvari(2)**2+pvari(3)**2))
    end if

    return 
end subroutine flux_vl_m


subroutine flux_vl(flux,u,norm)
    use globalvar
    implicit none 
!**********************************************************************
!
! Function: compute flux using conservative variables
!
!**********************************************************************

! Arguments
    real*8  :: flux(nmet),u(nmet),norm(ndimn)

! Local variables
    real*8  :: ps,ex,ey,ez,enorm
    integer :: i

    call compp(u,ps)

    if(3 .eq. ndimn) then
        flux(1)=norm(1)*u(2)+norm(2)*u(3)+norm(3)*u(4)
        flux(2)=flux(1)*u(2)/u(1)+norm(1)*ps
        flux(3)=flux(1)*u(3)/u(1)+norm(2)*ps
        flux(4)=flux(1)*u(4)/u(1)+norm(3)*ps
        flux(5)=flux(1)*(u(5)+ps)/u(1)
    else if(2 .eq. ndimn) then
        flux(1)=norm(1)*u(2)+norm(2)*u(3)
        flux(2)=flux(1)*u(2)/u(1)+norm(1)*ps
        flux(3)=flux(1)*u(3)/u(1)+norm(2)*ps
        flux(4)=flux(1)*(u(4)+ps)/u(1)
    end if

    return 
end subroutine flux_vl

        
        

