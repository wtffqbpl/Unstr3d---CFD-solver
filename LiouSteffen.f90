subroutine LiouSteffen(RHS)
    use globalvar
    implicit none 
!*************************************************************************************
!
! Function: Liou-Steffen scheme
!
!*************************************************************************************

! Arguments
    real*8  :: RHS(nmet,ncell)
! Local variables
    integer :: i,j,k,adcell,adface
    real*8  :: um(nmet),up(nmet),ul(nmet),ur(nmet),flx(nmet),gradi(ndimn,nmet), &
               gradj(ndimn,nmet)


    do i=1,ncell    ! cell-type

!   face-cycle (we can use local number of cell faces in mixed grids)
        do j=1,ncellface
            adcell = cellinfo%neighbour(j,i)
            um(:) = cellinfo%neighbour(:,i)
            um(:) = cellinfo%cvar(:,adcell)

    ! call evalgrad(adcell,gradj)
            call MUSCL(ul,ur,um,up)
    ! call limiter(ul,ur,um,up,gradi,gradj,i,adcell,j)
            call flux_ls_c(i,j,ul,ur,flx)

    ! accumulate residual
            do k=1,nmet
                RHS(k,i) = RHS(k,i)+flx(k)
             end do
        end do
    end do

    return
end subroutine LiouSteffen


subroutine flux_ls_c(icell,iface,ul,ur,flux)
    use globalvar
    implicit none 
!*************************************************************************************
!
! Function: compute flux accumulation
!
!*************************************************************************************
!
! Arguments
    real*8  :: ul(nmet),ur(nmet),flux(nmet)
    integer :: icell,iface

! Local variables
    integer :: i
    real*8  :: fluxm(nmet),fluxp(nmet),norm(ndimn),enorm,ex,ey,ez,uul, &
               uur,mal,mar,ccl,ccr,machm,machp,maface,ppl,ppr

    norm(:) = cellinfo%facenorm(:,iface,icell)

    if(3. eq. ndimn) then
        enorm = dsqrt(nor(1)**2+norm(2)**2+norm(3)**2)
        ex = norm(1)/enorm
        ey = norm(2)/enorm
        ez = norm(3)/enorm
        uul=(ul(2)*ex+ul(3)*ey+ul(4)*ez)/ul(1)
        uur=(ur(2)*ex+ur(3)*ey+ur(4)*ez)/ur(1)
    else if(2 .eq. ndimn) then
        enorm = dsqrt(norm(1)**2+norm(2)**2)
        ex = norm(1)/enorm
        ey = norm(2)/enorm
        uul=(ul(2)*ex+ul(3)*ey)/ul(1)
        uur=(ur(2)*ex+ur(3)*ey)/ur(1)
    end if

    call compp(ul,ppl)
    call compp(ur,ppr)

    ccl = sqrt(gamma*ppl/ul(1))     ! Left sound speed
    ccr = sqrt(gamma*ppr/ur(1))     ! Right sound speed

    mal = uul/ccl   ! Left Mach number
    mar = uur/ccr   ! Right Mach number

    if(dabs(mal) .le. 1.0D0) then
        machp = 0.25D0*(mal+1.0D0)**2
    else
        machp = 0.50D0*(mal+dabs(mal))
    end if

    if(dabs(mar) .le. 1.0D0) then
        machm =-0.25D0*(mar-1.0D0)**2
    else
        machm = 0.50D0*(mar-dabs(mar)
    end if

    maface = machp+machm


! Positive fluxes
    call flux_ls_p(fluxp,ul,ccl,maface,mal,norm)

! Negative fluxes
    call flux_ls_m(fluxp,ur,ccr,maface,mar,norm)

    do i=1,nmet
        flux(i) = fluxp(i)+fluxm(i)
    end do

    return 
end subroutine flux_ls_c


subroutine flux_ls_p(fluxp,ul,ccl,maface,mal,norm)
    use globalvar
    implicit none 
!*************************************************************************************
!
! Function: Compute positive flux
!
!*************************************************************************************

! Arguments:
    real*8  :: fluxp(nmet),ul(nmet),norm(ndimn),maface,mal,ccl

! Local variables
    real*8  :: enorm,exe,eye,eze,sm1,sm2,pp

    if(3 .eq. ndimn) then
        enorm = dsqrt(norm(1)**2+norm(2)**2+norm(3)**2)
        exe   = norm(1)/enorm
        eye   = norm(2)/enorm
        eze   = norm(3)/enorm
        sm1   = 0.50D0*(maface+dabs(maface))

        if(dabs(mal) .le. 1.0D0) then
            sm2 = 0.25D0*(1.0D0+mal)**2*(2.0D0-mal)
        else
            sm2 = 0.50D0*(mal+dabs(mal))/mal
        end if

        call compp(ul,pp)
        fluxp(1) = enorm*sm1*ul(1)*ccl
        fluxp(2) = fluxp(1)*ul(2)/ul(1)+enorm*exe*pp*sm2
        fluxp(3) = fluxp(1)*ul(3)/ul(1)+enorm*eye*pp*sm2
        fluxp(4) = fluxp(1)*ul(4)/ul(1)+enorm*eze*pp*sm2
        fluxp(5) = enorm*sm1*(ul(5)+pp)*cc1
    else if(2 .eq. ndimn) then
        enorm = dsqrt(norm(1)**2+norm(2)**2)
        exe   = norm(1)/enorm
        eye   = norm(2)/enorm
        sm1   = 0.50D0*(maface+dabs(maface))

        if(dabs(mal) .le. 1.0D0) then
            sm2 = 0.25D0*(1.0D0+mal)**2*(2.0D0-mal)
        else
            sm2 = 0.50D0*(mal+dabs(mal))/mal
        end if

        call compp(ul,pp)
        fluxp(1) = enorm*sm1*ul(1)*ccl
        fluxp(2) = fluxp(1)*ul(2)/ul(1)+enorm*exe*pp*sm2
        fluxp(3) = fluxp(1)*ul(3)/ul(1)+enorm*eye*pp*sm2
        fluxp(4) = enorm*sm1*(ul(4)+pp)*cc1
    end if

    return
end subroutine flux_ls_p


subroutine flux_ls_m(fluxm,ur,ccr,maface,mar,norm)
    use globalvar
    implicit none 
!*************************************************************************************
!
! Function: Compute positive flux
!
!*************************************************************************************

! Arguments
    real*8  :: fluxp(nmet),ur(nmet),norm(ndimn),maface,mar,ccr

! Local variables
    real*8  :: enorm,exe,eye,eze,sm1,sm2,pp

    if(3 .eq. ndimn) then
        enorm = dsqrt(norm(1)**2+norm(2)**2+norm(3)**2)
        exe   = norm(1)/enorm
        eye   = norm(2)/enorm
        eze   = norm(3)/enorm
        sm1   = 0.50D0*(maface-dabs(maface))

        if(dabs(mar) .le. 1.0D0) then
            sm2 = 0.25D0*(1.0D0-mar)**2*(2.0D0+mar)
        else
            sm2 = 0.50D0*(mar-dabs(mar))/mar
        end if

        call compp(ur,pp)
        fluxp(1) = enorm*sm1*ur(1)*ccr
        fluxp(2) = fluxp(1)*ur(2)/ul(1)+enorm*exe*pp*sm2
        fluxp(3) = fluxp(1)*ur(3)/ul(1)+enorm*eye*pp*sm2
        fluxp(4) = fluxp(1)*ur(4)/ul(1)+enorm*eze*pp*sm2
        fluxp(5) = enorm*sm1*(ur(5)+pp)*cc1
    else if(2 .eq. ndimn) then
        enorm = dsqrt(norm(1)**2+norm(2)**2)
        exe   = norm(1)/enorm
        eye   = norm(2)/enorm
        sm1   = 0.50D0*(maface-dabs(maface))

        if(dabs(mar) .le. 1.0D0) then
            sm2 = 0.25D0*(1.0D0-mar)**2*(2.0D0+mar)
        else
            sm2 = 0.50D0*(mar-dabs(mar))/mar
        end if

        call compp(ur,pp)
        fluxp(1) = enorm*sm1*ur(1)*ccr
        fluxp(2) = fluxp(1)*ur(2)/ur(1)+enorm*exe*pp*sm2
        fluxp(3) = fluxp(1)*ur(3)/ur(1)+enorm*eye*pp*sm2
        fluxp(4) = enorm*sm1*(ur(4)+pp)*cc1
    end if

    return
end subroutine flux_ls_m




