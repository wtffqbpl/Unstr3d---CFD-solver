subroutine RKmarching
    use globalvar
    implicit none 
!**********************************************************************
!
!  Function: Runge-Kutta (explicit multistage scheme)
!
!**********************************************************************

! Local variables
    real*8, allocatable :: RHS(:,:),RHSV(:,:),cvari0(:,:)
    real*8              :: cellvol,ddt
    integer             :: irk,i,n,adcell,adflag

    allocate(RHS(nmet,ncell),RHSV(nmet,ncell),cvari0(nmet,ncell))


!******************************
    RHSV(:,:) = 0.0D0
    do i=1,ncell
        do n=1,nmet
            cvari0(n,i)=cellinfo%cvar(n,i)
        end do
    end do
    cellinfo%cvar0(:,:)=cvari0(:,:)
!******************************

! Compute local timestep
    call dtcompute(ddt)

! Runge-Kutta iteration
    do irk=1,nrk

    ! Calculate convective flux
        call esolver(RHS)
    ! Calculate viscous flux
        if(keyvis .eq. 0) RHSV=0.0D0
        if(keyvis .eq. 1) call viscoslam(RHSV)
    !   if(keyvis .eq. 2) call viscousturb(RHSV)

    ! Compute convervative variables of i-stage
        do i=1,ncell
            cellvol = cellinfo%vol(i)
            if(timeid .eq. 1) ddt = cellinfo%deltt(i)
            do n=1,nmet
                cellinfo%cvar(n,i)=cvari0(n,i)-ark(irk)*ddt*(RHS(n,i)-RHSV(n,i))/cellvol
            end do
        end do

    ! Boundary conditions
        do i=ncell+1,nvitgho
            adcell=cellinfo%neighbour(1,i)
            do n=1,ncellface
                adflag=cellinfo%nflag(n,adcell)
                call bcond(i,adcell,adflag,n)
            end do
        end do

    end do

    deallocate(RHS,RHSV,cvari0)

    return
end subroutine RKmarching
