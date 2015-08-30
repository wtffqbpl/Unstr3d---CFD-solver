subroutine conver(iflag,time_begin,start_clock,icons,iter)
    use globalvar
    use ifport
    implicit none 
!*************************************************************************************
!
! Function: check convergence
!
!*************************************************************************************

! Arguments:
    integer :: iflag,iter,icons
    real    :: time_begin,start_clock

! Local variables
    real*8  :: drho,drmax,dr,dr1,drhou,drhov,drhow,drhoe
    real    :: time_end,end_clock
    integer :: i,icellmax

! Compute the residual
    drho = 0.0D0
    drhou= 0.0D0
    drhov= 0.0D0
    drhow= 0.0D0
    drhoe= 0.0D0
    drmax= 0.0D0


    if(3 .eq. ndimn) then
        do i=1,ncell
            dr      = cellinfo%cvar(1,i)-cellinfo%cvar0(1,i)
            drho    =drho+dr*dr
            dr1     = cellinfo%cvar(2,i)-cellinfo%cvar0(2,i)
            drhou   =drhou+dr1*dr1
            dr1     = cellinfo%cvar(3,i)-cellinfo%cvar0(3,i)
            drhov   =drhov+dr1*dr1
            dr1     = cellinfo%cvar(4,i)-cellinfo%cvar0(4,i)
            drhow   =drhow+dr1*dr1
            dr1     = cellinfo%cvar(5,i)-cellinfo%cvar0(5,i)
            drhoe   =drhoe+dr1*dr1

            if(abs(dr) .ge. drmax) then
                drmax=abs(dr)
                icellmax=i
            end if
        end do
        drho  =log10(dabs(drho ))
        drhou =log10(dabs(drhou))
        drhov =log10(dabs(drhov))
        drhow =log10(dabs(drhow))
        drhoe =log10(dabs(drhoe))
        drmax =log10(dabs(drmax))
        call CPU_TIME(time_end)
        end_clock=dclock()

        write(icons,1038)iter+laststep,drho,drhou,drhov,drhow,drhoe
    else if(2 .eq. ndimn) then
        do i=1,ncell
            dr      = cellinfo%cvar(1,i)-cellinfo%cvar0(1,i)
            drho    =drho+dr*dr
            dr1     = cellinfo%cvar(2,i)-cellinfo%cvar0(2,i)
            drhou   =drhou+dr1*dr1
            dr1     = cellinfo%cvar(3,i)-cellinfo%cvar0(3,i)
            drhov   =drhov+dr1*dr1
            dr1     = cellinfo%cvar(4,i)-cellinfo%cvar0(4,i)
            drhoe   =drhoe+dr1*dr1

            if(abs(dr) .ge. drmax) then
                drmax=abs(dr)
                icellmax=i
            end if
        end do
        drho  =log10(dabs(drho ))
        drhou =log10(dabs(drhou))
        drhov =log10(dabs(drhov))
        drhoe =log10(dabs(drhoe))
        drmax =log10(dabs(drmax))
        call CPU_TIME(time_end)
        end_clock=dclock()

        write(icons,1039)iter+laststep,drho,drhou,drhov,drhoe
    end if

    write(*,1037)iter+laststep,drho,drmax,icellmax,time_end-time_begin, &
                 end_clock-start_clock

    if((drho .lt. convtol) .and. (drmax .lt. maxtol)) iflag = 1


!*********************************** FORMAT ******************************************
1037 format(1X,I5,F6.2,7X,F6.2,4X,I5,4X,F7.3,8X,F7.3)
1038 format(1X,I5,F6.2,5X,F6.2,5X,F6.2,5X,F6.2,5X,F6.2)
1039 format(1X,I5,F6.2,5X,F6.2,5X,F6.2,5X,F6.2)


    return 
end subroutine conver
