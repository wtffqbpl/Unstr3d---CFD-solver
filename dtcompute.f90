subroutine decompute(dtmax)
    use globalvar
    implicit none
!**********************************************************************
!
! FunctionL Compute local time step
!
!**********************************************************************

! Arguments
    real*8  :: dtmax
! Local variables
    integer :: i,n
    real*8  :: norm(ndimn),enorm,ex,ey,vari(nmet),cc,deltx,delty,deltz,celvol

    dtmax=0.0D0
    if( 3 .eq. ndimn) then
        do i=1, ncell
            call ctop(vari(:),cellinfo%cvar(:,i))
            cc=dsqrt(gamma*vari(5)/vari(i))
            celvol=cellinfo%vol(i)
            deltx=0.0D0
            delty=0.0D0
            deltz=0.0D0
            do n=1,ncellface
                norm(:)=cellinfo%facenorm(:,n,i)
                deltx  =deltx+0.5D0*dabs(norm(1))
                delty  =delty+0.5D0*dabs(norm(2))
                deltz  =deltz+0.5D0*dabs(norm(3))
            end do
            deltx=(dabs(vari(2))+cc)*deltx
            delty=(dabs(vari(3))+cc)*delty
            deltz=(dabs(vari(4))+cc)*deltz
            cellinfo%eigen(1,1,i)=deltx
            cellinfo%eigen(2,1,i)=delty
            cellinfo%eigen(3,1,i)=deltz
            cellinfo%deltt(i)=cfl*celvol/(deltx+delty+deltz)
            if(cellinfo%deltt(i) .gt. dtmax) dtmax=cellinfo%deltt(i)
        end do
    else if(2 .eq. ndimn) then
        do i=1, ncell
            call ctop(vari(:),cellinfo%cvar(:,i))
            cc=dsqrt(gamma*vari(4)/vari(i))
            celvol=cellinfo%vol(i)
            deltx=0.0D0
            delty=0.0D0
            do n=1,ncellface
                norm(:)=cellinfo%facenorm(:,n,i)
                deltx  =deltx+0.5D0*dabs(norm(1))
                delty  =delty+0.5D0*dabs(norm(2))
            end do
            deltx=(dabs(vari(2))+cc)*deltx
            delty=(dabs(vari(3))+cc)*delty
            cellinfo%eigen(1,1,i)=deltx
            cellinfo%eigen(2,1,i)=delty
            cellinfo%deltt(i)=cfl*celvol/(deltx+delty)
            if(cellinfo%deltt(i) .gt. dtmax) dtmax=cellinfo%deltt(i)
        end do
    end if

    return 
end subroutine decompute
