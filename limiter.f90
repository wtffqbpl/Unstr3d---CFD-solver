subroutine limiter(ul,ur,um,up,gradi,gradj,icell,jcell,iface)
    use globalvar
    implicit none 
!*************************************************************************************
!
! Function: interpolate variables of UL and UR
!   on unstructured grids, the purpose of a limiter is to reduce the gradient
!   used to reconstruct the left and right state at the face of the control volume.
!   The limiter function must be zero at strong discontinuities, in order to obtain
!   a first-order upwind scheme which guarantees monotonicity.
!*************************************************************************************
! Arguments:
    real*8  :: ul(nmet),ur(nmet),um(nmet),up(nmet),gradi(ndimn,nmet),gradj(ndimn,nmet)
    integer :: icell,jcell,iface

! Local variables
    integer :: i
    real*8  :: temp1,temp2

    do i=1,nmet
        call Barth_Jespersen(temp1,gradi(:,i),icell,jcell,iface)
        ul(i)=um(i)
    end do

    return 
end subroutine limiter


subroutine Barth_Jespersen(resout,grad,icell,jcell,iface)
    use globalvar
    implicit none
!*************************************************************************************
!
! Function: Bath and Jespersen limiter
!
!*************************************************************************************
!
! Arguments
    real*8  :: resout,grad(2)
    integer :: icell,jcell,iface

! Local variables
    real*8  :: norm(ndimn),delta2

    resout = 0.0D0

    norm(:) = cellinfo%cetomid(:,iface,icell)
    delta2  = grad(1)*norm(1)+grad(2)*norm(2)

    return
end subroutine Barth_Jespersen



subroutine MUSCL(ul,ur,um,up)
    use globalvar
    implicit none 
!*************************************************************************************
!
! Function: Function: MUSCL limiter
!
!*************************************************************************************

! Arguments:
    real*8  :: ul(nmet),ur(nmet),um(nmet),up(nmet)

    ul = um
    ur = up

    return
end subroutine MUSCL
