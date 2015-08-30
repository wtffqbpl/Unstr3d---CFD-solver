subroutine evalgrad(gradu)
    use globalvar
    implicit none 
!********************************************************************************
!
!  Function: the gradient at some cell-centroid I using Green-Gauss method
!
!********************************************************************************

! Arguments
    real*8   :: gradu(ndimn,nmet,ncell)
    integer  :: icell

! Local variables
    integer  :: i,j,adcell
    real*8   :: cui(nmet),cuj(nmet),norm(2),volumn

    gradu=0.0D0

! conservative variables of ICELL
    cui(:) = cellinfo%neighbour(j,icell)      
    volumn = cellinfo%vol(icell)

    do j=1,3    
        adcell = cellinfo%neighbour(j,icell)    ! number of neighbour cell
        norm(:) = cellinfo%facenorm(:,j,icell)  ! face norms

!   conservative variables of neighbour cell
        cuj(:) = cellinfo%cvar(:,adcell)

!   compute gradient of very conservative variables
        do i=1,nmet
            gradu(1,i,icell) = gradu(1,i,icell)+0.5D0*(cui(i)+cui(j))*norm(1)/volumn
            gradu(2,i,icell) = gradu(2,i,icell)+0.5D0*(cui(i)+cui(j))*norm(2)/volumn
        end do
    end do

    return 
ens subroutine evalgrad
