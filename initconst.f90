subroutine initconst
    use globalvar
    implicit none 
!*********************************************************************************
!
!  Function: read mesh
!
!********************************************************************************
!
! Local variables


    pi    = 4.0D0*atan(1.0D0)
    rad   = 180.0D0 / pi
    alpha = alpha / rad

    gam1  = gamma - 1.0D0

    rhoinf = pinf / (rgas * tinf) 
    qinf   = machinf*sqrt(gamma*rgas*tinf)
    if( 3 .eq. ndimn) then
        uinf = qinf*cos(alpha)*cos(beta)
        winf = qinf*sin(alpha)*cos(beta)
        vinf = qinf*sin(beta)
    else if( 2 .eq. ndimn) then
        uinf = qinf*cos(alpha)
        vinf = qinf*sin(alpha)
    else
        write(*,*) 'DIMENSION is wrong!!!'
    end if

end subroutine initconst
