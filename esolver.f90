subroutine esolver(RHS)
    use globalvar
    implicit none 
!**********************************************************************
!
! Function: compute convective flux
!
!**********************************************************************

! Arguments
    real*8  :: RHS(nmet,ncell)

    RHS = 0.0D0

    select case(spasch)
        case(1)     !Steger-Warming Scheme
            call StegerWarming(RHS)
        case(2)     ! VanLeer scheme
            call VanLeer(RHS)
        case(3)     ! Liou-Steffen
            call LiouSteffen(RHS)
    end select

    return 
end subroutine esolver
