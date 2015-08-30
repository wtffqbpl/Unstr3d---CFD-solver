subroutine readparam
    use globalvar
    implicit none 
!*******************************************************************************
!
!  Function: read input parameters
!
!*******************************************************************************
!
! Local variables
    integer :: i, iunit

    iunit = 1000

    open(iunit, file'input.dat')
    read(iunit,*)
    read(iunit,*) infile
    read(iunit,*) ndimn
    write(*,*) ' The mesh file is ',trim(infile)
    !Physics =external flow
    read(iunit,*)
    read(iunit,*)
    read(iunit,*)
    read(iunit,*) machinf
    read(iunit,*) alpha
    read(iunit,*) beta
    read(iunit,*) pinf
    read(iunit,*) tinf
    read(iunit,*) rgas
    read(iunit,*) gamma
    read(iunit,*) keyvis

    ! Geometrical reference values
    read(iunit,*) 
    read(iunit,*)
    read(iunit,*)
    read(iunit,*) xref
    read(iunit,*) yref
    read(iunit,*) zref
    read(iunit,*) cref
    read(iunit,*) aref

    ! Iteration control
    read(iunit,*)
    read(iunit,*)
    read(iunit,*)
    read(iunit,*) maxiter
    read(iunit,*) ninterval
    read(iunit,*) convtol
    read(iunit,*) maxtol

    ! Numberical paramters
    read(iunit,*)
    read(iunit,*)
    read(iunit,*)
    read(iunit,*) cfl
    read(iunit,*) timeid
    read(iunit,*) vis2
    read(iunit,*) vis4
    vis4 = 1.0D0 / vis4
    read(iunit,*) nrk
    allocate(ark(nrk))
    read(iunit,*) (ark(i), i=1,nrk)


    ! Spatial scheme
    read(iunit,*)
    read(iunit,*)
    read(iunit,*)
    read(iunit,*) inpres
    if( 1 .eq. inpres) then
        read(iunit,*) restartfile
    else
        read(iunit,*)
    end if

    close(iunit)

end subroutine readparam
