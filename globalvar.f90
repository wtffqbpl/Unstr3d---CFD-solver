module globalvar
    implicit none 
!******************************************
!
!!
!
!*****************************
!
    type :: cell_type     !cell type
        integer, allocatable :: vertex(:,:), neighbour(:,:), nflag(:,:)
        real*8,  allocatable :: vol(:), facenorm(:,:,:), cvar(:,:), cvar0(:,:), pvar(:,:), &
                                centroid(:,:), cetomid(:,:,:), eigen(:,:,:),deltt(:)
    end type cell_type

    character                :: infile*80   ! input file name

    !parameters of geometry
    real*8, allocatable :: coords(:,:)
    type(cell_type)     :: cellinfo
    integer             :: ncell, nbsurf, npoin, ndimn, nlines,facetet(3,4), facetri(2,3), &
                           nmet, nvitgho, ncellface
    real*8              :: xref, yref, zref, cref, aref


    ! Temporal scheme
    real*8, allocatable :: ark(:)
    integer             :: nrk, maxiter, timeid, ninterval
    real*8              :: pi, rad, cfl, convtol, maxtol

    ! Parameters of physical conditions of the problem
    real*8              :: gamma,rgas, gam1,alpha,beta,qinf,pinf, tinf, machinf, rhoinf, &
                           uinf, vinf, winf
    integer             :: keyvis

    ! Spatial scheme
    integer             :: spasch
    real*8              :: vis2, vis4

    ! Initialization paramters
    integer             :: inpres, laststep
    character           :: restartfile*80


end module globalvar
