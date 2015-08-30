program unstr3d
    use globalvar
    use ifport
    implicit none   
!***********************************************************************************
!
!   An unstructured solver for 2D & 3D Euler equation (cell-centred scheme)
!                   Created by Ren Yuanjun
!                     (Version: V1.0)
!
!***********************************************************************************
!
! I/O channels;

! icons  = convergence history (Tecplot format )
! iflows = flow field variables

! Local variables
    integer    :: iter, iflag, icons, iflows
    character  :: filename*20, fieldname*80, concfile*20
    real       :: time_begin, time_end, start_clock, end_clock
    logical    :: alive

    iflag     = 0
    iflows    = 2100
    icons     = 2101
    concfile  = 'convergence.plt'
    filename  = 'flowfield.plt'
    fieldname = 'ddd'

    write(*,1000)

    !read input parameters
    write(*,*) "    Reading input paramters ....... "
    call Readparam

    ! Initialize some constants 
    call initconst

    ! Read mesh
    write(*,*) "    Reading mesh information....... "

    ! Calculate metrics, cell volumes
    write(*,*) "    Guessing initial solution ....... "
    call initflow

    ! Open file for convergence history
    if( 1 .eq. inpres ) then
        inquire( file=concfile, exist = alive )
        if( .not. alive ) then
            write(*,*) "There is no previous convergence file"
            stop
        end if
        open(icons, file=concfile, status='OLD', position='APPEND')
    else
        open(icons, file=concfile)
        if( 3 .eq. ndimn ) then
            write(icons, 1035) fieldname
        else if( 2. eq. ndimn ) then
            write(icons, 1055 ) fieldname
        else
            write(*,*) "DIMENSION is wrong!!!"
        end if
    end if
    
    write(*, 1036)

    !----------------------------START------------------------------------
    ! Iterating until convergence criterion is statisfied

    do iter=1,maxiter
        call CPU_TIME(time_begin)
        start_clock = dclock()
        call RKmarching 

        if( 0 .eq. mod(iter, niterval)) then
            call outfield(iflows, fieldname, filename)
        end if

        if(0 .eq. mod(iter,20)) then
            write(*,*)
            write(*,1036)
        end if

        ! Check convergence vertia
        call conver(iflag, time_begin, start_clock, icons, iter)
        if(iflag .eq. 1) then   ! Satisfies convergent certia
            call outfield(iflows, fieldname, filename)
            goto 2000
        end if
    end do


    ! output the results
2000 write(*, 1045)

    write(*,*) "    Finished! "


    !******************************FORMATS********************************
    1000 format(//,11X,59('*'),/,11X,'*',57X,'*',/,11X,'* ',&
                'SOLUTION OF 2D & 3D EULER EQUATIONS ON UNSTRUCTURED GRIDS  *', &
                /,11X,'*',/,11X,59('*'),//)
    1035 format('TITLE="',A,'"',/,&
                'VARIABLES= step, Rho, Rho*u, Rho*v, Rho*w, Rho*E')
    1055 format('TITLE="',A,'"',/, &
                'VARIABLES= step, Rho, Rho*u, Rho*v, Rho*E')
    1036 format(2X,'STEP', 5X,'RESID',7X,'RESMAX',4X,'INDEX',4X,'CPU_TIME',8X,'DCLOCK')
    1045 format(/,80('-'),//,' Writing plot files ...')


    stop
end program unstr3d
