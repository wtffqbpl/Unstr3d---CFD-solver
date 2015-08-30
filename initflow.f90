subroutine initflow
    use globalvar
    implicit none 
!*********************************************************************************
!
!   Function: read mesh and boundary conditions
!
!********************************************************************************

! Local variables
    integer     :: i,j,adcell,adflag,initout,nodest,nodend,noderd,iunit,lenstr, &
                   lenend,lentot,ttnodes,ttcells,iunit2,itemp
    character   :: fieldname*80,filename*20,tempchar*150
    real*8      :: temp, temp1(5)
    logical     :: alive


    iunit2   = 2501
    iunit    = 2500
    laststep = 0

    if( 0 .eq. inpres ) then
        do i=1, ncell
            cellinfo%cvar(1,i) = rhoinf
            cellinfo%cvar(2,i) = rhoinf*uinf
            cellinfo%cvar(3,i) = rhoinf*vinf
            if(3 .eq. ndimn) then
                cellinfo%cvar(4,i) = rhoinf*winf
                cellinfo%cvar(5,i) = pinf/gam1+0.5D0*rhoinf*qinf*qinf
            else if( 2 .eq. ndimn ) then
                cellinfo%cvar(4,i) = pinf/gam1+0.5D0*rhoinf*qinf*qinf
            else
                write(*,*) 'DIMENSION is wrong!!!'
                stop
            end if
        end do
    else if( 1 .eq. inpres ) then
        inquire(file=restartfile, exist=alive)
        if(.not. alive) then
            write(*,*) 'There is no restart file, please make sure restart file or filename'
            stop
        end if

        open(iunit,file=trim(restartfile))
        read(iunit,*)
        read(iunit,*)
        read(iunit,*)
        read(iunit,*) (temp, i1, npoin)
        read(iunit,*) (temp, i1, npoin)
        if(3 .eq. ndimn) then
            read(iunit,*) (temp, i1, npoin)
        end if
        do j=1,nmet
            read(iunit,*)(cellinfo%pvar(j,i), i=1,ncell)
        end do
        do i=1,ncell
            call ptoc(cellinfo%pvar(:,i),cellinfo%cvar(:,i))
        end do
        close(iunit)

        inquire(file'convergence.plt', exist=alive)
        if(.not. alive) then
            write(*,*) 'There is no convergence.plt file'
            stop
        end if
        open(iunit,file='convergence.plt')
        open(iunit2,file='convergence1.plt')
        read(iunit,'(A)') tempchar
        write(iunit2,'(A)') trim(tempchar)
        read(iunit,'(A)') tempchar
        write(iunit2,'(A)') trim(tempchar)
        if( 3 .eq. ndimn) then
            do while( .true. )
                read(iunit,*,iostat=itemp)laststep,temp1(1),temp1(2),temp1(3),temp1(4),temp1(5)
                if(itemp .ne. 0) exit
                write(iunit2,2551)laststep,temp1(1),temp1(2),temp1(3),temp1(4),temp1(5)
            end do
        else if(2 .eq. ndimn) then
            do while( .true. )
                read(iunit,*,iostat=itemp)laststep,temp1(1),temp1(2),temp1(3),temp1(4)
                if(itemp .ne. 0) exit
                write(iunit2,2551)laststep,temp1(1),temp1(2),temp1(3),temp1(4)
            end do
        end if
        close(iunit, status='delete')
        rewind(iunit2)
        open(iunit,file='convergence.plt')
        read(iunit,'(A)') tempchar
        write(iunit2,'(A)') trim(tempchar)
        read(iunit,'(A)') tempchar
        write(iunit2,'(A)') trim(tempchar)
        if(3 .eq. ndimn) then
            do i=1,laststep
                read(iunit2,*)itemp,temp1(1),temp1(2),temp1(3),temp1(4),temp1(5)
                write(iunit,2551)itemp,temp1(1),temp1(2),temp1(3),temp1(4),temp1(5)
            end do
        else if( 2 .eq. ndimn) then
            do i=1,laststep
                read(iunit2,*)itemp,temp1(1),temp1(2),temp1(3),temp1(4)
                write(iunit,2550)itemp,temp1(1),temp1(2),temp1(3),temp1(4)
            end do
        end if
        close(iunit2,status='delete')
        close(iunit)
    end if

    ! Update boundary conditions
    do i=ncell+1,nvitgho
        adcell=cellinfo%neighbour(1,i)
        do j=1,ncellface
            adflag = cellinfo%nflag(j,adcell)
            call bcond(i,adcell,adflag,j)
        end do
    end do

    ! Output initial flow field
    initout=1200
    fieldname='initial flow field'
    filename = 'initialflow.plt'
    call outfield(initout,fieldname,filename)

    !**********************************************************************
2500    format(1X,I5,F5.2,5X,F5,2,5X,F5.2,5X,F5.2)
2501    format(1X,I5,F5.2,5X,F5,2,5X,F5.2,5X,F5.2,5X,F5.2)

    return
end subroutine initflow
