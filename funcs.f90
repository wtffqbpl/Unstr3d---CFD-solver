subroutine ctop(up,uc)
    use globalvar
    implicit none 
!***********************************************************************
!
!  Function: Convert conservative variables to primitive variables
!
!**********************************************************************
!
! Arguments
    real*8   :: up(nmet),uc(nmet)

! Local variables
    integer  :: i
    real*8   :: temp1,temp2,temp3,temp4,temp5
    
    ! Converting
    temp1 = uc(1)
    temp2 = uc(2)
    temp3 = uc(3)
    up(1) = temp1
    up(2) = temp2/temp1
    up(3) = temp3/temp1
    if(3 .eq. ndimn) then
        temp4 = uc(4)
        temp5 = uc(5)
        uc(4) = temp4/temp1
        uc(5) = gam1*(temp5-0.50D0*(temp2*temp2+temp3*temp3+temp4*temp4)/temp1)
    else if( 2 .eq. ndimn) then
        temp4 = uc(4)
        uc(4) = gam1*(temp4-0.50D0*(temp2*temp2+temp3*temp3)/temp1)
    end if

    return 
end subroutine ctop


subroutine ptoc(up,uc)
    use globalvar
    implicit none 
!***********************************************************************
!
! Function: Convert primitive varibles to conservative variables
!
!**********************************************************************
! Arguments
    real*8  :: up(nmet),uc(nmet)


    ! Converting
    uc(1) = up(1)
    uc(2) = up(1)*up(2)
    uc(3) = up(1)*up(3)
    if(3 .eq. ndimn) then
        uc(4) = up(1)*up(4)
        uc(5) = up(5)/gam1+0.50D0*up(1)*(up(2)*up(2)+up(3)*up(3)+up(4)*up(4))
    else if(2 .eq. ndimn) then
        uc(4) = up(4)/gam1+0.50D0*up(1)*(up(2)*up(2)+up(3)*up(3))
    end if

    return
end subroutine ptoc



subroutine compp(u,ps)
    use globalvar
    implicit none 
!**********************************************************************
!
! Function: Compute pressure term using convervative variables
!
!**********************************************************************

! Arguments
    real*8  :: u(nmet),ps
    
    if(3 .eq. ndimn) then
        ps=gam1*(u(5)-0.5D0*(u(2)**2+u(3)**2+u(4)**2)/u(1))
    else if(2 .eq. ndimn) then
        ps=gam1*(u(5)-0.5D0*(u(2)**2+u(3)**2)/u(1))
    end if
    ps=dmax1(ps,1.0D-6)

    return 
end subroutine compp



subroutine progBar(icell,ncell)
    implicit none 
!**********************************************************************
!
! Function: Progress bar
!
!**********************************************************************
!
! Arguments
    integer :: icell,ncell
! Local variables
    integer :: i
    real    :: v,vtot,prog

    v = real(icell)
    vtot = real(ncell)
    prog=v/vtot*100D0
    write(*,'(7a,f6.1,a,$)') (char(8),i=1,7), prog,'%'

    return 
end subroutine progBar


subroutine outfield(iunit,fieldname,filename)
    use globalvar
    implicit none 
!**********************************************************************
!
! Fucntion: Progress bar
!
!**********************************************************************
!
! Arguments
    integer   :: iunit
    character :: fieldname*80,filename*20

    ! Local variables
    integer :: i,j,nodest,nodend,noderd,nodfrd
    real*8  :: cvari(nmet),pvari(nmet)

    open(iunit,file=trim(filename))
    do i=1,ncell
        cvari(:)=cellinfo%cvar(:,i)
        call ctop(pvari(:),cvari(:))
        cellinfo&pvar(:,i)=pvari(:)
    end do

    if(3 .eq. ndimn) then
        write(iunit,1301)trim(fieldname),trim(fieldname),npoin,ncell
        write(iunit,*)(coords(1,i),i=1,npoin)
        write(iunit,*)(coords(2,i),i=1,npoin)
        write(iunit,*)(coords(3,i),i=1,npoin)
    else if(2 .eq. ndimn) then
        write(iunit,1300)trim(fieldname),trim(fieldname),npoin,ncell
        write(iunit,*)(coords(1,i),i=1,npoin)
        write(iunit,*)(coords(2,i),i=1,npoin)
    end if

    do j=1,nmet
        write(iunit,*)(cellinfo%pvar(j,i),i=1,ncell)
    end do
    if(3 .eq. ndimn) then
        do i=1,ncell
            nodest=cellinfo%vertex(1,i)
            nodend=cellinfo%vertex(2,i)
            noderd=cellinfo%vertex(3,i)
            nodfrd=cellinfo%vertex(4,i)
            write(iunit,'(4(2X,I6))') nodest,nodend,noderd,nodfrd
        end do
    else if( 2 .eq. ndmin) then
        do i=1,ncell
            nodest=cellinfo%vertex(1,i)
            nodend=cellinfo%vertex(2,i)
            noderd=cellinfo%vertex(3,i)
            write(iunit,'(4(2X,I6))') nodest,nodend,noderd
        end do
    end if
    close(iunit)
    
!******************* FORMATS ******************************************
1300 format('TITLE ="',A,'"',/,&
            'VARIABLES="X", "Y", "Rho", "U", "V","P"',/, &
            'ZONE T="',A,'", DATAPACKING=BLOCK, NODES=', I6,', ELEMENTS=',I6, &
            ', ZONETYPE=FETRIANGLE, VARLOCATION=([3-6]=CELLCENTERED)')
1301 format('TITLE ="',A,'"',/,&
            'VARIABLES="X", "Y", "Z",  "Rho", "U", "V", "W", "P"',/, &
            'ZONE T="',A,'", DATAPACKING=BLOCK, NODES=', I6,', ELEMENTS=',I6, &
            ', ZONETYPE=FETETRAHEDRON, VARLOCATION=([4-8]=CELLCENTERED)')

    return 
end subroutine outfield
        

