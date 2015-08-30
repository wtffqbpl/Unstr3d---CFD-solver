subroutine readmsh
    use globalvar
    implicit none
!**************************************************************************************
!
! Function: read mesh and boundary conditions
!
!*************************************************************************************

! Local variables
    integer, allocatable    :: ncellm(:),nbsurfm(:)
    integer                 :: i,j,k,l,ntemp1,ntemp2,ntemp3,ioput,nzone
    logical                 :: mask(4)

    ioput = 1200
    open(ioput,file=trim(infile))

    read(ioput,*)
    read(ioput,*)
    read(ioput,*)
    read(ioput,*)
    read(ioput,*) ntemp1,nzone,npoin,nvp
    if(ntemp1 .ne. ndimn) then
        write(*,*) "The mesh dimension is different with input file, please check!!!"
        stop
    end if
    read(ioput,*)
    allocate(ncellm(nzone),nbsurfm(nzone))
    do i=1,nzone
        read(ioput,*) ntemp1,ncellm(i),nbsurfm(i)
    end do

! initialize some matries
    nmet = ndimn+2
    ncell = 0
    nbsurf = 0
    do i=1,nzone
        ncell = ncell+ncellm(i)
        nbsurf = nbsurf+nbsurfm(i)
    end do
    deallocate(ncellm,nbsurfm)
    nvitgho = ncell+nbsurf
    if(3 .eq. ndimn) then
        allocate(cellinfo%vertex(4,ncell)           ,cellinfo%neighbour(4,nvitgho),     &
                 cellinfo%nflag(4,nell)             ,cellinfo%facenorm(3,4,ncell),      &
                 cellinfo%vol(nvitgho)              ,cellinfo%cvar(nmet,nvitgho),       &
                 cellinfo%cvar0(nmet,ncell)         ,cellinfo%centroid(ndimn,ncell),    &
                 cellinfo%cetomid(2,3,ncell)        ,cellinfo%deltt(ncell),             &
                 cellinfo%eigen(ndimn,2,ncell)      ,cellinfo%pvar(nmet,ncell),         &
                 coords(ndimn,npoin))
        facetet(:,1) = [1,2,3]
        facetet(:,2) = [2,4,3]
        facetet(:,3) = [3,4,1]
        facetet(:,4) = [4,2,1]
        ncellface = 4
    else if(2 .eq. ndimn) then
        allocate(cellinfo%vertex(3,ncell)           ,cellinfo%neighbour(3,nvitgho),     &
                 cellinfo%nflag(3,nell)             ,cellinfo%facenorm(2,3,ncell),      &
                 cellinfo%vol(nvitgho)              ,cellinfo%cvar(nmet,nvitgho),       &
                 cellinfo%cvar0(nmet,ncell)         ,cellinfo%centroid(ndimn,ncell),    &
                 cellinfo%cetomid(2,3,ncell)        ,cellinfo%deltt(ncell),             &
                 cellinfo%eigen(ndimn,2,ncell)      ,cellinfo%pvar(nmet,ncell),         &
                 coords(ndimn,npoin))
        facetri(:,1) = [1,2]
        facetri(:,2) = [2,3]
        facetri(:,3) = [3,1]
        ncellface = 3
    end if

    cellinfo%nflag(:,:)     = 0
    cellinfo%neighbour(:,:) = 0

! read vertexes of every cell and its coordinates
    read(ioput,*)
    if(3. eq. ndimn) then
        do i=1,ncell
            read(ioput,*)ntemp1,ntemp2,cellinfo%vertex(1,i),cellinfo%vertex(2,i), &
                         cellinfo%vertex(3,i),cellinfo%vertex(4,i)
        end do
        read(ioput,*)
        do i=1,npoin
            read(ioput,*)coords(1,i),coords(2,i),coords(3,i)
        end do
    else if(2 .eq. ndimn) then
        do i=1,ncell
            read(ioput,*)ntemp1,ntemp2,cellinfo%vertex(1,i),cellinfo%vertex(2,i), &
                         cellinfo%vertex(3,i)
        end do
        read(ioput,*)
        do i=1,npoin
            read(ioput,*) coords(1,i),coords(2,i)
        end do
    end if

! read coordinates of every point
    
! find neighbour cells
    if(3 .eq. ndimn) then       ! 3-D meshes
        do i=1,ncell
            call progBar(i,ncell)
            if(cellinfo%neighbour(1,i)==0 .or. cellinfo%neighbour(2,i)==0 .or. &
               cellinfo%neighbour(3,i)==0 .or. cellinfo%neighbour(4,i)==0) then
                do j=1,ncell
                    if(i .ne. j) then
                        mask(:) = .FALSE.
                        do k=1,4
                            if(cellinfo%vertex(k,i)==cellinfo%vertex(1,j) .or. &
                               cellinfo%vertex(k,i)==cellinfo%vertex(2,j) .or. &
                               cellinfo%vertex(k,i)==cellinfo%vertex(3,j) .or. &
                               cellinfo%vertex(k,i)==cellinfo%vertex(4,j)) then
                                mask(k)=.TRUE.
                            end if
                        end do

                        do k=1,4
                            if(mask(facetet(1,k) .and. mask(facetet(2,k)) .and. mask(facetet(3,k))) then
                                cellinfo%neighbour(k,i)=j
                            end if
                        end do

                        if(cellinfo%neighbour(1,i)/=0 .and. cellinfo%neighbour(2,i)/=0 .and. &
                           cellinfo%neighbour(3,i)/=0 .and. cellinfo%neighbour(4,i)/=0) then
                            exit
                        end if
                    end if
                end do
            end if
        end do
    else if(2 .eq. ndimn) then
        do i=1,ncell
            call progBar(i,ncell)
            if(cellinfo%neighbour(1,i)==0 .or. cellinfo%neighbour(2,i)==0 .or. &
               cellinfo%neighbour(3,i)==0) then
                do j=1,ncell
                    if(i .ne. j) then
                        mask(:)=.FALSE.
                        do k=1,3
                            if(cellinfo%vertex(k,i)==cellinfo%vertex(1,j) .or. &
                               cellinfo%vertex(k,i)==cellinfo%vertex(2,j) .or. &
                               cellinfo%vertex(k,i)==cellinfo%vertex(3,j)) then
                                mask(k)=.TRUE.
                            end if
                        end do

                        do k=1,3
                            if(mask(facetri(1,k) .and. mask(facetri(2,k))) then
                                cellinfo%neighbour(k,i)=j
                            end if
                        end do

                        if(cellinfo%neighbour(1,i)/0 .and. cellinfo%neighbour(2,i)/=0 .and. &
                           cellinfo%neighbour(3,i)/=0) then
                            exit
                        end if
                    end if
                end do
            end if
        end do
    end if


! Boundary condition
! add dynamic arrays to store boundary cells
    write(*,*)
    read(ioput,*)
    do i=1,nbsurf
        read(ioput,*)ntemp1,ntemp2,ntemp3
        if(ntemp1 .eq. 1000) then       ! farfield
            cellinfo%neighbour(ntemp3,ntemp2) = i+ncell
            cellinfo%nflag(ntemp3,ntemp2) = 11
            cellinfo%neighbour(1,i+ncell) = ntemp2
        else if(ntemp .eq. 1001) then   ! symmetry plane
            cellinfo%neighbour(ntemp3,ntemp2) = i+ncell
            cellinfo%nflag(ntemp3,ntemp2) = 12
            cellinfo%neighbour(1,i+ncell) = ntemp2
        else if(ntemp .eq. 1001) then   ! inviscid surface
            cellinfo%neighbour(ntemp3,ntemp2) = i+ncell
            cellinfo%nflag(ntemp3,ntemp2) = 13
            cellinfo%neighbour(1,i+ncell) = ntemp2
        else if(ntemp .eq. 1001) then   ! no-slip wall
            cellinfo%neighbour(ntemp3,ntemp2) = i+ncell
            cellinfo%nflag(ntemp3,ntemp2) = 14
            cellinfo%neighbour(1,i+ncell) = ntemp2
        else 
            write(*,*) "No specified boundary type!!!"
        end if
    end do

    close(ioput)

    return 
end subroutine readmsh
