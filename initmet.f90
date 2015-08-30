subroutine initmet
    use globalvar
    implicit none 
!********************************************************************************
!
! Function: calculate metrics, cell volumns
!
!********************************************************************************

! Local variables
    integer :: i,j,vera,verb,verc,adcell
    real*8  :: verxa,verxb,verxc,verya,veryb,veryc,verza,verzb,verzc

! Compute face norms and cell volumns
    cellinfo%vol(:) = 0.0D0
    if( 3 .eq. ndmin) then      ! 3-D (Gauss' Theorem)
        do i=1,ncell
            do j=1,4
                vera = cellinfo%vertex(facetet(1,j),i)
                verb = cellinfo%vertex(facetet(2,j),i)
                verc = cellinfo%vertex(facetet(3,j),i)

                verxa = coords(1,verb)-coords(1,vera)
                verya = coords(2,verb)-coords(2,vera)
                verza = coords(3,verb)-coords(3,vera)

                verxb = coords(1,verc)-coords(1,vera)
                veryb = coords(2,verc)-coords(2,vera)
                verzb = coords(3,verc)-coords(3,vera)

                verxc = 0.50D0*(verya*verzb-veryb*verza)
                veryc = 0.50D0*(verxb*verza-verxa*verzb)
                verzc = 0.50D0*(verxa*veryb-verxb*verya)
    !   cell volume
                cellinfo%vol(i) = cellinfo%vol(i)+verxa*verxc+verya*veryc+verza*verzc

    !   face norm
                cellinfo%facenorm(1,j,i) = verxc
                cellinfo%facenorm(2,j,i) = veryc
                cellinfo%facenorm(3,j,i) = verzc
            end do
        end do
    else if( 2 .eq. ndmin ) then   ! 2-D (outer product, anther method: Green's Theorem)
        do i=1,ncell
            vera = cellinfo%vertex(1,i)
            verb = cellinfo%vertex(2,i)
            verc = cellinfo%vertex(3,i)
            
            verxa = coords(1,vera)
            verya = coords(2,vera)
            verxb = coords(1,verb)
            veryb = coords(2,verb)
            verxc = coords(1,verc)
            veryc = coords(2,verc)
            cellinfo%centroid(1,i)=(verxa+verxb+verxc)/3.0D0
            cellinfo%centroid(2,i)=(verya+veryb+veryc)/3.0D0
            cellinfo%vol(i) = 0.50D0*(verxb*veryc+verxa*veryb+verxc*verya- &
                                      verxb*verya-verxc*veryb-verxa*veryc)
        end do

        do i=1,ncell
            do i=1,3
                vera = cellinfo%vertex(facetri(1,j),i)
                verb = cellinfo%vertex(facetri(2,j),i)
                verxa = coords(1,vera)
                verya = coords(2,vera)
                verxb = coords(1,verb)
                veryb = coords(2,verb)

        ! compute vectors from the cell-centroid to the face-midpoint
                verxc = 0.50D0*(verxa+verxb)-cellinfo%centroid(1,i)
                veryc = 0.50D0*(verya+veryb)-cellinfo%centroid(2,i)
                cellinfo%cetomid(1,j,i) = verxc
                cellinfo%cetomid(2,j,i) = veryc

                verxc = verxb-verxa
                veryc = veryb-verya
                cellinfo%facenorm(1,j,i) = veryc
                cellinfo%facenorm(2,j,i) = verxc
            end do
        end do
    end if

! Add geometrical informination to ghost cells
    do i=ncell+1,nvitgho
        adcell = cellinfo%neighbour(1,i)
        cellinfo%vol(i) = cellinfo%vol(adcell)
    end do

    return 
end subroutine initmet
