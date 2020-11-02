module io 

implicit none

contains


  subroutine writearray(fileunit,array)
    integer, intent(in) :: fileunit
    real(8), dimension(:), intent(in) :: array

    integer :: k

    do k=1,size(array)
      write(fileunit,*) array(k)
    end do

  end subroutine writearray


  subroutine write2Darray(fileunit,array)
    integer, intent(in) :: fileunit
    real(8), dimension(:,:), intent(in) :: array

    integer :: k1,k2

    do k2=1,size(array,2)
      do k1=1,size(array,1)
        write(fileunit,*) array(k1,k2)
      end do
    end do

  end subroutine write2Darray

end module io
