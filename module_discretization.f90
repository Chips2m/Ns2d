module discr

contains


  subroutine gradient(x_c,y_c,fi,dx_fi,dy_fi,n_x,n_y)
    implicit none
    real(kind=8),dimension(:)  ,allocatable :: x_c,y_c
    real(kind=8),dimension(:,:),allocatable :: fi,dx_fi,dy_fi,p
    integer :: i,j,k,n_x,n_y
    real(kind=8) :: h_x,h_y

    do i=0,n_x
	do j=0,n_y+1
	dx_fi(i,j)=(p(i+1,j)-p(i,j))/(x_c(i+1)-x_c(i))
        end do
    end do

    do i=0,n_x
	do j=0,n_y+1
    
	dy_fi(i,j)=(p(i,j+1)-p(i,j))/(y_c(j+1)-y_c(j))
        end do
    end do  	
    
    

   end subroutine gradient



   subroutine divergence(x_f,y_f,d_v,n_x,n_y)
     implicit none
     real(kind=8),dimension(:)  ,allocatable :: x_f,y_f
     real(kind=8),dimension(:,:),allocatable :: d_v,u,v	
	  integer :: i,j,k,n_x,n_y
     real(kind=8) :: h_x,h_y

    do i=0,n_x
	do j=0,n_y+1
	d_v(i,j)=(  u(i,j) - u(i-1,j) ) / ( x_f(i+1) - x_f(i)  ) + ( v(i+1,j) - v(i,j)  ) / ( y_f(j+1) - y_f(j)  )
        
	end do
    end do
    end subroutine divergence



	subroutine linear(u,v,x_f,y_f,dx_l,dy_l,n_x,n_y)
     implicit none
     real(kind=8),dimension(:)  ,allocatable :: x_f,y_f
     real(kind=8),dimension(:,:),allocatable :: dx_l,dy_l,u,v
	  integer :: i,j,k,n_x,n_y
     real(kind=8) :: h_x,h_y
 
	     do i=0,n_x
		do j=0,n_y+1
		dx_l(i,j)=(  u(i+1,j) - ( 2 * u(i,j) ) + u(i-1,j) ) / ( x_f(i+1) - x_f(i) ) ** 2 
		
		end do
	    end do

	    do i=0,n_x
		do j=0,n_y+1
		dy_l(i,j)=( u(j,i+1) - ( 2 * u(i,j)  ) + u(j-1,i) ) / ( ( x_f(j+1) - x_f(j) ) ** 2 )
		
		end do
	    end do


	end subroutine linear 


end module discr
