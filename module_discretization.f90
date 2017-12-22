module module_discr
  use module_mesh
contains


  subroutine gradient(mesh,fi,dx_fi,dy_fi)
    implicit none
    type(tmesh) :: mesh
    real(kind=8),dimension(:,:),allocatable :: fi,dx_fi,dy_fi,p
    integer :: i,j,k,nx,ny
    real(kind=8) :: hx,hy

    nx = mesh%nx
    ny = mesh%ny
    hx = mesh%hx
    hy = mesh%hy
    do i=0,nx
       do j=0,ny+1
          dx_fi(i,j)= (fi(i+1,j)-fi(i,j))/hx
        end do
    end do

    do i=0,nx+1
       do j=0,ny
          dy_fi(i,j)=(fi(i,j+1)-fi(i,j))/hy
        end do
    end do  	
    
  end subroutine gradient

  subroutine divergence(mesh,u,v,div)
    implicit none
    type(tmesh) :: mesh
    real(kind=8),dimension(:,:),allocatable :: u,v,div
    integer :: i,j,nx,ny
    real(kind=8) :: hx,hy
    
    nx = mesh%nx
    ny = mesh%ny
    hx = mesh%hx
    hy = mesh%hy
    do i=1,nx
       do j=1,ny
          div(i,j)= (u(i,j)-u(i-1,j))/hx + (v(i,j)-v(i,j-1))/hy
        end do
     end do
     
   end subroutine divergence

   subroutine Linear(mesh,u,v,lu,lv)
    implicit none
    type(tmesh) :: mesh
    real(kind=8),dimension(:,:),allocatable :: u,v,lu,lv
    integer :: i,j,nx,ny
    real(kind=8) :: hx2,hy2
    
    nx = mesh%nx
    ny = mesh%ny
    hx2 = mesh%hx**2
    hy2 = mesh%hy**2
    do i=1,nx-1
       do j=1,ny
          lu(i,j)=&
               + (u(i-1,j  )-2*u(i,j)+u(i+1,j))/hx2 &
               + (u(i  ,j-1)-2*u(i,j)+u(i,j+1))/hy2
       end do
    end do
    do i=1,nx
       do j=1,ny-1
          lv(i,j)=&
               + (v(i-1,j  )-2*v(i,j)+v(i+1,j))/hx2 &
               + (v(i  ,j-1)-2*v(i,j)+v(i,j+1))/hy2
       end do
    end do
    
    
  end subroutine Linear


end module module_discr
