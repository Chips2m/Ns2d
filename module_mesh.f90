module module_mesh


  type tmesh
     integer :: nx,ny
     real(kind=8) ::hx,hy,xmin,xmax,ymin,ymax
     real(kind=8) ,dimension(:),allocatable :: xc,xf
     real(kind=8) ,dimension(:),allocatable :: yc,yf
  end type tmesh

contains

  subroutine tmesh_allocate(self,xmin,xmax,nx,ymin,ymax,ny)
    type(tmesh)  :: self
    integer      :: nx,ny
    real(kind=8) :: xmin,xmax,ymin,ymax
    print*,"tmesh_allocate"

    self%nx=nx
    self%ny=ny
    allocate(self%xc(0:nx+1)); self%xc=0.0
    allocate(self%yc(0:ny+1))
    allocate(self%xf(0:nx+1))
    allocate(self%yf(0:ny+1))


    hx=(xmax-xmin)/nx
    hy=(ymax-ymin)/ny

    
    
    do  i=0,nx+1
       self%xc(i)= (i-0.5)*hx 
    end do
    
    do i=0,nx+1
       self%xf(i)=i*hx
    end do
    
    do j=0,ny+1
       self%yc(j)= (j-0.5)*hy 
    end do
    do j=0,ny+1
       self%yf(j)=j*hy
    end do
    
    self%hy=hy
    self%hx=hx
    self%xmin=xmin
    self%ymax=ymax
    self%xmax=xmax
    self%ymax=ymax
    
    

    print*,"x"
    
    do i=0,nx+1
       print*,i,self%xc(i),self%xf(i)
    end do
    print*,"y"
    
    do j=0,ny+1
       print*,j,self%yc(j),self%yf(j)
    end do
    
    
end subroutine tmesh_allocate




  subroutine tmesh_allocate_var(self,dg01,dg02,dg03)
    type(tmesh)  :: self
    real(kind=8) ,dimension(:,:),allocatable :: dg01,dg02,dg03
    integer :: nx,ny
    nx = self%nx
    ny = self%ny
    allocate( dg01(0:nx+1,0:ny+1) ) ; dg01 =0
    allocate( dg02(0:nx+1,0:ny+1) ) ; dg02 =0
    allocate( dg03(0:nx+1,0:ny+1) ) ; dg03 = 0
    
  end subroutine tmesh_allocate_var

end module module_mesh
