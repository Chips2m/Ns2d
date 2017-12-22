module module_solver
  use module_mesh

  type tsolver
     integer :: nrows,nnz,NX
     integer,dimension(:),allocatable :: ia,ja,VV
     real(kind=8),dimension(:),allocatable :: fi,x,b
  end type tsolver
  
contains
  
  
  
  subroutine array2vec(mesh,u,fi)
    implicit none
    type(tmesh)    :: mesh
    real(kind=8),dimension(:,:),allocatable :: u
    real(kind=8),dimension(:),allocatable :: fi
    integer :: i,j,nx,ny,p
    
    p=0
    do j=1,mesh%ny
       do i=1,mesh%nx
          p=p+1
          fi(p)=u(i,j)
       end do
    end do
  end subroutine array2vec
  
  subroutine vec2array(mesh,fi,u)
    implicit none
    type(tmesh)    :: mesh
    real(kind=8),dimension(:),allocatable :: fi
    real(kind=8),dimension(:,:),allocatable :: u

    integer :: i,j,nx,ny,p
    
    u=0
    p=0
    do j=1,mesh%ny
       do i=1,mesh%nx
          p=p+1
          u(i,j)=fi(p)
       end do
    end do
  end subroutine vec2array
  

  subroutine Set_neumann_bc(mesh,u)
    implicit none
    type(tmesh)    :: mesh
    real(kind=8),dimension(:,:),allocatable :: u

    
    u(1:mesh%nx,       0)  = u(1:mesh%nx,     1)
    u(1:mesh%nx,mesh%nx+1) = u(1:mesh%nx,mesh%nx)
    u(        0,1:mesh%ny) = u(      1,1:mesh%ny)
    u(mesh%nx+1,1:mesh%ny) = u(mesh%nx,1:mesh%ny)
    
  end subroutine Set_neumann_bc
  

  subroutine tsolver_solve(self,mesh,s,fi)
    implicit none
    type(tsolver)  :: self
    type(tmesh)    :: mesh
    real(kind=8),dimension(:,:),allocatable :: fi,s
    integer iter,iprint
    real(kind=8) :: tol
    call array2vec(mesh,s,self%b)

    iter = 50
    tol = 1e-12
    iprint=0
    call dagmg(self%nrows,self%fi,self%ja,self%ia,self%b,self%x,0,iprint,10,iter,tol)
    
    call vec2array(mesh,self%x,fi)
    call Set_neumann_bc(mesh,fi)
    


  end subroutine tsolver_solve
  
  subroutine tsolver_allocate(self,mesh)
    implicit none
    type(tsolver)  :: self
    type(tmesh)    :: mesh
    integer        :: nx,ny
    integer :: i,j
    integer :: tmp
    integer :: nnz,nrows
    integer,dimension(5) ::col
    real(kind=8),dimension(5) ::val
    real(kind=8) :: ap,aw,ae,as,an
    integer       :: jp,jw,je,js,jn
    real(kind=8) :: hx,hy
    
    print*,'tsolver_allocate'
    nx = mesh%nx
    ny = mesh%ny
    nrows=nx*ny

    
    nnz=0
    do j=1,ny
       do i=1,nx
          tmp=5
          if (i==1)  tmp = tmp-1
          if (i==nx) tmp = tmp-1
          if (j==1)  tmp = tmp-1
          if (j==ny) tmp = tmp-1
          nnz=nnz+tmp
          
!          print*,i,j,nnz
       end do
    end do

    
    self%nnz=nnz
    self%nrows=nrows
    allocate(self%ia(1:nrows+1))
    allocate(self%ja(1:nnz))
    allocate(self%fi(1:nnz))

    allocate(self%x(1:nrows))
    allocate(self%b(1:nrows))
    
    hx=mesh%hx
    hy=mesh%hy
    nnz=0
    self%ia(1) = 1
    do j=1,ny
       do i=1,nx
          tmp=5
          
          jp = i + (j-1)*nx
          jw = jp-1
          je = jp+1
          js = jp-nx
          jn = jp+nx
          ap = -2/hx**2-2/hy**2
          aw = +1/hx**2
          ae = +1/hx**2
          as = +1/hy**2
          an = +1/hy**2
          
          if (i==1)  then
             jw = -1
             ap = ap + aw
             tmp = tmp-1
          end if
          if (i==nx) then
             je = -1
             ap = ap + ae
             tmp = tmp-1
          end if
          if (j==1)  then
             js = -1
             ap = ap + as
             tmp = tmp-1
          end if
          if (j==ny) then
             jn = -1
             ap = ap + an
             tmp = tmp-1
          end if
          
          if ((i==1).and.(j==1))  then
             ap = ap - 2*aw 
          end if
          col=(/jp,jw,je,js,jn/)
          val=(/ap,aw,ae,as,an/)
          
!!$          print*,PACK( col, col>0 )
!!$          print*,PACK( val, col>0 )
          
          self%ja(nnz+1:nnz+tmp) = PACK( col, col>0 )
          self%fi(nnz+1:nnz+tmp) = PACK( val, col>0 )
          nnz=nnz+tmp
          self%ia(jp+1) =  nnz+1
       end do
    end do
!!$    print*,self%ia
    
  end subroutine tsolver_allocate
  
  

  
end module module_solver
