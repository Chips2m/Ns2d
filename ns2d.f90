program main
  use module_mesh
  use module_discr
  implicit none
  type(tmesh) :: msh
  real(kind=8),dimension(:,:),allocatable :: u,v,p,fi,dpdx,dpdy,dfidx,dfidy
  real(kind=8),dimension(:,:),allocatable :: lu,lv,div,dg
  real(kind=8) :: nu,dt
  integer :: it

  
  call tmesh_allocate(self=msh,xmin=0.0,xmax=1.0,nx=10,ymin=0.0,ymax=1.0,ny=10)

  call tmesh_allocate_var(msh,dpdx,dpdy,dfidx)
  call tmesh_allocate_var(msh,dfidy,lu,lv)
  call tmesh_allocate_var(msh,fi,div,p)
  call tmesh_allocate_var(msh,u,v,dg)



  dt = 1e-3
  nu = 1
  do it=1,100
     print*,it

     !> prediction
     call gradient(msh,p,dpdx,dpdy)
     call Linear(msh,u,v,lu,lv)
     u = u + dt*(-dpdx+nu*Lu)
     v = v + dt*(-dpdy+nu*Lv)
     call divergence(msh,u,v,div)
     div = div/dt
     
     !> solution de poisson

     !> ===================
     !> correction de pression
     call gradient(msh,fi,dfidx,dfidy)
     u = u - dt*dfidx
     v = v - dt*dfidy
     p = p + fi
  end do

  
  

end program main
