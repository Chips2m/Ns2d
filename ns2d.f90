program main
  use module_mesh
  use module_discr
  use module_solver
  implicit none
  type(tmesh) :: msh
  type(tsolver) :: poisson
  real(kind=8),dimension(:,:),allocatable :: u,v,p,fi,dpdx,dpdy,dfidx,dfidy
  real(kind=8),dimension(:,:),allocatable :: lu,lv,div,dg
  real(kind=8) :: nu,dt
  integer :: it
  integer :: nx,ny
  
  nx = 20
  ny = nx
  call tmesh_allocate(self=msh,xmin=0d0,xmax=1d0,nx=nx,ymin=0d0,ymax=1d0,ny=ny)

  call tmesh_allocate_var(msh,dpdx,dpdy,dfidx)
  call tmesh_allocate_var(msh,dfidy,lu,lv)
  call tmesh_allocate_var(msh,fi,div,p)
  call tmesh_allocate_var(msh,u,v,dg)
  call tsolver_allocate(self=poisson,mesh=msh)

  
  dt = 0.2*msh%hx**2    !> pas de temps 
  nu = 1d0     !> viscosité du fluide
  do it=1,1000 !> avancement temporel

     !> etape de prediction
     call set_bc()
     call gradient(msh,p,dpdx,dpdy)
     call Linear(msh,u,v,lu,lv)
     u(1:nx-1,1:ny) = u(1:nx-1,1:ny) + dt*( -dpdx(1:nx-1,1:ny) + nu*Lu(1:nx-1,1:ny) )
     v(1:nx,1:ny-1) = v(1:nx,1:ny-1) + dt*( -dpdy(1:nx,1:ny-1) + nu*Lv(1:nx,1:ny-1) )
     
     !> divergence du champs de vitesse intermédiaire
     call divergence(msh,u,v,div)
     div = div/dt
     !> résolution de l'équation de poisson
     call tsolver_solve(poisson,msh,div,fi)
     
     !> ===================
     !> correction du champ de vitesse
     call gradient(msh,fi,dfidx,dfidy)
     u(1:nx-1,1:ny) = u(1:nx-1,1:ny) - dt*dfidx(1:nx-1,1:ny)
     v(1:nx,1:ny-1) = v(1:nx,1:ny-1) - dt*dfidy(1:nx,1:ny-1)
     p = p + fi
     !>  On affiche le maximum de la correction apportée max(dfidx,dfidy) et la divergence
     call divergence(msh,u,v,div)
     print*,it,max(maxval(abs(dfidx)),maxval(abs(dfidy))),maxval(abs(div))
  end do
  
  !> pour l'équipe maillage mon postraitement au format tecplot a vous de faire le votre!!!
  call export_tecplot()
contains

  subroutine set_bc()
    implicit none
    integer :: i,j
    do i=1,msh%nx
       u(i,0) = 2*1 - u(i,1)
    end do
  end subroutine set_bc
  
  subroutine export_tecplot()
    implicit none
    integer ::i,j
    write(22,*)'VARIABLES= X,Y,U,V,P'
    write(22,*)'ZONE I=',msh%nx,', J=',msh%ny,', DATAPACKING=POINT'
    do j=1,msh%ny
       do i=1,msh%nx
          write(22,*),msh%xc(i),msh%yc(j),&
               (u(i,j)+u(i-1,j))*0.5,&
               (v(i,j-1)+v(i,j))*0.5,p(i,j)
       end do
    end do
    
    
  end subroutine export_tecplot




end program main
