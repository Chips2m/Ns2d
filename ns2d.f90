program main
  use module_mesh
  implicit none
  type(tmesh) :: msh
  
  call tmesh_allocate(self=msh,xmin=0.0,xmax=1.0,nx=10,ymin=0.0,ymax=1.0,ny=10)
 

end program main
