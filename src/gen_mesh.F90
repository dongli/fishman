program gen_mesh

  use io_mod
  use delaunay_voronoi_mod
  use mpas_mesh_mod

  implicit none

  character(1024) grids_file_path
  integer nx
  real(8), allocatable :: x(:,:)

  call get_command_argument(1, grids_file_path)

  call io_init()
  call io_create_dataset('grids', file_path=grids_file_path, mode='input')
  call io_get_dim('grids', 'location_nv', size=nx)

  allocate(x(nx,3))

  call io_start_input('grids')
  call io_input('grids', 'vtx_p', x)
  call io_end_input('grids')

  call delaunay_voronoi_init(nx, x=x(:,1), y=x(:,2), z=x(:,3))
  call delaunay_triangulation()
  call voronoi_diagram(all=.true.)
  call delaunay_voronoi_output(mpas_mesh_output)

  deallocate(x)

end program gen_mesh
