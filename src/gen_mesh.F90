program gen_mesh

  use fiona
  use delaunay_voronoi_mod
  use mpas_mesh_mod

  implicit none

  character(1024) grids_file_path
  integer nx
  real(8), allocatable :: x(:,:)

  call get_command_argument(1, grids_file_path)

  call fiona_init()
  call fiona_open_dataset('grids', file_path=grids_file_path)
  call fiona_get_dim('grids', 'location_nv', size=nx)

  allocate(x(nx,3))

  call fiona_start_input('grids')
  call fiona_input('grids', 'vtx_p', x)
  call fiona_end_input('grids')

  call delaunay_voronoi_init(nx, x=x(:,1), y=x(:,2), z=x(:,3))
  call delaunay_triangulation()
  call voronoi_diagram(all=.true.)
  call delaunay_voronoi_output(mpas_mesh_output)

  deallocate(x)

end program gen_mesh
