module mpas_mesh_mod

  use array_mod
  use linked_list_mod
  use const_mod
  use log_mod
  use string_mod
  use io_mod
  use delaunay_voronoi_mod
  use sphere_geometry_mod

  implicit none

  private

  public mpas_mesh_run
  public mpas_mesh_output
  public mpas_mesh_type

  integer, parameter :: vertexDegree = 3
  integer, parameter :: maxEdges = 15
  integer, parameter :: maxEdges2 = maxEdges * 2

  type mpas_mesh_type
    ! Dimensions
    integer nCells
    integer nEdges
    integer nVertices
    ! Coordinates
    real(8), allocatable :: lonCell(:)
    real(8), allocatable :: latCell(:)
    real(8), allocatable :: xCell(:)
    real(8), allocatable :: yCell(:)
    real(8), allocatable :: zCell(:)
    real(8), allocatable :: lonEdge(:)
    real(8), allocatable :: latEdge(:)
    real(8), allocatable :: xEdge(:)
    real(8), allocatable :: yEdge(:)
    real(8), allocatable :: zEdge(:)
    real(8), allocatable :: lonVertex(:)
    real(8), allocatable :: latVertex(:)
    real(8), allocatable :: xVertex(:)
    real(8), allocatable :: yVertex(:)
    real(8), allocatable :: zVertex(:)
    ! Geometric measures
    real(8), allocatable :: dvEdge(:)               ! Distance in meters between the vertices that saddle a given edge
    real(8), allocatable :: dv1Edge(:)              ! Distance in meters between vertex 1 and edge point
    real(8), allocatable :: dv2Edge(:)              ! Distance in meters between vertex 2 and edge point
    real(8), allocatable :: dcEdge(:)               ! Distance in meters between the cells that saddle a given edge
    real(8), allocatable :: areaCell(:)             ! Area in square meters for a given cell of the primary mesh
    real(8), allocatable :: areaTriangle(:)         ! Area in square meters for a given triangle of the dual mesh
    real(8), allocatable :: kiteAreasOnVertex(:,:)  ! The intersection area of areaTriangle with each cell that radiates from a given vertex
    real(8), allocatable :: angleEdge(:)            ! Angle in radians an edge’s normal vector makes with the local eastward direction
    ! Indices
    integer, allocatable :: nEdgesOnCell(:)         ! Number of edges on a given cell
    integer, allocatable :: nEdgesOnEdge(:)         ! Number of edges on a given edge. Used to reconstruct tangential velocities.
    integer, allocatable :: indexToCellID(:)        ! Global cell ID for all cell centers
    integer, allocatable :: indexToEdgeID(:)        ! Global edge ID for all edge locations
    integer, allocatable :: indexToVertexID(:)      ! Global vertex ID for all cell vertices
    integer, allocatable :: cellsOnCell(:,:)        ! Cell indices that surround a given cell
    integer, allocatable :: cellsOnEdge(:,:)        ! Cell indices that saddle a given edge
    integer, allocatable :: cellsOnVertex(:,:)      ! Cell indices that radiate from a given vertex
    integer, allocatable :: edgesOnCell(:,:)        ! Edge indices that surround a given cell
    integer, allocatable :: edgesOnEdge(:,:)        ! Edge indices that are used to reconstruct tangential velocities
    integer, allocatable :: edgesOnVertex(:,:)      ! Edge indices that radiate from a given vertex
    integer, allocatable :: verticesOnCell(:,:)     ! Vertex indices that surround a given cell
    integer, allocatable :: verticesOnEdge(:,:)     ! Vertex indices that saddle a given edge
    ! Weights
    real(8), allocatable :: weightsOnEdge(:,:)      ! Weights to reconstruct tangential velocities
  contains
    final :: mpas_mesh_finalize
  end type mpas_mesh_type

contains

  subroutine mpas_mesh_run(DVT_array, DT_list, DE_array, VVT_array, VC_array, VE_array, mesh)

    type(array_type), intent(in) :: DVT_array
    type(linked_list_type), intent(in) :: DT_list
    type(array_type), intent(in) :: DE_array
    type(array_type), intent(in) :: VVT_array
    type(array_type), intent(in) :: VC_array
    type(array_type), intent(in) :: VE_array
    type(mpas_mesh_type), intent(inout) :: mesh

    integer iCell, iEdge, iVertex, iVertexOnVertex, i, j, i1, i2
    real(8) x(maxEdges), y(maxEdges), z(maxEdges)
    real(8) total_area
    type(delaunay_triangle_type), pointer :: DT
    type(delaunay_vertex_type), pointer :: DVT
    type(delaunay_edge_type), pointer :: DE
    type(voronoi_cell_type), pointer :: VC, adjVC
    type(voronoi_vertex_type), pointer :: VVT
    type(voronoi_edge_type), pointer :: VE, VE1, VE2
    type(linked_list_iterator_type) iterator
    type(point_type) xp

    mesh%nCells = VC_array%size
    mesh%nEdges = VE_array%size
    mesh%nVertices = VVT_array%size

    allocate(mesh%latCell          (mesh%nCells))
    allocate(mesh%lonCell          (mesh%nCells))
    allocate(mesh%xCell            (mesh%nCells))
    allocate(mesh%yCell            (mesh%nCells))
    allocate(mesh%zCell            (mesh%nCells))
    allocate(mesh%indexToCellID    (mesh%nCells))
    allocate(mesh%lonEdge          (mesh%nEdges))
    allocate(mesh%latEdge          (mesh%nEdges))
    allocate(mesh%xEdge            (mesh%nEdges))
    allocate(mesh%yEdge            (mesh%nEdges))
    allocate(mesh%zEdge            (mesh%nEdges))
    allocate(mesh%indexToEdgeID    (mesh%nEdges))
    allocate(mesh%lonVertex        (mesh%nVertices))
    allocate(mesh%latVertex        (mesh%nVertices))
    allocate(mesh%xVertex          (mesh%nVertices))
    allocate(mesh%yVertex          (mesh%nVertices))
    allocate(mesh%zVertex          (mesh%nVertices))
    allocate(mesh%indexToVertexID  (mesh%nVertices))
    allocate(mesh%nEdgesOnCell     (mesh%nCells))
    allocate(mesh%nEdgesOnEdge     (mesh%nEdges))
    allocate(mesh%cellsOnCell      (maxEdges,mesh%nCells))
    allocate(mesh%cellsOnEdge      (2,mesh%nEdges))
    allocate(mesh%edgesOnCell      (maxEdges,mesh%nCells))
    allocate(mesh%edgesOnEdge      (maxEdges2,mesh%nEdges))
    allocate(mesh%verticesOnCell   (maxEdges,mesh%nCells))
    allocate(mesh%verticesOnEdge   (2,mesh%nEdges))
    allocate(mesh%edgesOnVertex    (vertexDegree,mesh%nVertices))
    allocate(mesh%cellsOnVertex    (vertexDegree,mesh%nVertices))
    allocate(mesh%weightsOnEdge    (maxEdges2,mesh%nEdges))
    allocate(mesh%dvEdge           (mesh%nEdges))
    allocate(mesh%dv1Edge          (mesh%nEdges))
    allocate(mesh%dv2Edge          (mesh%nEdges))
    allocate(mesh%dcEdge           (mesh%nEdges))
    allocate(mesh%angleEdge        (mesh%nEdges))
    allocate(mesh%areaCell         (mesh%nCells))
    allocate(mesh%areaTriangle     (mesh%nVertices))
    allocate(mesh%kiteAreasOnVertex(vertexDegree,mesh%nVertices))

    do iCell = 1, mesh%nCells
      VC => get_VC(VC_array, iCell)
      mesh%lonCell(iCell) = VC%center%lon
      mesh%latCell(iCell) = VC%center%lat
      mesh%  xCell(iCell) = VC%center%x
      mesh%  yCell(iCell) = VC%center%y
      mesh%  zCell(iCell) = VC%center%z
      do i = 1, VC%VVT%size
        VVT => get_VVT(VC%VVT, i)
        mesh%verticesOnCell(i,iCell) = VVT%id
      end do
      do i = 1, VC%adjVC%size
        adjVC => get_VC(VC%adjVC, i)
        mesh%cellsOnCell(i,iCell) = adjVC%id
      end do
      do i = 1, VC%VE%size
        VE => get_VE(VC%VE, i)
        mesh%edgesOnCell(i,iCell) = VE%id
      end do
      mesh%nEdgesOnCell(iCell) = VC%VE%size
      mesh%indexToCellID(iCell) = VC%id
    end do

    do iVertex = 1, mesh%nVertices
      VVT => get_VVT(VVT_array, iVertex)
      mesh%lonVertex(iVertex) = VVT%lon
      mesh%latVertex(iVertex) = VVT%lat
      mesh%  xVertex(iVertex) = VVT%x
      mesh%  yVertex(iVertex) = VVT%y
      mesh%  zVertex(iVertex) = VVT%z
      do i = 1, VVT%DT%DVT%size
        DVT => get_DVT(VVT%DT%DVT, i)
        mesh%cellsOnVertex(i,iVertex) = DVT%id
      end do
      do i = 1, VVT%VE%size
        VE => get_VE(VVT%VE, i)
        mesh%edgesOnVertex(i,iVertex) = VE%id
      end do
      mesh%indexToVertexID(iVertex) = VVT%id
    end do

    do iEdge = 1, mesh%nEdges
      DE => get_DE(DE_array, iEdge)
      VE => get_VE(VE_array, iEdge)
      if (intersect(DE%DVT1, DE%DVT2, VE%VVT1, VE%VVT2, xp)) then
        mesh%lonEdge(iEdge) = xp%lon
        mesh%latEdge(iEdge) = xp%lat
        mesh%  xEdge(iEdge) = xp%x
        mesh%  yEdge(iEdge) = xp%y
        mesh%  zEdge(iEdge) = xp%z
      else
        call log_error('Unable to find intersection between DE and VE!')
      end if
      mesh%dvEdge(iEdge) = calc_arc_length(VE%VVT1, VE%VVT2)
      mesh%dcEdge(iEdge) = calc_arc_length(DE%DVT1, DE%DVT2)
      mesh%dv1Edge(iEdge) = calc_arc_length(VE%VVT1, [mesh%xEdge(iEdge),mesh%yEdge(iEdge),mesh%zEdge(iEdge)])
      mesh%dv2Edge(iEdge) = calc_arc_length(VE%VVT2, [mesh%xEdge(iEdge),mesh%yEdge(iEdge),mesh%zEdge(iEdge)])
      xp%lon = DE%DVT1%lon + 0.1 * pi
      if (xp%lon > pi2) xp%lon = xp%lon - pi2
      xp%lat = DE%DVT1%lat
      call cartesian_transform(xp)
      mesh%angleEdge(iEdge) = calc_sphere_angle(xp, DE%DVT1, DE%DVT2)
      mesh%verticesOnEdge(1,iEdge) = VE%VVT1%id
      mesh%verticesOnEdge(2,iEdge) = VE%VVT2%id
      mesh%cellsOnEdge(1,iEdge) = VE%VC1%id
      mesh%cellsOnEdge(2,iEdge) = VE%VC2%id
      mesh%nEdgesOnEdge(iEdge) = VE%VC1%VE%size + VE%VC2%VE%size - 2
      mesh%indexToEdgeID(iEdge) = VE%id
      j = 1
      i1 = VE%VC1%VE%index_ptr(VE)
      do i = i1 + 1, VE%VC1%VE%size
        VE1 => get_VE(VE%VC1%VE, i)
        mesh%edgesOnEdge(j,iEdge) = VE1%id
        j = j + 1
      end do
      do i = 1, i1 - 1
        VE1 => get_VE(VE%VC1%VE, i)
        mesh%edgesOnEdge(j,iEdge) = VE1%id
        j = j + 1
      end do
      i2 = VE%VC2%VE%index_ptr(VE)
      do i = i2 + 1, VE%VC2%VE%size
        VE2 => get_VE(VE%VC2%VE, i)
        mesh%edgesOnEdge(j,iEdge) = VE2%id
        j = j + 1
      end do
      do i = 1, i2 - 1
        VE2 => get_VE(VE%VC2%VE, i)
        mesh%edgesOnEdge(j,iEdge) = VE2%id
        j = j + 1
      end do
    end do

    total_area = 0.0
    do iCell = 1, mesh%nCells
      VC => get_VC(VC_array, iCell)
      do i = 1, VC%VVT%size
        VVT => get_VVT(VC%VVT, i)
        x(i) = VVT%x; y(i) = VVT%y; z(i) = VVT%z
      end do
      mesh%areaCell(iCell) = calc_area(x(:VC%VVT%size), y(:VC%VVT%size), z(:VC%VVT%size))
      total_area = total_area + mesh%areaCell(iCell)
    end do
    if (abs(total_area - 4 * pi) > 1.0d-10) then
      call log_error('Too large total cell area error!', __FILE__, __LINE__)
    end if

    total_area = 0.0
    do iVertex = 1, mesh%nVertices
      VVT => get_VVT(VVT_array, iVertex)
      do i = 1, VVT%DT%DVT%size
        DVT => get_DVT(VVT%DT%DVT, i)
        x(i) = DVT%x; y(i) = DVT%y; z(i) = DVT%z
      end do
      mesh%areaTriangle(iVertex) = calc_area(x(:VVT%DT%DVT%size), y(:VVT%DT%DVT%size), z(:VVT%DT%DVT%size))
      total_area = total_area + mesh%areaTriangle(iVertex)
    end do
    if (abs(total_area - 4 * pi) > 1.0d-10) then
      call log_error('Too large total triangle area error!', __FILE__, __LINE__)
    end if

    total_area = 0.0
    do iVertex = 1, mesh%nVertices
      VVT => get_VVT(VVT_array, iVertex)
      do i = 1, VVT%VE%size
        VE1 => get_VE(VVT%VE, i)
        VE2 => get_VE(VVT%VE, merge(1, i + 1, i == VVT%VE%size))
        VC => get_VC(VVT%VC, i)
        x(1) = mesh%xEdge(VE1%id); y(1) = mesh%yEdge(VE1%id); z(1) = mesh%zEdge(VE1%id)
        x(2) = mesh%xCell(VC %id); y(2) = mesh%yCell(VC %id); z(2) = mesh%zCell(VC %id)
        x(3) = mesh%xEdge(VE2%id); y(3) = mesh%yEdge(VE2%id); z(3) = mesh%zEdge(VE2%id)
        x(4) = VVT%x;              y(4) = VVT%y;              z(4) = VVT%z
        mesh%kiteAreasOnVertex(i,iVertex) = calc_area(x(:4), y(:4), z(:4))
        total_area = total_area + mesh%kiteAreasOnVertex(i,iVertex)
      end do
    end do
    if (abs(total_area - 4 * pi) > 1.0d-10) then
      call log_error('Too large total kite area error!', __FILE__, __LINE__)
    end if

    call calc_weightsOnEdge(VVT_array, VC_array, VE_array, mesh)

  end subroutine mpas_mesh_run

  subroutine mpas_mesh_output(mesh)

    type(mpas_mesh_type), intent(in) :: mesh

    call io_create_dataset(name='mpas_mesh', file_path='mpas_mesh.' // trim(to_string(mesh%nCells)) // '.nc')
    call io_add_att('on_a_sphere',   'YES',              dataset_name='mpas_mesh')
    call io_add_att('sphere_radius', 1.0d0,              dataset_name='mpas_mesh')
    call io_add_att('mesh_id',       'demo',             dataset_name='mpas_mesh')
    call io_add_att('mesh',     '1.0',              dataset_name='mpas_mesh')
    call io_add_dim('nCells',       size=mesh%nCells,    dataset_name='mpas_mesh')
    call io_add_dim('nEdges',       size=mesh%nEdges,    dataset_name='mpas_mesh')
    call io_add_dim('nVertices',    size=mesh%nVertices, dataset_name='mpas_mesh')
    call io_add_dim('vertexDegree', size=vertexDegree,   dataset_name='mpas_mesh')
    call io_add_dim('maxEdges',     size=maxEdges,       dataset_name='mpas_mesh')
    call io_add_dim('maxEdges2',    size=maxEdges2,      dataset_name='mpas_mesh')
    call io_add_dim('TWO',          size=2,              dataset_name='mpas_mesh')
    call io_add_dim('Time',                              dataset_name='mpas_mesh')
    call io_add_var('latCell',           long_name='Latitude of all cell centers',                                    units='radian', dim_names=['nCells      '],                 data_type='real(8)', dataset_name='mpas_mesh')
    call io_add_var('lonCell',           long_name='Longitude of all cell centers',                                   units='radian', dim_names=['nCells      '],                 data_type='real(8)', dataset_name='mpas_mesh')
    call io_add_var('xCell',             long_name='x axis position of all cell centers',                             units='m',      dim_names=['nCells      '],                 data_type='real(8)', dataset_name='mpas_mesh')
    call io_add_var('yCell',             long_name='y axis position of all cell centers',                             units='m',      dim_names=['nCells      '],                 data_type='real(8)', dataset_name='mpas_mesh')
    call io_add_var('zCell',             long_name='z axis position of all cell centers',                             units='m',      dim_names=['nCells      '],                 data_type='real(8)', dataset_name='mpas_mesh')
    call io_add_var('indexToCellID',     long_name='Global cell ID for all cell centers',                             units='',       dim_names=['nCells      '],                 data_type='integer', dataset_name='mpas_mesh')
    call io_add_var('latEdge',           long_name='Latitude of all edge locations',                                  units='radian', dim_names=['nEdges      '],                 data_type='real(8)', dataset_name='mpas_mesh')
    call io_add_var('lonEdge',           long_name='Longitude of all edge locations',                                 units='radian', dim_names=['nEdges      '],                 data_type='real(8)', dataset_name='mpas_mesh')
    call io_add_var('xEdge',             long_name='x axis position of all edge locations',                           units='m',      dim_names=['nEdges      '],                 data_type='real(8)', dataset_name='mpas_mesh')
    call io_add_var('yEdge',             long_name='y axis position of all edge locations',                           units='m',      dim_names=['nEdges      '],                 data_type='real(8)', dataset_name='mpas_mesh')
    call io_add_var('zEdge',             long_name='z axis position of all edge locations',                           units='m',      dim_names=['nEdges      '],                 data_type='real(8)', dataset_name='mpas_mesh')
    call io_add_var('indexToEdgeID',     long_name='Global edge ID for all edge locations',                           units='',       dim_names=['nEdges      '],                 data_type='integer', dataset_name='mpas_mesh')
    call io_add_var('latVertex',         long_name='Latitude of all cell vertices',                                   units='radian', dim_names=['nVertices   '],                 data_type='real(8)', dataset_name='mpas_mesh')
    call io_add_var('lonVertex',         long_name='Longitude of all cell vertices',                                  units='radian', dim_names=['nVertices   '],                 data_type='real(8)', dataset_name='mpas_mesh')
    call io_add_var('xVertex',           long_name='x axis position of all cell vertices',                            units='m',      dim_names=['nVertices   '],                 data_type='real(8)', dataset_name='mpas_mesh')
    call io_add_var('yVertex',           long_name='y axis position of all cell vertices',                            units='m',      dim_names=['nVertices   '],                 data_type='real(8)', dataset_name='mpas_mesh')
    call io_add_var('zVertex',           long_name='z axis position of all cell vertices',                            units='m',      dim_names=['nVertices   '],                 data_type='real(8)', dataset_name='mpas_mesh')
    call io_add_var('indexToVertexID',   long_name='Global vertex ID for all cell vertices',                          units='',       dim_names=['nVertices   '],                 data_type='integer', dataset_name='mpas_mesh')
    call io_add_var('nEdgesOnCell',      long_name='Number of edges on a given cell',                                 units='',       dim_names=['nCells      '],                 data_type='integer', dataset_name='mpas_mesh')
    call io_add_var('nEdgesOnEdge',      long_name='Number of edges on a given edge',                                 units='',       dim_names=['nEdges      '],                 data_type='integer', dataset_name='mpas_mesh')
    call io_add_var('cellsOnCell',       long_name='Cell indices that surround a given cell',                         units='',       dim_names=['maxEdges    ', 'nCells      '], data_type='integer', dataset_name='mpas_mesh')
    call io_add_var('cellsOnEdge',       long_name='Cell indices that saddle a given edge',                           units='',       dim_names=['TWO         ', 'nEdges      '], data_type='integer', dataset_name='mpas_mesh')
    call io_add_var('edgesOnCell',       long_name='Edge indices that surround a given cell',                         units='',       dim_names=['maxEdges    ', 'nCells      '], data_type='integer', dataset_name='mpas_mesh')
    call io_add_var('edgesOnEdge',       long_name='Edge indices that are used to reconstruct tangential velocities', units='',       dim_names=['maxEdges2   ', 'nEdges      '], data_type='integer', dataset_name='mpas_mesh')
    call io_add_var('verticesOnCell',    long_name='Vertex indices that surround a given cell',                       units='',       dim_names=['maxEdges    ', 'nCells      '], data_type='integer', dataset_name='mpas_mesh')
    call io_add_var('verticesOnEdge',    long_name='Vertex indices that saddle a given edge',                         units='',       dim_names=['TWO         ', 'nEdges      '], data_type='integer', dataset_name='mpas_mesh')
    call io_add_var('edgesOnVertex',     long_name='Edge indices that radiate from a given vertex',                   units='',       dim_names=['vertexDegree', 'nVertices   '], data_type='integer', dataset_name='mpas_mesh')
    call io_add_var('cellsOnVertex',     long_name='Cell indices that radiate from a given vertex',                   units='',       dim_names=['vertexDegree', 'nVertices   '], data_type='integer', dataset_name='mpas_mesh')
    call io_add_var('weightsOnEdge',     long_name='Weights used to reconstruct tangential velocities',               units='',       dim_names=['maxEdges2   ', 'nEdges      '], data_type='real(8)', dataset_name='mpas_mesh')
    call io_add_var('dvEdge',            long_name='Distance between the vertices that saddle a given edge',          units='m',      dim_names=['nEdges      '],                 data_type='real(8)', dataset_name='mpas_mesh')
    call io_add_var('dv1Edge',           long_name='Distance between the first vertex and edge location',             units='m',      dim_names=['nEdges      '],                 data_type='real(8)', dataset_name='mpas_mesh')
    call io_add_var('dv2Edge',           long_name='Distance between the second vertex and edge location',            units='m',      dim_names=['nEdges      '],                 data_type='real(8)', dataset_name='mpas_mesh')
    call io_add_var('dcEdge',            long_name='Distance between the cells that saddle a given edge',             units='m',      dim_names=['nEdges      '],                 data_type='real(8)', dataset_name='mpas_mesh')
    call io_add_var('angleEdge',         long_name='Angle between edge normal vector and local eastward direction',   units='radian', dim_names=['nEdges      '],                 data_type='real(8)', dataset_name='mpas_mesh')
    call io_add_var('areaCell',          long_name='Primary cell area',                                               units='m2',     dim_names=['nCells      '],                 data_type='real(8)', dataset_name='mpas_mesh')
    call io_add_var('areaTriangle',      long_name='Dual cell area',                                                  units='m2',     dim_names=['nVertices   '],                 data_type='real(8)', dataset_name='mpas_mesh')
    call io_add_var('kiteAreasOnVertex', long_name='Intersection area between areaCell and areaTriangle',             units='m2',     dim_names=['vertexDegree', 'nVertices   '], data_type='real(8)', dataset_name='mpas_mesh')
    call io_start_output(dataset_name='mpas_mesh')
    call io_output('latCell',           mesh%latCell,           dataset_name='mpas_mesh')
    call io_output('lonCell',           mesh%lonCell,           dataset_name='mpas_mesh')
    call io_output('xCell',             mesh%xCell,             dataset_name='mpas_mesh')
    call io_output('yCell',             mesh%yCell,             dataset_name='mpas_mesh')
    call io_output('zCell',             mesh%zCell,             dataset_name='mpas_mesh')
    call io_output('indexToCellID',     mesh%indexToCellID,     dataset_name='mpas_mesh')
    call io_output('lonEdge',           mesh%lonEdge,           dataset_name='mpas_mesh')
    call io_output('latEdge',           mesh%latEdge,           dataset_name='mpas_mesh')
    call io_output('xEdge',             mesh%xEdge,             dataset_name='mpas_mesh')
    call io_output('yEdge',             mesh%yEdge,             dataset_name='mpas_mesh')
    call io_output('zEdge',             mesh%zEdge,             dataset_name='mpas_mesh')
    call io_output('indexToEdgeID',     mesh%indexToEdgeID,     dataset_name='mpas_mesh')
    call io_output('lonVertex',         mesh%lonVertex,         dataset_name='mpas_mesh')
    call io_output('latVertex',         mesh%latVertex,         dataset_name='mpas_mesh')
    call io_output('xVertex',           mesh%xVertex,           dataset_name='mpas_mesh')
    call io_output('yVertex',           mesh%yVertex,           dataset_name='mpas_mesh')
    call io_output('zVertex',           mesh%zVertex,           dataset_name='mpas_mesh')
    call io_output('indexToVertexID',   mesh%indexToVertexID,   dataset_name='mpas_mesh')
    call io_output('nEdgesOnCell',      mesh%nEdgesOnCell,      dataset_name='mpas_mesh')
    call io_output('nEdgesOnEdge',      mesh%nEdgesOnEdge,      dataset_name='mpas_mesh')
    call io_output('cellsOnCell',       mesh%cellsOnCell,       dataset_name='mpas_mesh')
    call io_output('cellsOnEdge',       mesh%cellsOnEdge,       dataset_name='mpas_mesh')
    call io_output('edgesOnCell',       mesh%edgesOnCell,       dataset_name='mpas_mesh')
    call io_output('edgesOnEdge',       mesh%edgesOnEdge,       dataset_name='mpas_mesh')
    call io_output('verticesOnCell',    mesh%verticesOnCell,    dataset_name='mpas_mesh')
    call io_output('verticesOnEdge',    mesh%verticesOnEdge,    dataset_name='mpas_mesh')
    call io_output('edgesOnVertex',     mesh%edgesOnVertex,     dataset_name='mpas_mesh')
    call io_output('cellsOnVertex',     mesh%cellsOnVertex,     dataset_name='mpas_mesh')
    call io_output('weightsOnEdge',     mesh%weightsOnEdge,     dataset_name='mpas_mesh')
    call io_output('dvEdge',            mesh%dvEdge,            dataset_name='mpas_mesh')
    call io_output('dv1Edge',           mesh%dv1Edge,           dataset_name='mpas_mesh')
    call io_output('dv2Edge',           mesh%dv2Edge,           dataset_name='mpas_mesh')
    call io_output('dcEdge',            mesh%dcEdge,            dataset_name='mpas_mesh')
    call io_output('angleEdge',         mesh%angleEdge,         dataset_name='mpas_mesh')
    call io_output('areaCell',          mesh%areaCell,          dataset_name='mpas_mesh')
    call io_output('areaTriangle',      mesh%areaTriangle,      dataset_name='mpas_mesh')
    call io_output('kiteAreasOnVertex', mesh%kiteAreasOnVertex, dataset_name='mpas_mesh')
    call io_end_output(dataset_name='mpas_mesh')

  end subroutine mpas_mesh_output

  subroutine mpas_mesh_finalize(this)

    type(mpas_mesh_type), intent(inout) :: this

    if (allocated(this%lonCell))           deallocate(this%lonCell)
    if (allocated(this%latCell))           deallocate(this%latCell)
    if (allocated(this%xCell))             deallocate(this%xCell)
    if (allocated(this%yCell))             deallocate(this%yCell)
    if (allocated(this%zCell))             deallocate(this%zCell)
    if (allocated(this%indexToCellID))     deallocate(this%indexToCellID)
    if (allocated(this%lonEdge))           deallocate(this%lonEdge)
    if (allocated(this%latEdge))           deallocate(this%latEdge)
    if (allocated(this%xEdge))             deallocate(this%xEdge)
    if (allocated(this%yEdge))             deallocate(this%yEdge)
    if (allocated(this%zEdge))             deallocate(this%zEdge)
    if (allocated(this%indexToEdgeID))     deallocate(this%indexToEdgeID)
    if (allocated(this%lonVertex))         deallocate(this%lonVertex)
    if (allocated(this%latVertex))         deallocate(this%latVertex)
    if (allocated(this%xVertex))           deallocate(this%xVertex)
    if (allocated(this%yVertex))           deallocate(this%yVertex)
    if (allocated(this%zVertex))           deallocate(this%zVertex)
    if (allocated(this%indexToVertexID))   deallocate(this%indexToVertexID)
    if (allocated(this%nEdgesOnCell))      deallocate(this%nEdgesOnCell)
    if (allocated(this%nEdgesOnEdge))      deallocate(this%nEdgesOnEdge)
    if (allocated(this%cellsOnCell))       deallocate(this%cellsOnCell)
    if (allocated(this%cellsOnEdge))       deallocate(this%cellsOnEdge)
    if (allocated(this%edgesOnCell))       deallocate(this%edgesOnCell)
    if (allocated(this%edgesOnEdge))       deallocate(this%edgesOnEdge)
    if (allocated(this%verticesOnCell))    deallocate(this%verticesOnCell)
    if (allocated(this%verticesOnEdge))    deallocate(this%verticesOnEdge)
    if (allocated(this%edgesOnVertex))     deallocate(this%edgesOnVertex)
    if (allocated(this%cellsOnVertex))     deallocate(this%cellsOnVertex)
    if (allocated(this%weightsOnEdge))     deallocate(this%weightsOnEdge)
    if (allocated(this%dvEdge))            deallocate(this%dvEdge)
    if (allocated(this%dv1Edge))           deallocate(this%dv1Edge)
    if (allocated(this%dv2Edge))           deallocate(this%dv2Edge)
    if (allocated(this%dcEdge))            deallocate(this%dcEdge)
    if (allocated(this%angleEdge))         deallocate(this%angleEdge)
    if (allocated(this%areaCell))          deallocate(this%areaCell)
    if (allocated(this%areaTriangle))      deallocate(this%areaTriangle)
    if (allocated(this%kiteAreasOnVertex)) deallocate(this%kiteAreasOnVertex)


  end subroutine mpas_mesh_finalize

  subroutine calc_weightsOnEdge(VVT_array, VC_array, VE_array, mesh)

    type(array_type), intent(in) :: VVT_array
    type(array_type), intent(in) :: VC_array
    type(array_type), intent(in) :: VE_array
    type(mpas_mesh_type), intent(inout) :: mesh

    real(8), allocatable :: R(:,:)
    integer, allocatable :: n(:,:), t(:,:)
    integer iCell, iEdge, iEdgeOnEdge, iVertex, i, j, i0
    integer iLocalCell, iLocalVertex
    integer nEdgesOnCell
    integer iLocalEdge1, iLocalEdge2, iLocalEdge, iLocalEdgeOnVertex

    ! Calculate weightsOnEdge for reconstructing tangential velocities.
    allocate(R(maxEdges,mesh%nCells))
    allocate(n(maxEdges,mesh%nCells))
    allocate(t(vertexDegree,mesh%nVertices))
    ! Calculate divergence interpolation weights.
    R = 0.0
    do iVertex = 1, mesh%nVertices
      do iLocalCell = 1, vertexDegree
        iCell = mesh%cellsOnVertex(iLocalCell,iVertex)
        do iLocalVertex = 1, mesh%nEdgesOnCell(iCell)
          if (mesh%verticesOnCell(iLocalVertex,iCell) == iVertex) exit
        end do
        R(iLocalVertex,iCell) = R(iLocalVertex,iCell) + mesh%kiteAreasOnVertex(iLocalCell,iVertex)
      end do
    end do
    do iCell = 1, mesh%nCells
      nEdgesOnCell = mesh%nEdgesOnCell(iCell)
      R(:nEdgesOnCell,iCell) = R(:nEdgesOnCell,iCell) / sum(R(:nEdgesOnCell,iCell))
    end do
    ! Set normal vector indicator.
    n = 0
    do iCell = 1, mesh%nCells
      do i = 1, mesh%nEdgesOnCell(iCell)
        if (iCell == mesh%cellsOnEdge(1,mesh%edgesOnCell(i,iCell))) then
          n(i,iCell) =  1
        else if (iCell == mesh%cellsOnEdge(2,mesh%edgesOnCell(i,iCell))) then
          n(i,iCell) = -1
        end if
      end do
    end do
    ! Set tangential vector indicator.
    t = 0
    do iVertex = 1, mesh%nVertices
      do i = 1, vertexDegree
        if (iVertex == mesh%verticesOnEdge(1,mesh%edgesOnVertex(i,iVertex))) then
          t(i,iVertex) =  1
        else if (iVertex == mesh%verticesOnEdge(2,mesh%edgesOnVertex(i,iVertex))) then
          t(i,iVertex) = -1
        end if
      end do
    end do
    do iEdge = 1, mesh%nEdges
      ! Get the local index of edge on both side cells.
      do iLocalEdge = 1, mesh%nEdgesOnCell(mesh%cellsOnEdge(1,iEdge))
        if (mesh%edgesOnCell(iLocalEdge,mesh%cellsOnEdge(1,iEdge)) == iEdge) then
          iLocalEdge1 = iLocalEdge
          exit
        end if
      end do
      do iLocalEdge = 1, mesh%nEdgesOnCell(mesh%cellsOnEdge(2,iEdge))
        if (mesh%edgesOnCell(iLocalEdge,mesh%cellsOnEdge(2,iEdge)) == iEdge) then
          iLocalEdge2 = iLocalEdge
          exit
        end if
      end do
      ! Loop through all the involved edges.
      do i = 1, mesh%nEdgesOnEdge(iEdge)
        iEdgeOnEdge = mesh%edgesOnEdge(i,iEdge)
        if (any(mesh%cellsOnEdge(:,iEdgeOnEdge) == mesh%cellsOnEdge(1,iEdge))) then
          ! For the first cell
          i0 = 1
          iCell = mesh%cellsOnEdge(1,iEdge)
          iLocalEdge = iLocalEdge1
        else
          ! For the second cell
          if (i0 == 1) i0 = i
          iCell = mesh%cellsOnEdge(2,iEdge)
          iLocalEdge = iLocalEdge2
        end if
        iLocalVertex = iLocalEdge
        do j = i0, i
          iLocalVertex = iLocalVertex + 1
          if (iLocalVertex > mesh%nEdgesOnCell(iCell)) iLocalVertex = 1
          iVertex = mesh%verticesOnCell(iLocalVertex,iCell)
          mesh%weightsOnEdge(i,iEdge) = mesh%weightsOnEdge(i,iEdge) + R(iLocalVertex,iCell)
        end do
        do iLocalEdgeOnVertex = 1, vertexDegree
          if (mesh%edgesOnVertex(iLocalEdgeOnVertex,iVertex) == iEdgeOnEdge) exit
        end do
        mesh%weightsOnEdge(i,iEdge) = n(iLocalEdge,iCell) / t(iLocalEdgeOnVertex,iVertex) * (mesh%weightsOnEdge(i,iEdge) - 0.5)
        mesh%weightsOnEdge(i,iEdge) = mesh%weightsOnEdge(i,iEdge) * mesh%dcEdge(iEdge) / mesh%dvEdge(iEdgeOnEdge)
      end do
    end do
    deallocate(R)
    deallocate(n)
    deallocate(t)

  end subroutine calc_weightsOnEdge

end module mpas_mesh_mod
