module mpas_mesh_mod

  use flogger
  use fiona
  use container
  use string
  use const_mod
  use delaunay_voronoi_mod
  use sphere_geometry_mod

  implicit none

  private

  public mpas_mesh_output
  public mpas_mesh_type

  integer, parameter :: vertexDegree = 3
  integer, parameter :: maxEdges = 6
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
    procedure :: init => mpas_mesh_init
    final :: mpas_mesh_finalize
  end type mpas_mesh_type

contains

  subroutine mpas_mesh_output(DVT_array, DT_list, DE_array, VVT_array, VC_array, VE_array, tag)

    type(array_type), intent(in) :: DVT_array
    type(linked_list_type), intent(in) :: DT_list
    type(array_type), intent(in) :: DE_array
    type(array_type), intent(in) :: VVT_array
    type(array_type), intent(in) :: VC_array
    type(array_type), intent(in) :: VE_array
    character(*), intent(in), optional :: tag

    type(mpas_mesh_type) mesh

    call mpas_mesh_adapt(DVT_array, DT_list, DE_array, VVT_array, VC_array, VE_array, mesh)

    call fiona_create_dataset('mpas_mesh', file_path='mpas_mesh.' // trim(to_string(mesh%nCells)) // '.nc')
    call fiona_add_att('mpas_mesh', 'on_a_sphere',   'YES')
    call fiona_add_att('mpas_mesh', 'sphere_radius', 1.0d0)
    call fiona_add_att('mpas_mesh', 'mesh_id',       'demo')
    call fiona_add_att('mpas_mesh', 'mesh',          '1.0')
    call fiona_add_dim('mpas_mesh', 'nCells',       size=mesh%nCells)
    call fiona_add_dim('mpas_mesh', 'nEdges',       size=mesh%nEdges)
    call fiona_add_dim('mpas_mesh', 'nVertices',    size=mesh%nVertices)
    call fiona_add_dim('mpas_mesh', 'vertexDegree', size=vertexDegree)
    call fiona_add_dim('mpas_mesh', 'maxEdges',     size=maxEdges)
    call fiona_add_dim('mpas_mesh', 'maxEdges2',    size=maxEdges2)
    call fiona_add_dim('mpas_mesh', 'TWO',          size=2)
    call fiona_add_dim('mpas_mesh', 'Time')
    call fiona_add_var('mpas_mesh', 'latCell',           long_name='Latitude of all cell centers',                                    units='radian', dim_names=['nCells      '],                 data_type='real(8)')
    call fiona_add_var('mpas_mesh', 'lonCell',           long_name='Longitude of all cell centers',                                   units='radian', dim_names=['nCells      '],                 data_type='real(8)')
    call fiona_add_var('mpas_mesh', 'xCell',             long_name='x axis position of all cell centers',                             units='m',      dim_names=['nCells      '],                 data_type='real(8)')
    call fiona_add_var('mpas_mesh', 'yCell',             long_name='y axis position of all cell centers',                             units='m',      dim_names=['nCells      '],                 data_type='real(8)')
    call fiona_add_var('mpas_mesh', 'zCell',             long_name='z axis position of all cell centers',                             units='m',      dim_names=['nCells      '],                 data_type='real(8)')
    call fiona_add_var('mpas_mesh', 'indexToCellID',     long_name='Global cell ID for all cell centers',                             units='',       dim_names=['nCells      '],                 data_type='integer')
    call fiona_add_var('mpas_mesh', 'latEdge',           long_name='Latitude of all edge locations',                                  units='radian', dim_names=['nEdges      '],                 data_type='real(8)')
    call fiona_add_var('mpas_mesh', 'lonEdge',           long_name='Longitude of all edge locations',                                 units='radian', dim_names=['nEdges      '],                 data_type='real(8)')
    call fiona_add_var('mpas_mesh', 'xEdge',             long_name='x axis position of all edge locations',                           units='m',      dim_names=['nEdges      '],                 data_type='real(8)')
    call fiona_add_var('mpas_mesh', 'yEdge',             long_name='y axis position of all edge locations',                           units='m',      dim_names=['nEdges      '],                 data_type='real(8)')
    call fiona_add_var('mpas_mesh', 'zEdge',             long_name='z axis position of all edge locations',                           units='m',      dim_names=['nEdges      '],                 data_type='real(8)')
    call fiona_add_var('mpas_mesh', 'indexToEdgeID',     long_name='Global edge ID for all edge locations',                           units='',       dim_names=['nEdges      '],                 data_type='integer')
    call fiona_add_var('mpas_mesh', 'latVertex',         long_name='Latitude of all cell vertices',                                   units='radian', dim_names=['nVertices   '],                 data_type='real(8)')
    call fiona_add_var('mpas_mesh', 'lonVertex',         long_name='Longitude of all cell vertices',                                  units='radian', dim_names=['nVertices   '],                 data_type='real(8)')
    call fiona_add_var('mpas_mesh', 'xVertex',           long_name='x axis position of all cell vertices',                            units='m',      dim_names=['nVertices   '],                 data_type='real(8)')
    call fiona_add_var('mpas_mesh', 'yVertex',           long_name='y axis position of all cell vertices',                            units='m',      dim_names=['nVertices   '],                 data_type='real(8)')
    call fiona_add_var('mpas_mesh', 'zVertex',           long_name='z axis position of all cell vertices',                            units='m',      dim_names=['nVertices   '],                 data_type='real(8)')
    call fiona_add_var('mpas_mesh', 'indexToVertexID',   long_name='Global vertex ID for all cell vertices',                          units='',       dim_names=['nVertices   '],                 data_type='integer')
    call fiona_add_var('mpas_mesh', 'nEdgesOnCell',      long_name='Number of edges on a given cell',                                 units='',       dim_names=['nCells      '],                 data_type='integer')
    call fiona_add_var('mpas_mesh', 'nEdgesOnEdge',      long_name='Number of edges on a given edge',                                 units='',       dim_names=['nEdges      '],                 data_type='integer')
    call fiona_add_var('mpas_mesh', 'cellsOnCell',       long_name='Cell indices that surround a given cell',                         units='',       dim_names=['maxEdges    ', 'nCells      '], data_type='integer')
    call fiona_add_var('mpas_mesh', 'cellsOnEdge',       long_name='Cell indices that saddle a given edge',                           units='',       dim_names=['TWO         ', 'nEdges      '], data_type='integer')
    call fiona_add_var('mpas_mesh', 'edgesOnCell',       long_name='Edge indices that surround a given cell',                         units='',       dim_names=['maxEdges    ', 'nCells      '], data_type='integer')
    call fiona_add_var('mpas_mesh', 'edgesOnEdge',       long_name='Edge indices that are used to reconstruct tangential velocities', units='',       dim_names=['maxEdges2   ', 'nEdges      '], data_type='integer')
    call fiona_add_var('mpas_mesh', 'verticesOnCell',    long_name='Vertex indices that surround a given cell',                       units='',       dim_names=['maxEdges    ', 'nCells      '], data_type='integer')
    call fiona_add_var('mpas_mesh', 'verticesOnEdge',    long_name='Vertex indices that saddle a given edge',                         units='',       dim_names=['TWO         ', 'nEdges      '], data_type='integer')
    call fiona_add_var('mpas_mesh', 'edgesOnVertex',     long_name='Edge indices that radiate from a given vertex',                   units='',       dim_names=['vertexDegree', 'nVertices   '], data_type='integer')
    call fiona_add_var('mpas_mesh', 'cellsOnVertex',     long_name='Cell indices that radiate from a given vertex',                   units='',       dim_names=['vertexDegree', 'nVertices   '], data_type='integer')
    call fiona_add_var('mpas_mesh', 'weightsOnEdge',     long_name='Weights used to reconstruct tangential velocities',               units='',       dim_names=['maxEdges2   ', 'nEdges      '], data_type='real(8)')
    call fiona_add_var('mpas_mesh', 'dvEdge',            long_name='Distance between the vertices that saddle a given edge',          units='m',      dim_names=['nEdges      '],                 data_type='real(8)')
    call fiona_add_var('mpas_mesh', 'dv1Edge',           long_name='Distance between the first vertex and edge location',             units='m',      dim_names=['nEdges      '],                 data_type='real(8)')
    call fiona_add_var('mpas_mesh', 'dv2Edge',           long_name='Distance between the second vertex and edge location',            units='m',      dim_names=['nEdges      '],                 data_type='real(8)')
    call fiona_add_var('mpas_mesh', 'dcEdge',            long_name='Distance between the cells that saddle a given edge',             units='m',      dim_names=['nEdges      '],                 data_type='real(8)')
    call fiona_add_var('mpas_mesh', 'angleEdge',         long_name='Angle between edge normal vector and local eastward direction',   units='radian', dim_names=['nEdges      '],                 data_type='real(8)')
    call fiona_add_var('mpas_mesh', 'areaCell',          long_name='Primary cell area',                                               units='m2',     dim_names=['nCells      '],                 data_type='real(8)')
    call fiona_add_var('mpas_mesh', 'areaTriangle',      long_name='Dual cell area',                                                  units='m2',     dim_names=['nVertices   '],                 data_type='real(8)')
    call fiona_add_var('mpas_mesh', 'kiteAreasOnVertex', long_name='Intersection area between areaCell and areaTriangle',             units='m2',     dim_names=['vertexDegree', 'nVertices   '], data_type='real(8)')
    call fiona_start_output('mpas_mesh')
    call fiona_output('mpas_mesh', 'latCell',           mesh%latCell)
    call fiona_output('mpas_mesh', 'lonCell',           mesh%lonCell)
    call fiona_output('mpas_mesh', 'xCell',             mesh%xCell)
    call fiona_output('mpas_mesh', 'yCell',             mesh%yCell)
    call fiona_output('mpas_mesh', 'zCell',             mesh%zCell)
    call fiona_output('mpas_mesh', 'indexToCellID',     mesh%indexToCellID)
    call fiona_output('mpas_mesh', 'lonEdge',           mesh%lonEdge)
    call fiona_output('mpas_mesh', 'latEdge',           mesh%latEdge)
    call fiona_output('mpas_mesh', 'xEdge',             mesh%xEdge)
    call fiona_output('mpas_mesh', 'yEdge',             mesh%yEdge)
    call fiona_output('mpas_mesh', 'zEdge',             mesh%zEdge)
    call fiona_output('mpas_mesh', 'indexToEdgeID',     mesh%indexToEdgeID)
    call fiona_output('mpas_mesh', 'lonVertex',         mesh%lonVertex)
    call fiona_output('mpas_mesh', 'latVertex',         mesh%latVertex)
    call fiona_output('mpas_mesh', 'xVertex',           mesh%xVertex)
    call fiona_output('mpas_mesh', 'yVertex',           mesh%yVertex)
    call fiona_output('mpas_mesh', 'zVertex',           mesh%zVertex)
    call fiona_output('mpas_mesh', 'indexToVertexID',   mesh%indexToVertexID)
    call fiona_output('mpas_mesh', 'nEdgesOnCell',      mesh%nEdgesOnCell)
    call fiona_output('mpas_mesh', 'nEdgesOnEdge',      mesh%nEdgesOnEdge)
    call fiona_output('mpas_mesh', 'cellsOnCell',       mesh%cellsOnCell)
    call fiona_output('mpas_mesh', 'cellsOnEdge',       mesh%cellsOnEdge)
    call fiona_output('mpas_mesh', 'edgesOnCell',       mesh%edgesOnCell)
    call fiona_output('mpas_mesh', 'edgesOnEdge',       mesh%edgesOnEdge)
    call fiona_output('mpas_mesh', 'verticesOnCell',    mesh%verticesOnCell)
    call fiona_output('mpas_mesh', 'verticesOnEdge',    mesh%verticesOnEdge)
    call fiona_output('mpas_mesh', 'edgesOnVertex',     mesh%edgesOnVertex)
    call fiona_output('mpas_mesh', 'cellsOnVertex',     mesh%cellsOnVertex)
    call fiona_output('mpas_mesh', 'weightsOnEdge',     mesh%weightsOnEdge)
    call fiona_output('mpas_mesh', 'dvEdge',            mesh%dvEdge)
    call fiona_output('mpas_mesh', 'dv1Edge',           mesh%dv1Edge)
    call fiona_output('mpas_mesh', 'dv2Edge',           mesh%dv2Edge)
    call fiona_output('mpas_mesh', 'dcEdge',            mesh%dcEdge)
    call fiona_output('mpas_mesh', 'angleEdge',         mesh%angleEdge)
    call fiona_output('mpas_mesh', 'areaCell',          mesh%areaCell)
    call fiona_output('mpas_mesh', 'areaTriangle',      mesh%areaTriangle)
    call fiona_output('mpas_mesh', 'kiteAreasOnVertex', mesh%kiteAreasOnVertex)
    call fiona_end_output(dataset_name='mpas_mesh')

  end subroutine mpas_mesh_output

  subroutine mpas_mesh_adapt(DVT_array, DT_list, DE_array, VVT_array, VC_array, VE_array, mesh)

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

    call mesh%init(VC_array%size, VE_array%size, VVT_array%size)

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
      do i = 1, VVT%VC%size
        VC => get_VC(VVT%VC, i)
        mesh%cellsOnVertex(i,iVertex) = VC%id
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

    call calc_weightsOnEdge(mesh)

  end subroutine mpas_mesh_adapt

  subroutine mpas_mesh_init(this, nCells, nEdges, nVertices)

    class(mpas_mesh_type), intent(inout) :: this
    integer, intent(in) :: nCells
    integer, intent(in) :: nEdges
    integer, intent(in) :: nVertices

    this%nCells = nCells
    this%nEdges = nEdges
    this%nVertices = nVertices

    allocate(this%latCell          (this%nCells))
    allocate(this%lonCell          (this%nCells))
    allocate(this%xCell            (this%nCells))
    allocate(this%yCell            (this%nCells))
    allocate(this%zCell            (this%nCells))
    allocate(this%indexToCellID    (this%nCells))
    allocate(this%lonEdge          (this%nEdges))
    allocate(this%latEdge          (this%nEdges))
    allocate(this%xEdge            (this%nEdges))
    allocate(this%yEdge            (this%nEdges))
    allocate(this%zEdge            (this%nEdges))
    allocate(this%indexToEdgeID    (this%nEdges))
    allocate(this%lonVertex        (this%nVertices))
    allocate(this%latVertex        (this%nVertices))
    allocate(this%xVertex          (this%nVertices))
    allocate(this%yVertex          (this%nVertices))
    allocate(this%zVertex          (this%nVertices))
    allocate(this%indexToVertexID  (this%nVertices))
    allocate(this%nEdgesOnCell     (this%nCells))
    allocate(this%nEdgesOnEdge     (this%nEdges))
    allocate(this%cellsOnCell      (maxEdges,this%nCells))
    allocate(this%cellsOnEdge      (2,this%nEdges))
    allocate(this%edgesOnCell      (maxEdges,this%nCells))
    allocate(this%edgesOnEdge      (maxEdges2,this%nEdges))
    allocate(this%verticesOnCell   (maxEdges,this%nCells))
    allocate(this%verticesOnEdge   (2,this%nEdges))
    allocate(this%edgesOnVertex    (vertexDegree,this%nVertices))
    allocate(this%cellsOnVertex    (vertexDegree,this%nVertices))
    allocate(this%weightsOnEdge    (maxEdges2,this%nEdges))
    allocate(this%dvEdge           (this%nEdges))
    allocate(this%dv1Edge          (this%nEdges))
    allocate(this%dv2Edge          (this%nEdges))
    allocate(this%dcEdge           (this%nEdges))
    allocate(this%angleEdge        (this%nEdges))
    allocate(this%areaCell         (this%nCells))
    allocate(this%areaTriangle     (this%nVertices))
    allocate(this%kiteAreasOnVertex(vertexDegree,this%nVertices))

  end subroutine mpas_mesh_init

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

  subroutine calc_weightsOnEdge(mesh)

    type(mpas_mesh_type), intent(inout) :: mesh

    real(8), allocatable :: R(:,:)
    integer, allocatable :: n(:,:), t(:,:)
    integer iCell1, iCell2, iCell, iEdge, iEdgeOnEdge, iVertex, i, j
    integer iLocalCell, iLocalVertex
    integer nEdgesOnCell1, nEdgesOnCell2, nEdgesOnCell, jEnd
    integer iLocalEdge1, iLocalEdge2, iLocalEdge, iLocalEdgeOnVertex

    ! Calculate weightsOnEdge for reconstructing tangential velocities.
    allocate(R(maxEdges,mesh%nCells))
    allocate(n(maxEdges,mesh%nCells))
    allocate(t(vertexDegree,mesh%nVertices))
    ! Calculate divergence interpolation weights.
    do iVertex = 1, mesh%nVertices
      do iLocalCell = 1, vertexDegree
        iCell = mesh%cellsOnVertex(iLocalCell,iVertex)
        do iLocalVertex = 1, mesh%nEdgesOnCell(iCell)
          if (mesh%verticesOnCell(iLocalVertex,iCell) == iVertex) exit
        end do
        if (iLocalVertex > mesh%nEdgesOnCell(iCell)) call log_error('Internal error!', __FILE__, __LINE__)
        R(iLocalVertex,iCell) = mesh%kiteAreasOnVertex(iLocalCell,iVertex) / mesh%areaCell(iCell)
      end do
    end do
#ifndef NDEBUG
    do iCell = 1, mesh%nCells
      if (abs(sum(R(:mesh%nEdgesOnCell(iCell),iCell)) - 1.0) > 1.0d-10) then
        call log_error('Internal error!', __FILE__, __LINE__)
      end if
    end do
#endif
    ! Set normal vector indicator.
    do iCell = 1, mesh%nCells
      do iLocalEdge = 1, mesh%nEdgesOnCell(iCell)
        iEdge = mesh%edgesOnCell(iLocalEdge,iCell)
        if (iCell == mesh%cellsOnEdge(1,iEdge)) then
          n(iLocalEdge,iCell) =  1
        else if (iCell == mesh%cellsOnEdge(2,iEdge)) then
          n(iLocalEdge,iCell) = -1
        else
          call log_error('Internal error!', __FILE__, __LINE__)
        end if
      end do
    end do
    ! Set tangential vector indicator.
    do iVertex = 1, mesh%nVertices
      do iLocalEdge = 1, vertexDegree
        iEdge = mesh%edgesOnVertex(iLocalEdge,iVertex)
        if (iVertex == mesh%verticesOnEdge(1,iEdge)) then
          t(iLocalEdge,iVertex) =  1
        else if (iVertex == mesh%verticesOnEdge(2,iEdge)) then
          t(iLocalEdge,iVertex) = -1
        else
          call log_error('Internal error!', __FILE__, __LINE__)
        end if
      end do
    end do
    ! Loop through the target edge.
    do iEdge = 1, mesh%nEdges
      iCell1 = mesh%cellsOnEdge(1,iEdge)
      iCell2 = mesh%cellsOnEdge(2,iEdge)
      nEdgesOnCell1 = mesh%nEdgesOnCell(iCell1)
      nEdgesOnCell2 = mesh%nEdgesOnCell(iCell2)
      ! Get the local index of edge on both side cells.
      iLocalEdge1 = 0
      do iLocalEdge = 1, nEdgesOnCell1
        if (mesh%edgesOnCell(iLocalEdge,iCell1) == iEdge) then
          iLocalEdge1 = iLocalEdge
          exit
        end if
      end do
      if (iLocalEdge1 == 0) call log_error('Internal error!', __FILE__, __LINE__)
      iLocalEdge2 = 0
      do iLocalEdge = 1, nEdgesOnCell2
        if (mesh%edgesOnCell(iLocalEdge,iCell2) == iEdge) then
          iLocalEdge2 = iLocalEdge
          exit
        end if
      end do
      if (iLocalEdge2 == 0) call log_error('Internal error!', __FILE__, __LINE__)
      ! Loop through all the involved edges that contribute to the target edge.
      do i = 1, mesh%nEdgesOnEdge(iEdge)
        iEdgeOnEdge = mesh%edgesOnEdge(i,iEdge)
        if (any(mesh%cellsOnEdge(:,iEdgeOnEdge) == iCell1)) then
          ! For the first cell
          iCell = iCell1
          nEdgesOnCell = nEdgesOnCell1
          jEnd = nEdgesOnCell1 - i
          iLocalEdge = iLocalEdge1 + i
        else if (any(mesh%cellsOnEdge(:,iEdgeOnEdge) == iCell2)) then
          ! For the second cell
          iCell = iCell2
          nEdgesOnCell = nEdgesOnCell2
          jEnd = nEdgesOnCell2 - i + nEdgesOnCell1 - 1
          iLocalEdge = iLocalEdge2 + i - nEdgesOnCell1 + 1
        end if
        if (iLocalEdge > nEdgesOnCell) iLocalEdge = iLocalEdge - nEdgesOnCell
        iLocalVertex = iLocalEdge
        do j = 1, jEnd
          iLocalVertex = iLocalVertex + 1
          if (iLocalVertex > nEdgesOnCell) iLocalVertex = 1
          iVertex = mesh%verticesOnCell(iLocalVertex,iCell)
          mesh%weightsOnEdge(i,iEdge) = mesh%weightsOnEdge(i,iEdge) + R(iLocalVertex,iCell)
        end do
        ! Get the local index of the target edge on the final traversed vertex.
        do iLocalEdgeOnVertex = 1, vertexDegree
          if (mesh%edgesOnVertex(iLocalEdgeOnVertex,iVertex) == iEdge) exit
        end do
        if (iLocalEdgeOnVertex > vertexDegree) call log_error('Internal error!', __FILE__, __LINE__)
        mesh%weightsOnEdge(i,iEdge) = n(iLocalEdge,iCell) / t(iLocalEdgeOnVertex,iVertex) * (mesh%weightsOnEdge(i,iEdge) - 0.5)
        mesh%weightsOnEdge(i,iEdge) = mesh%weightsOnEdge(i,iEdge) * mesh%dvEdge(iEdgeOnEdge) / mesh%dcEdge(iEdge)
      end do
    end do
    deallocate(R)
    deallocate(n)
    deallocate(t)

  end subroutine calc_weightsOnEdge

end module mpas_mesh_mod
