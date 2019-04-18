program gen_mesh

  use const_mod
  use log_mod
  use io_mod
  use string_mod
  use array_mod
  use linked_list_mod
  use delaunay_voronoi_mod
  use sphere_geometry_mod

  implicit none

  integer, parameter :: nx = 40962
  real(8), allocatable :: x(:,:)

  allocate(x(nx,3))
  
  call io_init()
  call io_create_dataset(name='grids', file_path='./point.' // trim(to_string(nx)) // '.nc', mode='input')
  call io_start_input('grids')
  call io_input('vtx_p', x, dataset_name='grids')
  call io_end_input('grids')

  call delaunay_voronoi_init(nx, x=x(:,1), y=x(:,2), z=x(:,3))
  call delaunay_triangulation()
  call voronoi_diagram()
  call delaunay_voronoi_output(output)

  deallocate(x)

contains

  subroutine output(DVT_array, DT_list, VVT_array, VC_array, tag)

    type(array_type), intent(in) :: DVT_array
    type(linked_list_type), intent(in) :: DT_list
    type(array_type), intent(in) :: VVT_array
    type(array_type), intent(in) :: VC_array
    character(*), intent(in), optional :: tag

    integer i, j, k, num_DE, num_VE
    real(8), allocatable :: lon_DVT(:)
    real(8), allocatable :: lat_DVT(:)
    real(8), allocatable :: lon_VVT(:)
    real(8), allocatable :: lat_VVT(:)
    integer, allocatable :: DT_DVT_idx(:,:)
    integer, allocatable :: DE_DVT_idx(:,:)
    integer, allocatable :: VC_VVT_idx(:,:)
    integer, allocatable :: VE_VVT_idx(:,:)
    type(linked_list_iterator_type) iterator
    type(delaunay_vertex_type), pointer :: DVT, DVT1, DVT2
    type(delaunay_triangle_type), pointer :: adjDT
    type(voronoi_vertex_type), pointer :: VVT, VVT1, VVT2
    type(voronoi_cell_type), pointer :: adjVC, VC

    num_DE = euler_formula(num_cell=DT_list%size, num_vertex=DVT_array%size)
    num_VE = euler_formula(num_cell=VC_array%size, num_vertex=VVT_array%size)

    allocate(lon_DVT(DVT_array%size))
    allocate(lat_DVT(DVT_array%size))
    allocate(lon_VVT(VVT_array%size))
    allocate(lat_VVT(VVT_array%size))
    allocate(DT_DVT_idx(3,DT_list%size))
    allocate(DE_DVT_idx(2,num_DE))
    allocate(VC_VVT_idx(6,VC_array%size))
    allocate(VE_VVT_idx(2,num_VE))

    do i = 1, DVT_array%size
      select type (DVT => DVT_array%value_at(i))
      type is (delaunay_vertex_type)
        lon_DVT(i) = DVT%lon * deg
        lat_DVT(i) = DVT%lat * deg
      end select
    end do

    do i = 1, VVT_array%size
      select type (VVT => VVT_array%value_at(i))
      type is (voronoi_vertex_type)
        lon_VVT(i) = VVT%lon * deg
        lat_VVT(i) = VVT%lat * deg
      end select
    end do

    i = 1
    j = 1
    iterator = linked_list_iterator(DT_list)
    do while (.not. iterator%ended())
      select type (DT => iterator%value)
      type is (delaunay_triangle_type)
        do k = 1, 3
          DVT => get_DVT(DT%DVT, k)
          DT_DVT_idx(k,i) = DVT%id
          adjDT => get_DT(DT%adjDT, k)
          if (adjDT%id > DT%id) then
            DVT1 => get_DVT(DT%DVT, ip1(k))
            DVT2 => get_DVT(DT%DVT, im1(k))
            DE_DVT_idx(1,j) = DVT1%id
            DE_DVT_idx(2,j) = DVT2%id
            j = j + 1 
          end if
        end do
      end select
      i = i + 1
      call iterator%next()
    end do

    j = 1
    do i = 1, VC_array%size
      VC => get_VC(VC_array, i)
      do k = 1, VC%VVT%size
        VVT => get_VVT(VC%VVT, k)
        VC_VVT_idx(k,i) = VVT%id
        adjVC => get_VC(VC%adjVC, k)
        if (adjVC%id > VC%id) then
          VVT1 => get_VVT(VC%VVT, k)
          VVT2 => get_VVT(VC%VVT, merge(1, k + 1, k == VC%VVT%size))
          VE_VVT_idx(1,j) = VVT1%id
          VE_VVT_idx(2,j) = VVT2%id
          j = j + 1 
        end if
      end do
    end do

    call io_create_dataset(name='delaunay', file_path='mesh.' // trim(to_string(nx)) // '.nc', mode='output')
    call io_add_dim('num_DVT', size=DVT_array%size, dataset_name='delaunay')
    call io_add_dim('num_DT',  size=DT_list%size,   dataset_name='delaunay')
    call io_add_dim('num_DE',  size=num_DE,         dataset_name='delaunay')
    call io_add_dim('num_VVT', size=VVT_array%size, dataset_name='delaunay')
    call io_add_dim('num_VC',  size=VC_array%size,  dataset_name='delaunay')
    call io_add_dim('num_VE',  size=num_DE,         dataset_name='delaunay')
    call io_add_dim('TWO',     size=2,              dataset_name='delaunay')
    call io_add_dim('THREE',   size=3,              dataset_name='delaunay')
    call io_add_dim('SIX',     size=6,              dataset_name='delaunay')
    call io_add_var('lon_DVT',    long_name='longitude of DVT',  units='degree_east',  dim_names=['num_DVT'],            data_type='real(8)', dataset_name='delaunay')
    call io_add_var('lat_DVT',    long_name='latitude of DVT',   units='degree_north', dim_names=['num_DVT'],            data_type='real(8)', dataset_name='delaunay')
    call io_add_var('lon_VVT',    long_name='longitude of VVT',  units='degree_east',  dim_names=['num_VVT'],            data_type='real(8)', dataset_name='delaunay')
    call io_add_var('lat_VVT',    long_name='latitude of VVT',   units='degree_north', dim_names=['num_VVT'],            data_type='real(8)', dataset_name='delaunay')
    call io_add_var('DT_DVT_idx', long_name='DVT indices of DT', units='1',            dim_names=['THREE  ', 'num_DT '], data_type='integer', dataset_name='delaunay')
    call io_add_var('DE_DVT_idx', long_name='DVT indices of DE', units='1',            dim_names=['TWO    ', 'num_DE '], data_type='integer', dataset_name='delaunay')
    call io_add_var('VC_VVT_idx', long_name='VVT indices of VC', units='1',            dim_names=['SIX    ', 'num_VC '], data_type='integer', dataset_name='delaunay')
    call io_add_var('VE_VVT_idx', long_name='VVT indices of VE', units='1',            dim_names=['TWO    ', 'num_VE '], data_type='integer', dataset_name='delaunay')
    call io_start_output(dataset_name='delaunay', tag=tag)
    call io_output('lon_DVT',    lon_DVT, dataset_name='delaunay')
    call io_output('lat_DVT',    lat_DVT, dataset_name='delaunay')
    call io_output('lon_VVT',    lon_VVT, dataset_name='delaunay')
    call io_output('lat_VVT',    lat_VVT, dataset_name='delaunay')
    call io_output('DT_DVT_idx', DT_DVT_idx, dataset_name='delaunay')
    call io_output('DE_DVT_idx', DE_DVT_idx, dataset_name='delaunay')
    call io_output('VC_VVT_idx', VC_VVT_idx, dataset_name='delaunay')
    call io_output('VE_VVT_idx', VE_VVT_idx, dataset_name='delaunay')
    call io_end_output('delaunay')

  end subroutine output

end program gen_mesh
