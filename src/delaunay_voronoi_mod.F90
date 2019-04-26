module delaunay_voronoi_mod

  use const_mod
  use log_mod
  use string_mod
  use random_number_mod
  use math_mod
  use sphere_geometry_mod
  use linked_list_mod
  use array_mod

  implicit none

  private

  public delaunay_voronoi_init
  public delaunay_triangulation
  public voronoi_diagram
  public delaunay_voronoi_output
  public delaunay_vertex_type
  public delaunay_triangle_type
  public delaunay_edge_type
  public voronoi_vertex_type
  public voronoi_cell_type
  public voronoi_edge_type
  public get_DVT
  public get_DT
  public get_DE
  public get_VVT
  public get_VC
  public get_VE

  type, extends(point_type) :: delaunay_vertex_type
    integer :: id = -1
    type(linked_list_type) incDT   ! Incident Delaunay triangles
    type(linked_list_type) linkDVT ! Linked Delaunay vertices
    ! For triangulation
    type(delaunay_triangle_type), pointer :: cntDT1    ! Delaunay triangle that contains this vertex.
    type(delaunay_triangle_type), pointer :: cntDT2    ! Delaunay triangle that contains this vertex.
    type(linked_list_item_type), pointer :: stub_item1 ! For fast incident Delaunay triangle deletion.
    type(linked_list_item_type), pointer :: stub_item2 ! For fast incident Delaunay triangle deletion.
    integer edge_idx
  contains
    procedure :: init => delaunay_vertex_init
    procedure :: print => delaunay_vertex_print
  end type delaunay_vertex_type

  type delaunay_triangle_type
    integer :: id = -1
    type(array_type) DVT                         ! Delaunay vertices
    type(array_type) DE                          ! Delaunay edges
    type(array_type) adjDT                       ! Adjacent Delaunay triangles
    type(linked_list_type) subDT                 ! Subdivided Delaunay triangles
    type(linked_list_type) incDVT                ! Included Delaunay vertices
    type(voronoi_vertex_type), pointer :: center ! Center of circumcircle
    real(8) radius
    type(linked_list_item_type), pointer :: stub ! For fast Delaunay triangle deletion.
  contains
    procedure :: init => delaunay_triangle_init
    procedure :: add_DVT => delaunay_triangle_add_DVT
    procedure :: add_adjDT => delaunay_triangle_add_adjDT
    procedure :: print => delaunay_triangle_print
  end type delaunay_triangle_type

  type delaunay_edge_type
    integer :: id = -1
    type(delaunay_vertex_type), pointer :: DVT1
    type(delaunay_vertex_type), pointer :: DVT2
    type(delaunay_triangle_type), pointer :: DT1
    type(delaunay_triangle_type), pointer :: DT2
  contains
    procedure :: init => delaunay_edge_init
  end type delaunay_edge_type

  type, extends(point_type) :: voronoi_vertex_type
    integer :: id = -1
    type(delaunay_triangle_type), pointer :: DT
    type(array_type) VE
    type(array_type) VC
  contains
    procedure :: init => voronoi_vertex_init
    procedure :: print => voronoi_vertex_print
  end type voronoi_vertex_type

  type voronoi_edge_type
    integer :: id = -1
    type(voronoi_vertex_type), pointer :: VVT1
    type(voronoi_vertex_type), pointer :: VVT2
    type(voronoi_cell_type), pointer :: VC1
    type(voronoi_cell_type), pointer :: VC2
  contains
    procedure :: init => voronoi_edge_init
    procedure :: print => voronoi_edge_print
  end type voronoi_edge_type

  type voronoi_cell_type
    integer :: id = -1
    type(delaunay_vertex_type), pointer :: center
    type(array_type) VVT
    type(array_type) VE
    type(array_type) adjVC
  contains
    procedure :: init => voronoi_cell_init
    procedure :: print => voronoi_cell_print
  end type voronoi_cell_type

  ! Main data structure
  type(array_type) global_DVT_array
  type(linked_list_type) virtual_DVT_list ! Virtual Delaunay vertex array
  type(linked_list_type) global_DT_list
  type(array_type) global_DE_array
  type(array_type) global_VVT_array
  type(array_type) global_VE_array
  type(array_type) global_VC_array
  type(linked_list_type) tmpDT_list ! Temporal Delaunay triangles
  type(linked_list_type) obsDT_list ! Obsolete Delaunay triangles

  integer, parameter :: ORIENT_LEFT      = 1
  integer, parameter :: ORIENT_RIGHT     = 2
  integer, parameter :: ORIENT_ON        = 3
  integer, parameter :: INSIDE_TRIANGLE  = 4
  integer, parameter :: OUTSIDE_TRIANGLE = -4
  integer, parameter :: INSIDE_CIRCLE    = 1
  integer, parameter :: OUTSIDE_CIRCLE   = 2
  integer, parameter :: ON_CIRCLE        = 3

  integer, parameter :: im1(3) = [3, 1, 2], ip1(3) = [2, 3, 1]

  interface get_DVT
    module procedure get_DVT_from_list
    module procedure get_DVT_from_array
    module procedure get_DVT_from_value
  end interface get_DVT

  interface get_DT
    module procedure get_DT_from_list
    module procedure get_DT_from_array
    module procedure get_DT_from_value
  end interface get_DT

  interface get_DE
    module procedure get_DE_from_array
  end interface get_DE

  interface get_VVT
    module procedure get_VVT_from_array
    module procedure get_VVT_from_value
  end interface get_VVT

  interface get_VC
    module procedure get_VC_from_array
  end interface get_VC

  interface get_VE
    module procedure get_VE_from_array
    module procedure get_VE_from_value
  end interface get_VE

  interface orient
    module procedure orient1
    module procedure orient2
  end interface orient

  interface in_circle
    module procedure in_circle1
    module procedure in_circle2
    module procedure in_circle3
  end interface in_circle

  interface
    subroutine output_interface(DVT_array, DT_list, DE_array, VVT_array, VC_array, VE_array, tag)
      import array_type
      import linked_list_type
      type(array_type), intent(in) :: DVT_array
      type(linked_list_type), intent(in) :: DT_list
      type(array_type), intent(in) :: DE_array
      type(array_type), intent(in) :: VVT_array
      type(array_type), intent(in) :: VC_array
      type(array_type), intent(in) :: VE_array
      character(*), intent(in), optional :: tag
    end subroutine output_interface
  end interface

  ! For tagging DT.
  integer :: id_DT = 0

contains

  subroutine delaunay_voronoi_init(num_point, lon, lat, x, y, z)

    integer, intent(in) :: num_point
    real(8), intent(in), optional :: lon(:)
    real(8), intent(in), optional :: lat(:)
    real(8), intent(in), optional :: x  (:)
    real(8), intent(in), optional :: y  (:)
    real(8), intent(in), optional :: z  (:)

    type(delaunay_vertex_type) DVT
    type(voronoi_cell_type) VC
    integer i

    global_DVT_array = array(num_point)
    global_VC_array  = array(num_point)

    if (present(lon) .and. present(lat)) then
      do i = 1, num_point
        call DVT%init(id=i, lon=lon(i), lat=lat(i))
        call global_DVT_array%append(DVT)
        call VC%init(id=i)
        VC%center => get_DVT(global_DVT_array%last_value())
        call global_VC_array%append(VC)
      end do
    else if (present(x) .and. present(y) .and. present(z)) then
      do i = 1, num_point
        call DVT%init(id=i, x=x(i), y=y(i), z=z(i))
        call global_DVT_array%append(DVT)
        call VC%init(id=i)
        VC%center => get_DVT(global_DVT_array%last_value())
        call global_VC_array%append(VC)
      end do
    end if

    call random_number_init()

  end subroutine delaunay_voronoi_init

  subroutine delaunay_triangulation()

    integer i, j, ret, idx(3)
    type(delaunay_vertex_type), pointer :: DVT1, DVT2, DVT3
    type(delaunay_triangle_type), pointer :: DT1, DT2
    type(linked_list_iterator_type) DVT_iterator, DT_iterator
    logical inserted(global_DVT_array%size), found

    inserted = .false.

    call get_three_DVT_indices(idx)
    DVT1 => get_DVT(global_DVT_array, idx(1))
    DVT2 => get_DVT(global_DVT_array, idx(2))
    DVT3 => get_DVT(global_DVT_array, idx(3))

    ! Initialize the first eight triangles.
    call log_notice('Initialize the first eight triangles with selected vertices ' // trim(to_string(idx)) // '.')
    call init_DTs(DVT1, DVT2, DVT3)
    inserted(idx) = .true.

    ! Initialize the point-in-triangle relation between the rest vertices.
    DVT_loop: do i = 1, global_DVT_array%size
      if (inserted(i)) cycle
      DVT1 => get_DVT(global_DVT_array, i)
      DT_iterator = linked_list_iterator(global_DT_list)
      DT_loop: do while (.not. DT_iterator%ended())
        DT1 => get_DT(DT_iterator%value)
        ret = in_triangle_relaxed(DT1, DVT1)
        if (ret == INSIDE_TRIANGLE) then
          ! DVT is included in some DT .................................. CASE 1
          call record_included_DVT(DVT1, DT1)
          exit DT_loop
        else if (ret == OUTSIDE_TRIANGLE) then
        else if (ret > 0) then
          ! DVT is on some edge of DT ................................... CASE 2
          DT2 => get_DT(DT1%adjDT, ret)
          DVT1%edge_idx = ret
          call record_included_DVT(DVT1, DT1, DT2)
          exit DT_loop
        else if (ret < 0) then
          ! DVT coincides with some vertex of DT ........................ CASE 3
          DVT2 => get_DVT(DT1%DVT, -ret)
          if (DVT2%id < 0) then
            ! The vertex is virtual ................................... CASE 3-1
            ! It is ok, just replace it with DVT1.
            call extract_incident_DT_and_linked_DVT(DVT2)
            ! Replace the DVT from DVT2 to DVT1 for each incident DT.
            do j = 1, DVT2%incDT%size
              DT2 => get_DT(DVT2%incDT, j)
              call DT2%DVT%replace_ptr(DVT2, DVT1)
            end do
            ! Copy the first incident DT and link DVT from DVT2 to DVT1.
            call add_incident_DT(DVT1, get_DT(DVT2%incDT, 1))
            inserted(i) = .true.
            ! Delete DVT2 from virtual DVT 
            call virtual_DVT_list%delete_ptr(DVT2)
            exit DT_loop
          else
            ! The vertex is real ...................................... CASE 3-2
            ! Complain this degenerate case.
            call log_error('DVT ' // trim(to_string(DVT1%id)) // ' coincides to ' // &
                           'DVT ' // trim(to_string(-ret)) // ' of ' // &
                           'DT ' // trim(to_string(DT1%id)) // '.', __FILE__, __LINE__)
          end if
        end if
        call DT_iterator%next()
      end do DT_loop
    end do DVT_loop

    ! Insert the rest vertices one at a time.
    do i = 1, global_DVT_array%size
      if (inserted(i)) cycle
      inserted(i) = .true.
      DVT1 => get_DVT(global_DVT_array, i)
      ! Update the Delauny triangulation.
      call insert_DVT(DVT1)
      ! Update the point-in-triangle relation "locally"
      DT_iterator = linked_list_iterator(obsDT_list)
      do while (.not. DT_iterator%ended())
        DT1 => get_DT(DT_iterator%value)
        DVT_iterator = linked_list_iterator(DT1%incDVT)
        do while (.not. DVT_iterator%ended())
          DVT2 => get_DVT(DVT_iterator%value)
          call update_point_in_triangle(DVT2, DT1, found)
          if (.not. found) then
            call log_error('Failed to update DVT (' // trim(to_string(DVT2%id)) // &
              ') after delete DT (' // trim(to_string(DT1%id)) // '!')
          end if
          call DVT_iterator%next()
        end do
        call DT_iterator%next()
      end do
      ! Delete obsolete and temporal triangles
      call delete_obsolete_DT()
      call delete_temporal_DT()
    end do

    ! Extract the full list of incident DTs and link DVTs.
    do i = 1, global_DVT_array%size
      DVT1 => get_DVT(global_DVT_array, i)
      call extract_incident_DT_and_linked_DVT(DVT1)
    end do
    do i = 1, virtual_DVT_list%size
      DVT1 => get_DVT(virtual_DVT_list, i)
      call extract_incident_DT_and_linked_DVT(DVT1)
    end do

    ! Delete virtual DVTs.
    do i = 1, virtual_DVT_list%size
      DVT1 => get_DVT(virtual_DVT_list, i)
      call delete_DVT(DVT1)
    end do
    call virtual_DVT_list%clear()

    ! Reset id of DVT and DT.
    do i = 1, global_DVT_array%size
      DVT1 => get_DVT(global_DVT_array, i)
      DVT1%id = i
    end do
    i = 1
    DT_iterator = linked_list_iterator(global_DT_list)
    do while (.not. DT_iterator%ended())
      DT1 => get_DT(DT_iterator%value)
      DT1%id = i
      i = i + 1
      call DT_iterator%next()
    end do

    call log_notice('Delaunay triangulation is done.')

  end subroutine delaunay_triangulation

  subroutine voronoi_diagram(all)

    logical, intent(in), optional :: all ! If extract all topological information or not.

    type(linked_list_iterator_type) DT_iterator, DVT_iterator
    type(delaunay_triangle_type), pointer :: DT, DT1, DT2, adjDT
    type(delaunay_vertex_type), pointer :: DVT, linkDVT
    type(delaunay_edge_type) DE
    type(voronoi_vertex_type) VVT
    type(voronoi_vertex_type), pointer :: VVT1, VVT2
    type(voronoi_edge_type) VE
    type(voronoi_edge_type), pointer :: VE1, VE2
    type(voronoi_cell_type), pointer :: VC, VC1, VC2, adjVC
    integer i, j, k
    integer j_first_VE, j_VE
    integer j_first_VC, j_VC
    integer VE_idx(6)
    integer VC_idx(3)
    logical matched

    ! Calculate VVT which is the circumcenter of DT.
    global_VVT_array = array(global_DT_list%size)
    DT_iterator = linked_list_iterator(global_DT_list)
    i = 1
    do while (.not. DT_iterator%ended())
      DT => get_DT(DT_iterator%value)
      call VVT%init(id=i)
      call global_VVT_array%append(VVT)
      VVT1 => get_VVT(global_VVT_array%last_value())
      VVT1%DT => DT
      call calc_circumcircle(DT, VVT1)
      call DT_iterator%next()
      i = i + 1
    end do

    ! Connect VC and VVT.
    do i = 1, global_VC_array%size
      VC => get_VC(global_VC_array, i)
      DVT => VC%center
      ! Loop on incDT of DVT.
      DT_iterator = linked_list_iterator(DVT%incDT)
      DVT_iterator = linked_list_iterator(DVT%linkDVT)
      do while (.not. DT_iterator%ended(cyclic=.false.))
        linkDVT => get_DVT(DVT_iterator%value)
        adjVC => get_VC(global_VC_array, linkDVT%id)
        call VC%adjVC%append_ptr(adjVC)
        DT => get_DT(DT_iterator%value)
        call VC%VVT%append_ptr(get_VVT(global_VVT_array, DT%id))
        call DT_iterator%next()
        call DVT_iterator%next()
      end do
    end do

    if (present(all) .and. all) then
      ! Create Delaunay edges.
      global_DE_array = array(euler_formula(num_cell=global_DT_list%size, num_vertex=global_DVT_array%size))
      k = 1
      do i = 1, global_DVT_array%size
        DVT => get_DVT(global_DVT_array, i)
        DVT_iterator = linked_list_iterator(DVT%linkDVT)
        DT_iterator = linked_list_iterator(DVT%incDT)
        do while (.not. DVT_iterator%ended(cyclic=.false.))
          linkDVT => get_DVT(DVT_iterator%value)
          if (linkDVT%id > DVT%id) then
            call DE%init(id=k)
            DE%DVT1 => DVT
            DE%DVT2 => linkDVT
            DE%DT1 => get_DT(DT_iterator%value)
            DE%DT2 => get_DT(DT_iterator%item%prev%value)
            call global_DE_array%append(DE)
            k = k + 1
          end if
          call DVT_iterator%next()
          call DT_iterator%next()
        end do
      end do
      ! Create Voronoi edges.
      global_VE_array = array(euler_formula(num_cell=global_VC_array%size, num_vertex=global_VVT_array%size))
      k = 1
      do i = 1, global_VC_array%size
        VC => get_VC(global_VC_array, i)
        do j = 1, VC%adjVC%size
          adjVC => get_VC(VC%adjVC, j)
          ! Avoid add VE repeatedly.
          if (adjVC%id > VC%id) then
            call VE%init(id=k)
            VE%VVT1 => get_VVT(VC%VVT, j)
            VE%VVT2 => get_VVT(VC%VVT, merge(j + 1, 1, j /= VC%adjVC%size))
            VE%VC1 => VC
            VE%VC2 => adjVC
            ! Link VVT and VC.
            call VE%VVT1%VC%append_ptr(VC)
            call VE%VVT2%VC%append_ptr(adjVC)
            call global_VE_array%append(VE)
            ! We need to append the pointer to the inserted VE.
            call VE%VC1%VE%append_ptr(global_VE_array%last_value())
            call VE%VC2%VE%append_ptr(global_VE_array%last_value())
            call VE%VVT1%VE%append_ptr(global_VE_array%last_value())
            call VE%VVT2%VE%append_ptr(global_VE_array%last_value())
            k = k + 1
          end if
        end do
      end do
      ! Fix the order of VE in VC.
      do i = 1, global_VC_array%size
        VC => get_VC(global_VC_array, i)
        ! Find out the edge whose vertices contains the first vertex of VC.
        VVT1 => get_VVT(VC%VVT, 1)
        j_first_VE = -1
        do j = 1, VC%VE%size
          VE1 => get_VE(VC%VE, j)
          ! Assume VC is on the left side of VE.
          if ((VE1%VC1%id == VC%id .and. VE1%VVT1%id == VVT1%id) .or. (VE1%VC2%id == VC%id .and. VE1%VVT2%id == VVT1%id)) then
            j_first_VE = j
            VE_idx(1) = VE1%id
            exit
          end if
        end do
        if (j_first_VE == -1) call log_error('Internal error!', __FILE__, __LINE__)
        j_VE = 2
        do j = 1, VC%VE%size
          if (j == j_first_VE) cycle
          matched = .false.
          do k = 1, VC%VE%size
            VE2 => get_VE(VC%VE, k)
            if (VE1%id == VE2%id) cycle
            ! Keep VC on the left side of VE.
            if (VE1%VC1%id == VC%id) then
              matched = (VE2%VC1%id == VC%id .and. VE1%VVT2%id == VE2%VVT1%id) .or. (VE2%VC2%id == VC%id .and. VE1%VVT2%id == VE2%VVT2%id)
            else if (VE1%VC2%id == VC%id) then
              matched = (VE2%VC1%id == VC%id .and. VE1%VVT1%id == VE2%VVT1%id) .or. (VE2%VC2%id == VC%id .and. VE1%VVT1%id == VE2%VVT2%id)
            end if
            if (matched) then
              VE_idx(j_VE) = VE2%id
              j_VE = j_VE + 1
              VE1 => VE2
              exit
            end if
          end do
        end do
        ! Replace the VE pointers.
        do j = 1, VC%VE%size
          call VC%VE%replace_ptr_at(j, get_VE(global_VE_array, VE_idx(j)))
        end do
      end do
      ! Fix the order of VE in VVT.
      do i = 1, global_VVT_array%size
        VVT1 => get_VVT(global_VVT_array, i)
        VE1 => get_VE(VVT1%VE, 1)
        VE_idx(1) = VE1%id
        do j = 2, VVT1%VE%size
          matched = .false.
          do k = 2, VVT1%VE%size
            VE2 => get_VE(VVT1%VE, k)
            if (VE1%id == VE2%id) cycle
            ! Keep VC on the left side of VE.
            if (VE1%VVT1%id == VVT1%id) then
              matched = (VE2%VVT1%id == VVT1%id .and. VE1%VC1%id == VE2%VC2%id) .or. (VE2%VVT2%id == VVT1%id .and. VE1%VC1%id == VE2%VC1%id)
            else if (VE1%VVT2%id == VVT1%id) then
              matched = (VE2%VVT1%id == VVT1%id .and. VE1%VC2%id == VE2%VC2%id) .or. (VE2%VVT2%id == VVT1%id .and. VE1%VC2%id == VE2%VC1%id)
            end if
            if (matched) then
              VE_idx(j) = VE2%id
              VE1 => VE2
              exit
            end if
          end do
        end do
        ! Replace the VE pointers.
        do j = 1, VVT1%VE%size
          call VVT1%VE%replace_ptr_at(j, get_VE(global_VE_array, VE_idx(j)))
        end do
      end do
      ! Fix the order of VC in VVT.
      do i = 1, global_VVT_array%size
        VVT1 => get_VVT(global_VVT_array, i)
        ! Find out the edge whose vertices contains the first vertex of VC.
        VE1 => get_VE(VVT1%VE, 1)
        j_first_VC = -1
        VC_idx = -1
        do j = 1, VVT1%VC%size
          VC1 => get_VC(VVT1%VC, j)
          if ((VE1%VVT1%id == VVT1%id .and. VE1%VC1%id == VC1%id) .or. (VE1%VVT2%id == VVT1%id .and. VE1%VC2%id == VC1%id)) then
            j_first_VC = j
            VC_idx(1) = VC1%id
            exit
          end if
        end do
        if (j_first_VC == -1) call log_error('Internal error!', __FILE__, __LINE__)
        j_VC = 2
        VE1 => get_VE(VVT1%VE, j_VC)
        do j = 1, VVT1%VC%size
          if (j == j_first_VC) cycle
          do k = 1, VVT1%VC%size
            VC2 => get_VC(VVT1%VC, k)
            if (VC1%id == VC2%id) cycle
            ! Keep VC on the left side of VE.
            if ((VE1%VVT1%id == VVT1%id .and. VE1%VC1%id == VC2%id) .or. (VE1%VVT2%id == VVT1%id .and. VE1%VC2%id == VC2%id)) then
              VC_idx(j_VC) = VC2%id
              j_VC = j_VC + 1
              if (j_VC <= VVT1%VC%size) then
                VC1 => VC2
                VE1 => get_VE(VVT1%VE, j_VC)
              end if
              exit
            end if
          end do
        end do
        ! Replace the VE pointers.
        do j = 1, VVT1%VC%size
          call VVT1%VC%replace_ptr_at(j, get_VC(global_VC_array, VC_idx(j)))
        end do
      end do
    end if

  end subroutine voronoi_diagram

  subroutine delaunay_voronoi_output(output)

    procedure(output_interface) output

    call output(global_DVT_array, global_DT_list, global_DE_array, global_VVT_array, global_VC_array, global_VE_array)

  end subroutine delaunay_voronoi_output

  subroutine get_three_DVT_indices(idx)

    integer, intent(out) :: idx(3)

    logical :: success = .false.

    do while (.not. success)
      call random_number_get(1, global_DVT_array%size, idx)
      if (idx(1) /= idx(2) .and. &
          idx(1) /= idx(3) .and. &
          idx(2) /= idx(3)) then
        success = .true.
      end if
    end do

  end subroutine get_three_DVT_indices

  subroutine init_DTs(DVT1, DVT2, DVT3)

    type(delaunay_vertex_type), intent(inout) :: DVT1
    type(delaunay_vertex_type), intent(inout) :: DVT2
    type(delaunay_vertex_type), intent(inout) :: DVT3

    type(delaunay_vertex_type) local_DVT
    type(delaunay_vertex_type), pointer :: DVT4, DVT5, DVT6
    type(delaunay_triangle_type), pointer :: DT
    type(array_type) DVT_array
    type(array_type) DT_array
    real(8) lon, lat
    integer i, j, map(24)

    ! Create virtual vertices which are 'antipodal' to the corresponding
    ! vertices of the first three inserted ones.
    call virtual_DVT_list%append(local_DVT)
    DVT4 => get_DVT(virtual_DVT_list%last_value())
    call inverse_rotation_transform(DVT1%lon, DVT1%lat, lon, lat, 0.0d0, -pi05)
    call DVT4%init(id=-1, lon=lon, lat=lat)
    call virtual_DVT_list%append(local_DVT)
    DVT5 => get_DVT(virtual_DVT_list%last_value())
    call inverse_rotation_transform(DVT2%lon, DVT2%lat, lon, lat, 0.0d0, -pi05)
    call DVT5%init(id=-2, lon=lon, lat=lat)
    call virtual_DVT_list%append(local_DVT)
    DVT6 => get_DVT(virtual_DVT_list%last_value())
    call inverse_rotation_transform(DVT3%lon, DVT3%lat, lon, lat, 0.0d0, -pi05)
    call DVT6%init(id=-3, lon=lon, lat=lat)

    call DVT_array%append_ptr(DVT1)
    call DVT_array%append_ptr(DVT2)
    call DVT_array%append_ptr(DVT3)
    call DVT_array%append_ptr(DVT4)
    call DVT_array%append_ptr(DVT5)
    call DVT_array%append_ptr(DVT6)

    ! Create the first eight triangles.
    do i = 1, 8
      call create_DT(DT)
      call DT%init(id=i)
      call DT_array%append_ptr(DT)
    end do

    ! Set neightbors for each triangle.
    map = [5,2,4, 6,3,1, 7,4,2, 8,1,3, 1,8,6, 2,5,7, 3,6,8, 4,7,5]
    j = 1
    do i = 1, 8
      DT => get_DT(DT_array, i)
      call DT%add_adjDT(get_DT(DT_array, map(j)), get_DT(DT_array, map(j+1)), get_DT(DT_array, map(j+2)))
      j = j + 3
    end do

    ! Set vertices for each triangle.
    select case (orient(DVT1, DVT2, DVT3))
    case (ORIENT_LEFT)
      map = [1,2,3, 1,3,5, 1,5,6, 1,6,2, 4,3,2, 4,5,3, 4,6,5, 4,2,6]
    case (ORIENT_RIGHT)
      map = [1,3,2, 1,2,6, 1,6,5, 1,5,3, 4,2,3, 4,6,2, 4,5,6, 4,3,5]
    case (ORIENT_ON)
      call log_error('The first three points are collinear.' // &
                     'Change the order of points to ensure the non-collinearity ' // &
                     'of the first three ones. If you are sure they are not, ' // &
                     'this may be caused by the failure of floating-point' // &
                     'calculation.')
    end select
    j = 1
    do i = 1, 8
      DT => get_DT(DT_array, i)
      call DT%add_DVT(get_DVT(DVT_array, map(j)), get_DVT(DVT_array, map(j+1)), get_DVT(DVT_array, map(j+2)))
      j = j + 3
    end do

    ! Set one incident triangle for each vertex.
    ! Note: Any one is ok, the full list will be extracted at the end of
    !       the Delaunay triangulation.
    call add_incident_DT(DVT1, get_DT(DT_array, 1))
    call add_incident_DT(DVT2, get_DT(DT_array, 1))
    call add_incident_DT(DVT3, get_DT(DT_array, 1))
    call add_incident_DT(DVT4, get_DT(DT_array, 7))
    call add_incident_DT(DVT5, get_DT(DT_array, 7))
    call add_incident_DT(DVT6, get_DT(DT_array, 7))

  end subroutine init_DTs

  ! ************************************************************************** !
  ! extract_incident_DT_and_linked_DVT                                         !
  ! Purpose:                                                                   !
  !   During the Delaunay triangulation, only one incident triangle of each    !
  !   vertex is recorded. This subroutine extract all the incident triangles   !
  !   and link vertices to form a ring in counter-clockwise order.             !
  ! ************************************************************************** !

  subroutine extract_incident_DT_and_linked_DVT(DVT)

    type(delaunay_vertex_type), intent(inout), target :: DVT

    type(delaunay_triangle_type), pointer :: incDT0, incDT1
    type(delaunay_vertex_type), pointer :: linkDVT0, linkDVT1

    integer i

    ! Get the first incident DT which should be set already.
    incDT1 => get_DT(DVT%incDT, 1)
    incDT0 => incDT1
    linkDVT0 => null()
    incident_DT_loop: do
      ! Loop each DVT of incident DT.
      DVT_loop: do i = 1, 3
        linkDVT1 => get_DVT(incDT1%DVT, i)
        if (associated(linkDVT1, DVT)) then
          ! linkDVT is the target DVT, so shift to the previous DVT.
          linkDVT1 => get_DVT(incDT1%DVT, im1(i))
          if (.not. associated(linkDVT0)) then
            ! Record the first one.
            linkDVT0 => linkDVT1
          else if (associated(linkDVT1, linkDVT0)) then
            ! The ring has formed.
            call DVT%linkDVT%cycle()
            call DVT%incDT%cycle()
            exit incident_DT_loop
          end if
          call DVT%linkDVT%append_ptr(linkDVT1)
          exit DVT_loop
        end if
      end do DVT_loop
      ! Shift to next incident DT.
      incDT1 => get_DT(incDT1%adjDT, ip1(i))
      ! Avoid add duplicated incident DT.
      if (.not. associated(incDT1, incDT0)) call add_incident_DT(DVT, incDT1)
    end do incident_DT_loop

  end subroutine extract_incident_DT_and_linked_DVT

  ! ************************************************************************** !
  ! insert_DVT                                                                 !
  ! Purpose:                                                                   !
  !   Insert one vertex into the Delaunay triangulation. It is the core of     !
  !   the incremental algorithm.                                               !
  ! ************************************************************************** !

  subroutine insert_DVT(DVT)

    type(delaunay_vertex_type), intent(inout), target :: DVT

    if (associated(DVT%cntDT2)) then
      ! DVT is on the edge shared by two adjacent triangles.
      call flip24(DVT%cntDT1, DVT%cntDT2, DVT)
    else
      ! DVT is only contained by one triangle.
      call flip13(DVT%cntDT1, DVT)
    end if

  end subroutine insert_DVT

  ! ************************************************************************** !
  ! delete_DVT                                                                 !
  ! Purpose:                                                                   !
  !   The inverse of insert_DVT. It is used for deleting virtual vertices.     !
  ! ************************************************************************** !

  subroutine delete_DVT(DVT)

    type(delaunay_vertex_type), intent(inout), target :: DVT

    type(linked_list_item_type), pointer :: incDT_item
    type(linked_list_item_type), pointer :: linkDVT_item1, linkDVT_item2
    type(delaunay_triangle_type), pointer :: DT1, DT2, DT3, newDT1, newDT2
    type(delaunay_vertex_type), pointer :: DVT1, DVT2, DVT3, DVT4
    integer i, j, k
    logical empty

    ! Delete incident DTs to three ones.
    incDT_item => DVT%incDT%first_item
    linkDVT_item1 => DVT%linkDVT%first_item
    do while (DVT%incDT%size > 3)
      DT1 => get_DT(incDT_item%value)
      DT2 => get_DT(incDT_item%next%value)
      DVT1 => get_DVT(linkDVT_item1%value)
      DVT2 => get_DVT(linkDVT_item1%next%value)
      DVT3 => get_DVT(linkDVT_item1%next%next%value)
      ! Check 1: Is potential DT convex?
      if (orient(DVT1, DVT2, DVT3) == ORIENT_RIGHT) then
        incDT_item => incDT_item%next
        linkDVT_item1 => linkDVT_item1%next
        cycle ! NOT PASS
      end if
      ! Check 2: Does potential DT contain DVT?
      if (orient(DVT3, DVT1, DVT) == ORIENT_LEFT) then
        incDT_item => incDT_item%next
        linkDVT_item1 => linkDVT_item1%next
        cycle ! NOT PASS
      end if
      ! Check 3: Does potential DT satisfy empty-circumcircle rule?
      empty = .true.
      linkDVT_item2 => linkDVT_item1%next%next%next
      do i = 1, DVT%incDT%size - 3
        DVT4 => get_DVT(linkDVT_item2%value)
        select case (in_circle(DVT1, DVT2, DVT3, DVT4))
        case (INSIDE_CIRCLE)
          empty = .false.
          exit
        case (ON_CIRCLE)
          ! call log_warning('Encounter cocircular vertices.', __FILE__, __LINE__)
        end select
        linkDVT_item2 => linkDVT_item2%next
      end do
      if (.not. empty) then
        incDT_item => incDT_item%next
        linkDVT_item1 => linkDVT_item1%next
        cycle ! NOT PASS
      end if
      ! So far, the potential DT is a real DT.
      i = DT1%DVT%index_ptr(DVT)
      j = DT2%DVT%index_ptr(DVT)
      call flip22(DT1, DT2, newDT1, newDT2, [im1(i), i, ip1(i), im1(j)])
      call global_DT_list%remove_item(DT1%stub)
      call global_DT_list%remove_item(DT2%stub)
    end do

    ! Do the final flip31.
    DT1 => get_DT(DVT%incDT%first_item%value)
    DT2 => get_DT(DVT%incDT%first_item%next%value)
    DT3 => get_DT(DVT%incDT%first_item%next%next%value)
    i = DT1%DVT%index_ptr(DVT)
    j = DT2%DVT%index_ptr(DVT)
    k = DT3%DVT%index_ptr(DVT)
    call flip31(DT1, DT2, DT3, newDT1, [ip1(i), im1(i), i, j, k])
    call global_DT_list%remove_item(DT1%stub)
    call global_DT_list%remove_item(DT2%stub)
    call global_DT_list%remove_item(DT3%stub)

  end subroutine delete_DVT

  recursive subroutine validate_DT(newDT)

    type(delaunay_triangle_type), intent(inout), target :: newDT

    type(delaunay_triangle_type), pointer :: opsDT, newerDT1, newerDT2
    type(delaunay_vertex_type), pointer :: opsDVT
    integer i

    ! Assume 3rd adjDT is the existing DT.
    opsDT => get_DT(newDT%adjDT, 3)
    ! Get the opposite DVT to the newly inserted DVT with respect to
    ! the shared edge by newDT and opsDT.
    do i = 1, 3
      if (associated(get_DT(opsDT%adjDT, i), newDT)) then
        opsDVT => get_DVT(opsDT%DVT, i)
        exit
      end if
    end do
    ! Check if the opposite DVT is in the circumcircle of DT.
    select case (in_circle(newDT, opsDVT))
    case (OUTSIDE_CIRCLE)
      return ! THE ONLY EXIT
    case (INSIDE_CIRCLE)
      call flip22(newDT, opsDT, newerDT1, newerDT2, [1, 2, 3, i])
      ! Recursively validate the two newer DT.
      call validate_DT(newerDT1)
      call validate_DT(newerDT2)
      call record_temporal_DT(newDT)
      call newDT%subDT%append_ptr(newerDT1)
      call newDT%subDT%append_ptr(newerDT2)
      call record_obsolete_DT(opsDT)
      call opsDT%subDT%append_ptr(newerDT1)
      call opsDT%subDT%append_ptr(newerDT2)
    case (ON_CIRCLE)
      ! call log_warning('Encounter cocircular vertices!', __FILE__, __LINE__)
    end select

  end subroutine validate_DT

  recursive subroutine update_point_in_triangle(DVT, DT, found)

    type(delaunay_vertex_type), intent(inout) :: DVT
    type(delaunay_triangle_type), intent(inout), target :: DT
    logical, intent(inout) :: found

    type(delaunay_vertex_type), pointer :: DVT1, DVT2, DVT3
    type(delaunay_triangle_type), pointer :: DT1
    integer i, on_plane(2)

    if (DT%subDT%size > 0) then
      ! DT has been subdivided, go to next level ....................... ENTRY 1
      do i = 1, DT%subDT%size
        DT1 => get_DT(DT%subDT, i)
        call update_point_in_triangle(DVT, DT1, found)
        if (found) return
      end do
    else
      ! DT is at the last level without being subdivided ............... ENTRY 2
      on_plane = 0
      DVT1 => get_DVT(DT%DVT, 1)
      DVT2 => get_DVT(DT%DVT, 2)
      DVT3 => get_DVT(DT%DVT, 3)
      ! Check edge 3->1
      select case (orient(DVT3, DVT1, DVT))
      case (ORIENT_RIGHT)
        found = .false.
        return
      case (ORIENT_ON)
        on_plane(1) = 2
      end select
      ! Check edge 2->3
      select case (orient(DVT2, DVT3, DVT))
      case (ORIENT_RIGHT)
        found = .false.
        return
      case (ORIENT_ON)
        on_plane(2) = 1
      end select
      ! Summary
      found = .true.
      ! Note: Delete DVT from its contained DTs to avoid duplicate update.
      if (associated(DVT%cntDT2)) then
        call delete_included_DVT(DVT, DVT%cntDT1, DVT%cntDT2)
      else
        call delete_included_DVT(DVT, DVT%cntDT1)
      end if
      if (on_plane(1) == 0 .and. on_plane(2) == 0) then
        call record_included_DVT(DVT, DT)
      else if (on_plane(1) == 2 .and. on_plane(2) == 0) then
        DVT%edge_idx = on_plane(1)
        DT1 => get_DT(DT%adjDT, 2)
        call record_included_DVT(DVT, DT, DT1)
      else if (on_plane(1) == 0 .and. on_plane(2) == 1) then
        DVT%edge_idx = on_plane(2)
        DT1 => get_DT(DT%adjDT, 1)
        call record_included_DVT(DVT, DT, DT1)
      else
        call log_error('DVT ' // trim(to_string(DVT%id)) // &
          ' coincides with vertex 3 of DT ' // trim(to_string(DT%id)) // '!', &
          __FILE__, __LINE__)
      end if
    end if

  end subroutine update_point_in_triangle

  ! ************************************************************************** !
  ! Flip??                                                                     !
  ! Purpose:                                                                   !
  !   These flips are the basic operations in the action of insertion and      !
  !   deletion of vertices.                                                    !
  ! ************************************************************************** !

  subroutine flip13(oldDT, DVT)

    type(delaunay_triangle_type), intent(inout) :: oldDT
    type(delaunay_vertex_type), intent(inout) :: DVT

    type(delaunay_vertex_type), pointer :: DVT1, DVT2, DVT3
    type(delaunay_triangle_type), pointer :: adjDT1, adjDT2, adjDT3
    type(delaunay_triangle_type), pointer :: newDT1, newDT2, newDT3

    DVT1 => get_DVT(oldDT%DVT, 1)
    DVT2 => get_DVT(oldDT%DVT, 2)
    DVT3 => get_DVT(oldDT%DVT, 3)
    adjDT1 => get_DT(oldDT%adjDT, 3)
    adjDT2 => get_DT(oldDT%adjDT, 1)
    adjDT3 => get_DT(oldDT%adjDT, 2)
    ! Subdivide the old triangle into three new ones.
    call create_DT(newDT1)
    call create_DT(newDT2)
    call create_DT(newDT3)
    ! Set up the topology of the three triangles.
    call newDT1%add_DVT(DVT1, DVT2, DVT)
    call newDT1%add_adjDT(newDT2, newDT3, adjDT1)
    call newDT2%add_DVT(DVT2, DVT3, DVT)
    call newDT2%add_adjDT(newDT3, newDT1, adjDT2)
    call newDT3%add_DVT(DVT3, DVT1, DVT)
    call newDT3%add_adjDT(newDT1, newDT2, adjDT3)
    ! Link the newly inserted DVT to one of its incident DT
    call add_incident_DT(DVT, newDT1)
    ! Make change to the old triangle's adjDT.
    call adjDT1%adjDT%replace_ptr(oldDT, newDT1)
    call adjDT2%adjDT%replace_ptr(oldDT, newDT2)
    call adjDT3%adjDT%replace_ptr(oldDT, newDT3)
    call DVT1%incDT%replace_ptr(oldDT, newDT1, fatal=.false.)
    call DVT2%incDT%replace_ptr(oldDT, newDT2, fatal=.false.)
    call DVT3%incDT%replace_ptr(oldDT, newDT3, fatal=.false.)
    ! Validate the three triangles.
    call validate_DT(newDT1)
    call validate_DT(newDT2)
    call validate_DT(newDT3)
    ! Record the obsolete triangle.
    call record_obsolete_DT(oldDT)
    call oldDT%subDT%append_ptr(newDT1)
    call oldDT%subDT%append_ptr(newDT2)
    call oldDT%subDT%append_ptr(newDT3)
    ! Remove DVT from point-in-triangle relation.
    call delete_included_DVT(DVT, oldDT)

  end subroutine flip13

  subroutine flip24(oldDT1, oldDT2, DVT)

    type(delaunay_triangle_type), intent(inout), target :: oldDT1
    type(delaunay_triangle_type), intent(inout), target :: oldDT2
    type(delaunay_vertex_type), intent(inout) :: DVT

    type(delaunay_vertex_type), pointer :: DVT1, DVT2, DVT3, DVT4
    type(delaunay_triangle_type), pointer :: adjDT1, adjDT2, adjDT3, adjDT4
    type(delaunay_triangle_type), pointer :: newDT1, newDT2, newDT3, newDT4, DT
    integer i

    do i = 1, 3
      if (associated(get_DT(oldDT2%adjDT, i), oldDT1)) exit
    end do
    DVT1 => get_DVT(oldDT1%DVT, DVT%edge_idx)
    DVT2 => get_DVT(oldDT1%DVT, ip1(DVT%edge_idx))
    DVT3 => get_DVT(oldDT1%DVT, im1(DVT%edge_idx))
    DVT4 => get_DVT(oldDT2%DVT, i)
    adjDT1 => get_DT(oldDT1%adjDT, ip1(DVT%edge_idx))
    adjDT2 => get_DT(oldDT1%adjDT, im1(DVT%edge_idx))
    adjDT3 => get_DT(oldDT2%adjDT, ip1(i))
    adjDT4 => get_DT(oldDT2%adjDT, im1(i))
    ! Subdivide the two triangles into four new ones
    call create_DT(newDT1)
    call create_DT(newDT2)
    call create_DT(newDT3)
    call create_DT(newDT4)
    ! Set up the topology of the four triangles
    ! For oldDT1's newDTs
    call newDT1%add_DVT(DVT3, DVT1, DVT)
    call newDT1%add_adjDT(newDT2, newDT4, adjDT1)
    call newDT2%add_DVT(DVT1, DVT2, DVT)
    call newDT2%add_adjDT(newDT3, newDT1, adjDT2)
    ! For oldDT2's newDTs
    call newDT3%add_DVT(DVT2, DVT4, DVT)
    call newDT3%add_adjDT(newDT4, newDT2, adjDT3)
    call newDT4%add_DVT(DVT4, DVT3, DVT)
    call newDT4%add_adjDT(newDT1, newDT3, adjDT4)
    ! Link newly inserted DVT to its incident DTs.
    call add_incident_DT(DVT, newDT1)
    ! Make change to the oldDTs' adjDT and DVT.
    ! - adjDT
    do i = 1, 3
      DT => get_DT(adjDT1%adjDT, i)
      if (associated(DT, oldDT1)) then
        call adjDT1%adjDT%insert_ptr_at(i, newDT1)
        exit
      end if
    end do
    do i = 1, 3
      DT => get_DT(adjDT2%adjDT, i)
      if (associated(DT, oldDT1)) then
        call adjDT2%adjDT%insert_ptr_at(i, newDT2)
        exit
      end if
    end do
    do i = 1, 3
      DT => get_DT(adjDT3%adjDT, i)
      if (associated(DT, oldDT2)) then
        call adjDT3%adjDT%insert_ptr_at(i, newDT3)
        exit
      end if
    end do
    do i = 1, 3
      DT => get_DT(adjDT4%adjDT, i)
      if (associated(DT, oldDT2)) then
        call adjDT4%adjDT%insert_ptr_at(i, newDT4)
        exit
      end if
    end do
    ! - DVT
    call DVT1%incDT%replace_ptr(oldDT1, newDT1, fatal=.false.)
    call DVT2%incDT%replace_ptr(oldDT1, newDT2, old_value2=oldDT2, fatal=.false.)
    call DVT3%incDT%replace_ptr(oldDT1, newDT4, old_value2=oldDT2, fatal=.false.)
    call DVT4%incDT%replace_ptr(oldDT2, newDT3, fatal=.false.)
    ! Validate the new DTs.
    call validate_DT(newDT1)
    call validate_DT(newDT2)
    call validate_DT(newDT3)
    call validate_DT(newDT4)
    ! Record the obsolete DTs.
    call record_obsolete_DT(oldDT1)
    call oldDT1%subDT%append_ptr(newDT1)
    call oldDT1%subDT%append_ptr(newDT2)
    call record_obsolete_DT(oldDT2)
    call oldDT2%subDT%append_ptr(newDT3)
    call oldDT2%subDT%append_ptr(newDT4)
    ! Remove DVT from point-in-triangle relation.
    call delete_included_DVT(DVT, oldDT1, oldDT2)

  end subroutine flip24

  subroutine flip22(oldDT1, oldDT2, newDT1, newDT2, idx_map)

    type(delaunay_triangle_type), intent(in), target :: oldDT1
    type(delaunay_triangle_type), intent(in), target :: oldDT2
    type(delaunay_triangle_type), intent(out), pointer :: newDT1
    type(delaunay_triangle_type), intent(out), pointer :: newDT2
    integer, intent(in) :: idx_map(4)

    type(delaunay_vertex_type), pointer :: DVT1, DVT2, DVT3, DVT4
    type(delaunay_triangle_type), pointer :: adjDT1, adjDT2, adjDT3, adjDT4

    DVT1 => get_DVT(oldDT1%DVT, idx_map(1))
    DVT2 => get_DVT(oldDT1%DVT, idx_map(2))
    DVT3 => get_DVT(oldDT1%DVT, idx_map(3))
    DVT4 => get_DVT(oldDT2%DVT, idx_map(4))
    adjDT1 => get_DT(oldDT1%adjDT, idx_map(1))
    adjDT2 => get_DT(oldDT1%adjDT, idx_map(2))
    adjDT3 => get_DT(oldDT2%adjDT, ip1(idx_map(4)))
    adjDT4 => get_DT(oldDT2%adjDT, im1(idx_map(4)))
    ! Create two new DTs and set up their topology
    call create_DT(newDT1)
    call create_DT(newDT2)
    call newDT1%add_DVT(DVT1, DVT4, DVT3)
    call newDT1%add_adjDT(newDT2, adjDT2, adjDT3)
    call newDT2%add_DVT(DVT4, DVT2, DVT3)
    call newDT2%add_adjDT(adjDT1, newDT1, adjDT4)
    ! Make change to the old triangle's adjDTs and DVTs
    call adjDT1%adjDT%replace_ptr(oldDT1, newDT2)
    call adjDT2%adjDT%replace_ptr(oldDT1, newDT1)
    call adjDT3%adjDT%replace_ptr(oldDT2, newDT1)
    call adjDT4%adjDT%replace_ptr(oldDT2, newDT2)
    if (DVT1%incDT%cyclic()) then
      ! The full list of incident DTs has been extracted.
      ! This may be called from delete_DVT.
      call merge_incident_DT(DVT1, oldDT1, oldDT2, newDT1)
      call merge_incident_DT(DVT2, oldDT1, oldDT2, newDT2)
      call split_incident_DT(DVT3, oldDT1, newDT1, newDT2)
      call split_incident_DT(DVT4, oldDT2, newDT2, newDT1)
      call delete_link_DVT(DVT1, DVT2)
      call delete_link_DVT(DVT2, DVT1)
      call add_link_DVT(DVT3, DVT4, DVT1, DVT2)
      call add_link_DVT(DVT4, DVT3, DVT2, DVT1)
    else
      ! The full list of incident DTs has not been extracted.
      ! This may be called from validate_DT.
      call DVT1%incDT%insert_ptr_at(1, newDT1)
      call DVT2%incDT%insert_ptr_at(1, newDT2)
      call DVT3%incDT%insert_ptr_at(1, newDT2)
      call DVT4%incDT%insert_ptr_at(1, newDT1)
    end if

  end subroutine flip22

  subroutine flip31(oldDT1, oldDT2, oldDT3, newDT, idx_map)

    type(delaunay_triangle_type), intent(in), target :: oldDT1, oldDT2, oldDT3
    type(delaunay_triangle_type), intent(out), pointer :: newDT
    integer, intent(in) :: idx_map(5)

    type(delaunay_vertex_type), pointer :: DVT1, DVT2, DVT3, DVT4
    type(delaunay_triangle_type), pointer :: adjDT1, adjDT2, adjDT3

    DVT1 => get_DVT(oldDT1%DVT, idx_map(1))
    DVT2 => get_DVT(oldDT1%DVT, idx_map(2))
    DVT3 => get_DVT(oldDT1%DVT, idx_map(3))
    DVT4 => get_DVT(oldDT2%DVT, im1(idx_map(4)))
    adjDT1 => get_DT(oldDT1%adjDT, idx_map(3))
    adjDT2 => get_DT(oldDT2%adjDT, idx_map(4))
    adjDT3 => get_DT(oldDT3%adjDT, idx_map(5))
    ! Create a new DT and set up its topology.
    call create_DT(newDT)
    call newDT%add_DVT(DVT1, DVT2, DVT4)
    call newDT%add_adjDT(adjDT2, adjDT3, adjDT1)
    ! Make change to the old triangle's adjDTs and DVTs.
    call adjDT1%adjDT%replace_ptr(oldDT1, newDT)
    call adjDT2%adjDT%replace_ptr(oldDT2, newDT)
    call adjDT3%adjDT%replace_ptr(oldDT3, newDT)
    ! This should be called after the full list of incident DTs
    ! has been extracted. (Called in delete_DVT).
    call merge_incident_DT(DVT1, oldDT1, oldDT3, newDT)
    call merge_incident_DT(DVT2, oldDT2, oldDT1, newDT)
    call merge_incident_DT(DVT4, oldDT3, oldDT2, newDT)
    call delete_link_DVT(DVT1, DVT3)
    call delete_link_DVT(DVT2, DVT3)
    call delete_link_DVT(DVT4, DVT3)

  end subroutine flip31

  subroutine merge_incident_DT(DVT, oldDT1, oldDT2, newDT)

    type(delaunay_vertex_type), intent(inout), target :: DVT
    type(delaunay_triangle_type), intent(in), target :: oldDT1
    type(delaunay_triangle_type), intent(in), target :: oldDT2
    type(delaunay_triangle_type), intent(in), target :: newDT

    call DVT%incDT%replace_ptr(oldDT1, newDT)
    call DVT%incDT%delete_ptr(oldDT2)

  end subroutine merge_incident_DT

  subroutine split_incident_DT(DVT, oldDT, newDT1, newDT2)

    type(delaunay_vertex_type), intent(inout), target :: DVT
    type(delaunay_triangle_type), intent(in), target :: oldDT
    type(delaunay_triangle_type), intent(in), target :: newDT1
    type(delaunay_triangle_type), intent(in), target :: newDT2

    call DVT%incDT%replace_ptr(oldDT, newDT1)
    call DVT%incDT%insert_ptr_after(newDT1, newDT2)

  end subroutine split_incident_DT

  subroutine add_incident_DT(DVT, DT)

    type(delaunay_vertex_type), intent(inout) :: DVT
    type(delaunay_triangle_type), intent(in) :: DT

    call DVT%incDT%append_ptr(DT)

  end subroutine add_incident_DT

  subroutine delete_link_DVT(DVT, delDVT)

    type(delaunay_vertex_type), intent(inout), target :: DVT
    type(delaunay_vertex_type), intent(in), target :: delDVT

    call DVT%linkDVT%delete_ptr(delDVT)

  end subroutine delete_link_DVT

  subroutine add_link_DVT(DVT, addDVT, DVT1, DVT2)

    type(delaunay_vertex_type), intent(inout) :: DVT
    type(delaunay_vertex_type), intent(in) :: addDVT
    type(delaunay_vertex_type), intent(in), optional :: DVT1
    type(delaunay_vertex_type), intent(in), optional :: DVT2

    if (present(DVT1)) then
      call DVT%linkDVT%insert_ptr_after(DVT1, addDVT)
      ! TODO: Check if next link DVT is DVT2
    else
      call DVT%linkDVT%append_ptr(addDVT)
    end if

  end subroutine add_link_DVT

  subroutine record_included_DVT(DVT, DT1, DT2)

    type(delaunay_vertex_type), intent(inout) :: DVT
    type(delaunay_triangle_type), intent(inout), target :: DT1
    type(delaunay_triangle_type), intent(inout), target, optional :: DT2

    call DT1%incDVT%append_ptr(DVT)
    DVT%cntDT1 => DT1
    DVT%stub_item1 => DT1%incDVT%last_item
    if (present(DT2)) then
      call DT2%incDVT%append_ptr(DVT)
      DVT%cntDT2 => DT2
      DVT%stub_item2 => DT2%incDVT%last_item
    else
      DVT%cntDT2 => null()
      DVT%stub_item2 => null()
    end if

  end subroutine record_included_DVT

  subroutine delete_included_DVT(DVT, DT1, DT2)

    type(delaunay_vertex_type), intent(inout) :: DVT
    type(delaunay_triangle_type), intent(inout) :: DT1
    type(delaunay_triangle_type), intent(inout), optional :: DT2

    if (associated(DVT%stub_item1)) call DT1%incDVT%remove_item(DVT%stub_item1)
    if (present(DT2) .and. associated(DVT%stub_item2)) call DT2%incDVT%remove_item(DVT%stub_item2)
    DVT%cntDT1 => null()
    DVT%cntDT2 => null()

  end subroutine delete_included_DVT

  subroutine record_obsolete_DT(DT)

    type(delaunay_triangle_type), intent(in) :: DT

    call obsDT_list%append_ptr(DT)

  end subroutine record_obsolete_DT

  subroutine delete_obsolete_DT()

    type(linked_list_iterator_type) iterator

    iterator = linked_list_iterator(obsDT_list)
    do while (.not. iterator%ended())
      select type (DT => iterator%value)
      type is (delaunay_triangle_type)
        call global_DT_list%remove_item(DT%stub)
        call obsDT_list%delete_and_next(iterator)
      end select
    end do

  end subroutine delete_obsolete_DT

  subroutine record_temporal_DT(DT)

    type(delaunay_triangle_type), intent(in) :: DT

    call tmpDT_list%append_ptr(DT)

  end subroutine record_temporal_DT

  subroutine delete_temporal_DT()

    type(linked_list_iterator_type) iterator

    iterator = linked_list_iterator(tmpDT_list)
    do while (.not. iterator%ended())
      select type (DT => iterator%value)
      type is (delaunay_triangle_type)
        call global_DT_list%remove_item(DT%stub)
        call tmpDT_list%delete_and_next(iterator)
      end select
    end do

  end subroutine delete_temporal_DT

  subroutine calc_circumcircle(DT, center)

    type(delaunay_triangle_type), intent(inout) :: DT
    type(voronoi_vertex_type), intent(inout), target :: center

    type(delaunay_vertex_type), pointer :: DVT1, DVT2, DVT3
    real(8) e2(3), e3(3), n(3), l, tmp
    integer i

    DVT1 => get_DVT(DT%DVT, 1)
    DVT2 => get_DVT(DT%DVT, 2)
    DVT3 => get_DVT(DT%DVT, 3)
    e3 = [DVT2%x - DVT1%x, DVT2%y - DVT1%y, DVT2%z - DVT1%z]
    e2 = [DVT3%x - DVT1%x, DVT3%y - DVT1%y, DVT3%z - DVT1%z]
    n  = cross_product(e3, e2)
    l  = sqrt(sum(n**2))
    if (l < eps) then
      call log_error('Vertices ' // trim(to_string(DVT1%id)) // ', ' // trim(to_string(DVT2%id)) // ', ' // trim(to_string(DVT3%id)) // ' are collinear!')
    end if
    center%x = n(1) / l
    center%y = n(2) / l
    center%z = n(3) / l
    call inverse_cartesian_transform(center)
    DT%center => center

  end subroutine calc_circumcircle

  ! ----------------------------------------------------------------------------
  !                        Type bounding procedures
  ! ----------------------------------------------------------------------------

  subroutine delaunay_vertex_init(this, id, lon, lat, x, y, z)

    class(delaunay_vertex_type), intent(inout) :: this
    integer, intent(in) :: id
    real(8), intent(in), optional :: lon
    real(8), intent(in), optional :: lat
    real(8), intent(in), optional :: x
    real(8), intent(in), optional :: y
    real(8), intent(in), optional :: z

    this%id = id
    if (present(lon) .and. present(lat)) then
      this%lon = lon
      this%lat = lat
      call cartesian_transform(this)
    else if (present(x) .and. present(y) .and. present(z)) then
      this%x = x
      this%y = y
      this%z = z
      call inverse_cartesian_transform(this)
    end if

  end subroutine delaunay_vertex_init

  subroutine delaunay_vertex_print(this)

    class(delaunay_vertex_type), intent(in) :: this

    type(delaunay_triangle_type), pointer :: DT
    type(delaunay_vertex_type), pointer :: DVT
    integer i

    write(*, *) 'Delaunay vertex (' // trim(to_string(this%id)) // '):'
    write(*, *) '  LON: ', this%lon * deg, ' LAT: ', this%lat * deg
    write(*, *) '  X: ', this%x, ' Y: ', this%y, ' Z: ', this%z
    write(*, '(A)', advance='no') '  Incident triangles: '
    do i = 1, this%incDT%size
      DT => get_DT(this%incDT, i)
      write(*, '(A)', advance='no') trim(to_string(DT%id)) // ', '
    end do
    write(*, *)
    write(*, '(A)', advance='no') '  Link vertices: '
    do i = 1, this%linkDVT%size
      DVT => get_DVT(this%linkDVT, i)
      write(*, '(A)', advance='no') trim(to_string(DVT%id)) // ', '
    end do
    write(*, *)
    write(*, '(A)', advance='no') '  Include triangle: '
    if (associated(this%cntDT1)) then
      write(*, '(A)', advance='no') trim(to_string(this%cntDT1%id)) // ', '
    end if
    if (associated(this%cntDT2)) then
      write(*, '(A)', advance='no') trim(to_string(this%cntDT2%id)) // ', '
    end if
    write(*, *)

  end subroutine delaunay_vertex_print

  subroutine delaunay_triangle_init(this, id)

    class(delaunay_triangle_type), intent(inout) :: this
    integer, intent(in) :: id

    this%id = id

  end subroutine delaunay_triangle_init

  subroutine delaunay_triangle_add_DVT(this, DVT1, DVT2, DVT3)

    class(delaunay_triangle_type), intent(inout) :: this
    type(delaunay_vertex_type), intent(in) :: DVT1
    type(delaunay_vertex_type), intent(in) :: DVT2
    type(delaunay_vertex_type), intent(in) :: DVT3

    call this%DVT%append_ptr(DVT1)
    call this%DVT%append_ptr(DVT2)
    call this%DVT%append_ptr(DVT3)

  end subroutine delaunay_triangle_add_DVT

  subroutine delaunay_triangle_add_adjDT(this, adjDT1, adjDT2, adjDT3)

    class(delaunay_triangle_type), intent(inout) :: this
    type(delaunay_triangle_type), intent(in) :: adjDT1
    type(delaunay_triangle_type), intent(in) :: adjDT2
    type(delaunay_triangle_type), intent(in) :: adjDT3

    call this%adjDT%append_ptr(adjDT1)
    call this%adjDT%append_ptr(adjDT2)
    call this%adjDT%append_ptr(adjDT3)

  end subroutine delaunay_triangle_add_adjDT

  subroutine delaunay_triangle_print(this, details)

    class(delaunay_triangle_type), intent(in) :: this
    logical, intent(in), optional :: details

    type(delaunay_vertex_type), pointer :: DVT
    type(delaunay_triangle_type), pointer :: DT
    integer i

    write(*, *) 'Delaunay triangle (' // trim(to_string(this%id)) // '): '
    if (this%DVT%size > 0) then
      write(*, '(A)', advance='no') '   Vertices: '
      do i = 1, this%DVT%size
        DVT => get_DVT(this%DVT, i)
        write(*, '(A)', advance='no') trim(to_string(DVT%id)) // ', '
      end do
      write(*, *)
      if (present(details) .and. details) then
        do i = 1, this%DVT%size
          DVT => get_DVT(this%DVT, i)
          call DVT%print()
        end do
        write(*, *) 'For NCL:'
        write(*, '(A)', advance='no') '  Lon: (/'
        do i = 1, this%DVT%size
          DVT => get_DVT(this%DVT, i)
          write(*, '(F10.5, ",")', advance='no') DVT%lon * deg
        end do
        DVT => get_DVT(this%DVT, 1)
        write(*, '(F10.5, A)') DVT%lon * deg, '/)'
        write(*, '(A)', advance='no') '  Lat: (/'
        do i = 1, this%DVT%size
          DVT => get_DVT(this%DVT, i)
          write(*, '(F10.5, ",")', advance='no') DVT%lat * deg
        end do
        DVT => get_DVT(this%DVT, 1)
        write(*, '(F10.5, A)') DVT%lat * deg, '/)'
      end if
    else
      write(*, *) '  Vertices: Not connected with vertices yet.'
    end if
    if (this%adjDT%size > 0) then
      write(*, '(A)', advance='no') '   Neightbor triangles: '
      do i = 1, this%adjDT%size
        DT => get_DT(this%adjDT, i)
        write(*, '(A)', advance='no') trim(to_string(DT%id)) // ', '
      end do
      write(*, *)
    else
      write(*, *) '  Neightbor triangles: Not connected with triangles yet.'
    end if

  end subroutine delaunay_triangle_print

  subroutine delaunay_edge_init(this, id)

    class(delaunay_edge_type), intent(inout) :: this
    integer, intent(in) :: id

    this%id = id

  end subroutine delaunay_edge_init

  subroutine voronoi_vertex_init(this, id)

    class(voronoi_vertex_type), intent(inout) :: this
    integer, intent(in) :: id

    this%id = id
    this%VE = array(3)

  end subroutine voronoi_vertex_init

  subroutine voronoi_vertex_print(this)

    class(voronoi_vertex_type), intent(in) :: this

    type(voronoi_edge_type), pointer :: VE
    type(voronoi_cell_type), pointer :: VC
    integer i

    write(*, *) 'Voronoi vertex (' // trim(to_string(this%id)) // '): '
    if (this%VE%size > 0) then
      write(*, '(A)', advance='no') '   Edges: '
      do i = 1, this%VE%size
        VE => get_VE(this%VE, i)
        write(*, '(A)', advance='no') trim(to_string(VE%id)) // ', '
      end do
      write(*, *)
    else
      write(*, *) '  Edges: Not connected with edges yet.'
    end if
    if (this%VC%size > 0) then
      write(*, '(A)', advance='no') '   Cells: '
      do i = 1, this%VC%size
        VC => get_VC(this%VC, i)
        write(*, '(A)', advance='no') trim(to_string(VC%id)) // ', '
      end do
      write(*, *)
    else
      write(*, *) '  Cells: Not connected with cells yet.'
    end if

  end subroutine voronoi_vertex_print

  subroutine voronoi_edge_init(this, id)

    class(voronoi_edge_type), intent(inout) :: this
    integer, intent(in) :: id

    this%id = id

  end subroutine voronoi_edge_init

  subroutine voronoi_edge_print(this)

    class(voronoi_edge_type), intent(in) :: this

    write(*, *) 'Voronoi edge (' // trim(to_string(this%id)) // '): '
    write(*, '(A)', advance='no') '   Vertices: '
    if (associated(this%VVT1)) then
      write(*, '(A)', advance='no') trim(to_string(this%VVT1%id)) // ', '
    else
      write(*, '(A)', advance='no') 'NaN, '
    end if
    if (associated(this%VVT2)) then
      write(*, '(A)', advance='no') trim(to_string(this%VVT2%id)) // ', '
    else
      write(*, '(A)', advance='no') 'NaN'
    end if
    write(*, *)
    write(*, '(A)', advance='no') '   Cells: '
    if (associated(this%VC1)) then
      write(*, '(A)', advance='no') trim(to_string(this%VC1%id)) // ', '
    else
      write(*, '(A)', advance='no') 'NaN, '
    end if
    if (associated(this%VC2)) then
      write(*, '(A)', advance='no') trim(to_string(this%VC2%id)) // ', '
    else
      write(*, '(A)', advance='no') 'NaN'
    end if
    write(*, *)

  end subroutine voronoi_edge_print

  subroutine voronoi_cell_init(this, id)

    class(voronoi_cell_type), intent(inout) :: this
    integer, intent(in) :: id

    this%id = id
    this%VVT = array(6)
    this%VE = array(6)
    this%adjVC = array(6)

  end subroutine voronoi_cell_init

  subroutine voronoi_cell_print(this)

    class(voronoi_cell_type), intent(in) :: this

    type(voronoi_vertex_type), pointer :: VVT
    type(voronoi_cell_type), pointer :: VC
    type(voronoi_edge_type), pointer :: VE
    integer i

    write(*, *) 'Voronoi cell (' // trim(to_string(this%id)) // '): '
    if (this%VVT%size > 0) then
      write(*, '(A)', advance='no') '   Vertices: '
      do i = 1, this%VVT%size
        VVT => get_VVT(this%VVT, i)
        write(*, '(A)', advance='no') trim(to_string(VVT%id)) // ', '
      end do
      write(*, *)
    else
      write(*, *) '  Vertices: Not connected with vertices yet.'
    end if
    if (this%adjVC%size > 0) then
      write(*, '(A)', advance='no') '   Neightbor cells: '
      do i = 1, this%adjVC%size
        VC => get_VC(this%adjVC, i)
        write(*, '(A)', advance='no') trim(to_string(VC%id)) // ', '
      end do
      write(*, *)
    else
      write(*, *) '  Neightbor cells: Not connected with cells yet.'
    end if
    if (this%VE%size > 0) then
      write(*, '(A)', advance='no') '   Edges: '
      do i = 1, this%VE%size
        VE => get_VE(this%VE, i)
        write(*, '(A)', advance='no') trim(to_string(VE%id)) // ', '
      end do
      write(*, *)
    else
      write(*, *) '  Edges: Not connected with edges yet.'
    end if

  end subroutine voronoi_cell_print

  function get_DVT_from_list(list, index) result(res)

    type(linked_list_type), intent(in) :: list
    integer, intent(in) :: index
    type(delaunay_vertex_type), pointer :: res

    select type (val => list%value_at(index))
    type is (delaunay_vertex_type)
      res => val
    class default
      res => null()
    end select

  end function get_DVT_from_list

  function get_DVT_from_array(array, index) result(res)

    type(array_type), intent(in) :: array
    integer, intent(in) :: index
    type(delaunay_vertex_type), pointer :: res

    select type (val => array%value_at(index))
    type is (delaunay_vertex_type)
      res => val
    class default
      res => null()
    end select

  end function get_DVT_from_array

  function get_DVT_from_item(item) result(res)

    type(linked_list_item_type), intent(in) :: item
    type(delaunay_vertex_type), pointer :: res

    select type (val => item%value)
    type is (delaunay_vertex_type)
      res => val
    class default
      res => null()
    end select

  end function get_DVT_from_item

  function get_DVT_from_value(value) result(res)

    class(*), intent(in), pointer :: value
    type(delaunay_vertex_type), pointer :: res

    select type (value)
    type is (delaunay_vertex_type)
      res => value
    class default
      res => null()
    end select

  end function get_DVT_from_value

  subroutine create_DT(DT)

    type(delaunay_triangle_type), intent(out), pointer :: DT

    type(delaunay_triangle_type) local_DT

    ! Let global_DT_list manages memory to
    id_DT = id_DT + 1
    call global_DT_list%append(local_DT)
    select type (val => global_DT_list%last_value())
    type is (delaunay_triangle_type)
      DT => val
      call DT%init(id=id_DT)
      DT%stub => global_DT_list%last_item
    end select

  end subroutine create_DT

  function get_DT_from_list(list, index) result(res)

    type(linked_list_type), intent(in) :: list
    integer, intent(in) :: index
    type(delaunay_triangle_type), pointer :: res

    select type (val => list%value_at(index))
    type is (delaunay_triangle_type)
      res => val
    class default
      res => null()
    end select

  end function get_DT_from_list

  function get_DT_from_array(array, index) result(res)

    type(array_type), intent(in) :: array
    integer, intent(in) :: index
    type(delaunay_triangle_type), pointer :: res

    select type (val => array%value_at(index))
    type is (delaunay_triangle_type)
      res => val
    class default
      res => null()
    end select

  end function get_DT_from_array

  function get_DT_from_value(value) result(res)

    class(*), intent(in), target :: value
    type(delaunay_triangle_type), pointer :: res

    select type (val => value)
    type is (delaunay_triangle_type)
      res => val
    class default
      res => null()
    end select

  end function get_DT_from_value

  function get_DE_from_array(array, index) result(res)

    type(array_type), intent(in) :: array
    integer, intent(in) :: index
    type(delaunay_edge_type), pointer :: res

    select type (val => array%value_at(index))
    type is (delaunay_edge_type)
      res => val
    class default
      res => null()
    end select

  end function get_DE_from_array

  function get_VVT_from_array(array, index) result(res)

    type(array_type), intent(in) :: array
    integer, intent(in) :: Index
    type(voronoi_vertex_type), pointer :: res

    select type (val => array%value_at(index))
    type is (voronoi_vertex_type)
      res => val
    class default
      res => null()
    end select

  end function get_VVT_from_array

  function get_VVT_from_value(value) result(res)

    class(*), intent(in), target :: value
    type(voronoi_vertex_type), pointer :: res

    select type (val => value)
    type is (voronoi_vertex_type)
      res => val
    class default
      res => null()
    end select

  end function get_VVT_from_value

  function get_VC_from_array(array, index) result(res)

    type(array_type), intent(in) :: array
    integer, intent(in) :: index
    type(voronoi_cell_type), pointer :: res

    select type (val => array%value_at(index))
    type is (voronoi_cell_type)
      res => val
    class default
      res => null()
    end select

  end function get_VC_from_array

  function get_VE_from_array(array, index) result(res)

    type(array_type), intent(in) :: array
    integer, intent(in) :: index
    type(voronoi_edge_type), pointer :: res

    select type (val => array%value_at(index))
    type is (voronoi_edge_type)
      res => val
    class default
      res => null()
    end select

  end function get_VE_from_array

  function get_VE_from_value(value) result(res)

    class(*), intent(in), target :: value
    type(voronoi_edge_type), pointer :: res

    select type (val => value)
    type is (voronoi_edge_type)
      res => val
    class default
      res => null()
    end select

  end function get_VE_from_value

  ! ----------------------------------------------------------------------------
  !                            Geometric predicates
  ! ----------------------------------------------------------------------------

  ! ************************************************************************** !
  ! orient subroutine                                                          !
  ! Purpose:                                                                   !
  !   This is the one of two basic geometric tests to give the relation        !
  !   between a point (x0,y0,z0) and a plane containing (x1,y1,z1),            !
  !   (x2,y2,z2), and (0,0,0). There are three types of relation, i.e. left,   !
  !   right, and on, where left is defined relative to an observer at          !
  !   (x1,y1,z1) facing (x2,y2,z2).                                            !
  ! Return value:                                                              !
  !   ORIENT_LEFT                                                              !
  !   ORIENT_RIGHT                                                             !
  !   ORIENT_ON                                                                !
  ! ************************************************************************** !

  integer function orient1(x1, y1, z1, x2, y2, z2, x0, y0, z0) result(res)

    real(8), intent(in) :: x1, y1, z1, x2, y2, z2, x0, y0, z0

    real(8) det

    det = x0 * (y1 * z2 - y2 * z1) - y0 * (x1 * z2 - x2 * z1) + z0 * (x1 * y2 - x2 * y1)

    if (det > eps) then
      res = ORIENT_LEFT
    else if (-det > eps) then
      res = ORIENT_RIGHT
    else
      res = ORIENT_ON
    end if

  end function orient1

  integer function orient2(p1, p2, p0) result(res)

    class(point_type), intent(in) :: p1, p2, p0

    res = orient1(p1%x, p1%y, p1%z, p2%x, p2%y, p2%z, p0%x, p0%y, p0%z)

  end function orient2

  ! ************************************************************************** !
  ! in_triangle(_relaxed) subroutine                                           !
  ! Purpose:                                                                   !
  !   This is a derived geometric test from "Orient". It gives the relation    !
  !   between a point (DVT) and a triangle (DT). There are four types of       !
  !   relation, i.e. inside, outside, on one edge and on one vertex.           !
  !   For relaxed InTriangle, the overlap with virtual DVTs situation will     !
  !   be checked.                                                              !
  ! Return value:                                                              !
  !   INSIDE_TRIANGLE                                                          !
  !   OUTSIDE_TRIANGLE                                                         !
  !    1 ~  3 for on one edge, where the value indicates the edge index        !
  !   -1 ~ -3 for on one vertex, where the abs(value) indicates the vertex     !
  !           index                                                            !
  ! ************************************************************************** !

  integer function in_triangle(DT, DVT) result(res)

    type(delaunay_triangle_type), intent(in) :: DT
    type(delaunay_vertex_type), intent(in) :: DVT

    real(8) x(3), y(3), z(3), x0, y0, z0
    integer i, k, ret, on_plane(2)
    type(delaunay_vertex_type), pointer :: DVT1

    ! Copy for short-hand
    do i = 1, DT%DVT%size
      DVT1 => get_DVT(DT%DVT, i)
      x(i) = DVT1%x; y(i) = DVT1%y; z(i) = DVT1%z
    end do
    x0 = DVT%x; y0 = DVT%y; z0 = DVT%z

    k = 0
    do i = 1, 3
      ret = orient(x(i), y(i), z(i), x(ip1(i)), y(ip1(i)), z(ip1(i)), x0, y0, z0)
      if (ret == ORIENT_RIGHT) then
        res = OUTSIDE_TRIANGLE
        return
      else if (ret == ORIENT_ON) then
        k = k + 1; on_plane(k) = im1(i)
      end if
    end do

    if (k == 0) then
      res = INSIDE_TRIANGLE
    else if (k == 1) then ! on the edge
      res = on_plane(k)
    else if (k == 2) then ! on the vertex
      if (on_plane(1) == 1 .and. on_plane(2) == 2) then
        res = -3
      else if (on_plane(1) == 2 .and. on_plane(2) == 1) then
        res = -3
      else if (on_plane(1) == 2 .and. on_plane(2) == 3) then
        res = -1
      else if (on_plane(1) == 3 .and. on_plane(2) == 2) then
        res = -1
      else if (on_plane(1) == 1 .and. on_plane(2) == 3) then
        res = -2
      else if (on_plane(1) == 3 .and. on_plane(2) == 1) then
        res = -2
      end if
    end if

  end function in_triangle

  integer function in_triangle_relaxed(DT, DVT) result(res)

    type(delaunay_triangle_type), intent(in) :: DT
    type(delaunay_vertex_type), intent(in) :: DVT

    real(8), parameter :: eps = 1.0e-4
    real(8) x(3), y(3), z(3), x0, y0, z0
    integer i, j, k, ret, on_plane(2)
    type(delaunay_vertex_type), pointer :: DVT1

    ! Copy for short-hand
    do i = 1, DT%DVT%size
      DVT1 => get_DVT(DT%DVT, i)
      if (DVT1%id < 0 .and. abs(DVT1%lon - DVT%lon) < eps .and. abs(DVT1%lat - DVT%lat) < eps) then
        res = -i
        return
      end if
      x(i) = DVT1%x; y(i) = DVT1%y; z(i) = DVT1%z
    end do
    x0 = DVT%x; y0 = DVT%y; z0 = DVT%z

    k = 0
    do i = 1, 3
      j = ip1(i)
      ret = orient(x(i), y(i), z(i), x(j), y(j), z(j), x0, y0, z0)
      if (ret == ORIENT_RIGHT) then
        res = OUTSIDE_TRIANGLE
        return
      else if (ret == ORIENT_ON) then
        k = k + 1; on_plane(k) = im1(i)
      end if
    end do

    if (k == 0) then
      res = INSIDE_TRIANGLE
    else if (k == 1) then ! on the edge
      res = on_plane(k)
    else if (k == 2) then ! on the vertex
      if (on_plane(1) == 1 .and. on_plane(2) == 2) then
        res = -3
      else if (on_plane(1) == 2 .and. on_plane(2) == 1) then
        res = -3
      else if (on_plane(1) == 2 .and. on_plane(2) == 3) then
        res = -1
      else if (on_plane(1) == 3 .and. on_plane(2) == 2) then
        res = -1
      else if (on_plane(1) == 1 .and. on_plane(2) == 3) then
        res = -2
      else if (on_plane(1) == 3 .and. on_plane(2) == 1) then
        res = -2
      end if
    end if

  end function in_triangle_relaxed

  ! ************************************************************************** !
  ! in_circle                                                                  !
  ! Purpose:                                                                   !
  !   This is the one of two basic geometric tests to check whether or not a   !
  !   point (DVT) is inside the circumcircle (here is spherical circle) of a   !
  !   Delaunay triangle (DT).                                                  !
  !   In Cartesian coordinate system, this test turns into checking whether    !
  !   or not DVT is above the plane defined by the three vertices of DT.       !
  ! ************************************************************************** !

  integer function in_circle1(DT, DVT) result(res)

    type(delaunay_triangle_type), intent(in) :: DT
    type(delaunay_vertex_type), intent(in) :: DVT

    real(8) dx(3), dy(3), dz(3), x0, y0, z0, det
    integer i
    type(delaunay_vertex_type), pointer :: DVT1

    x0 = DVT%x
    y0 = DVT%y
    z0 = DVT%z
    do i = 1, 3
      DVT1 => get_DVT(DT%DVT, i)
      dx(i) = DVT1%x - x0
      dy(i) = DVT1%y - y0
      dz(i) = DVT1%z - z0
    end do

    det = dx(3)*(dy(2)*dz(1) - dy(1)*dz(2)) &
        - dy(3)*(dx(2)*dz(1) - dx(1)*dz(2)) &
        + dz(3)*(dx(2)*dy(1) - dx(1)*dy(2))

    if (det > eps) then
      res = INSIDE_CIRCLE
    else if (-det > eps) then
      res = OUTSIDE_CIRCLE
    else
      res = ON_CIRCLE
    end if

  end function in_circle1

  integer function in_circle2(x1, y1, z1, x2, y2, z2, x3, y3, z3, x0, y0, z0) result(res)

    real(8), intent(in) :: x1, y1, z1, x2, y2, z2, x3, y3, z3
    real(8), intent(in) :: x0, y0, z0

    real(8) dx1, dy1, dz1, dx2, dy2, dz2, dx3, dy3, dz3, det

    dx1 = x1 - x0; dy1 = y1 - y0; dz1 = z1 - z0
    dx2 = x2 - x0; dy2 = y2 - y0; dz2 = z2 - z0
    dx3 = x3 - x0; dy3 = y3 - y0; dz3 = z3 - z0

    det = dx3*(dy2*dz1 - dy1*dz2) &
        - dy3*(dx2*dz1 - dx1*dz2) &
        + dz3*(dx2*dy1 - dx1*dy2)

    if (det > eps) then
      res = INSIDE_CIRCLE
    else if (-det > eps) then
      res = OUTSIDE_CIRCLE
    else
      res = ON_CIRCLE
    end if

  end function in_circle2

  integer function in_circle3(p1, p2, p3, p0) result(res)

    class(point_type), intent(in) :: p1, p2, p3, p0

    res = in_circle2(p1%x, p1%y, p1%z, p2%x, p2%y, p2%z, p3%x, p3%y, p3%z, p0%x, p0%y, p0%z)

  end function in_circle3

end module delaunay_voronoi_mod
