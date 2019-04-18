module sphere_geometry_mod

  use const_mod
  use math_mod

  implicit none

  private

  public euler_formula
  public cartesian_transform
  public inverse_cartesian_transform
  public rotation_transform
  public inverse_rotation_transform
  public calc_distance
  public cross_product
  public norm_vector
  public calc_sphere_angle
  public calc_plane_angle
  public calc_arc_length
  public point_type

  type point_type
    real(8) lon
    real(8) lat
    real(8) x
    real(8) y
    real(8) z
  contains
    procedure :: copy_coord => point_copy_coord
  end type point_type

  interface cartesian_transform
    module procedure cartesian_transform_1
    module procedure cartesian_transform_2
  end interface cartesian_transform

  interface inverse_cartesian_transform
    module procedure inverse_cartesian_transform_1
    module procedure inverse_cartesian_transform_2
  end interface inverse_cartesian_transform

contains

  integer function euler_formula(num_cell, num_vertex, num_edge) result(res)

    integer, intent(in), optional :: num_cell
    integer, intent(in), optional :: num_vertex
    integer, intent(in), optional :: num_edge

    if (present(num_cell) .and. present(num_vertex)) then
      res = num_cell + num_vertex - 2
    else if (present(num_cell) .and. present(num_edge)) then
      res = num_edge - num_cell + 2
    else if (present(num_vertex) .and. present(num_edge)) then
      res = num_edge - num_vertex + 2
    end if

  end function euler_formula

  subroutine cartesian_transform_1(lon, lat, x, y, z)

    real(8), intent(in)  :: lon, lat
    real(8), intent(out) :: x, y, z

    real(8) cos_lat

    cos_lat = cos(lat)
    x = cos_lat * cos(lon)
    y = cos_lat * sin(lon)
    z = sin(lat)

  end subroutine cartesian_transform_1

  subroutine cartesian_transform_2(point)

    class(point_type), intent(inout) :: point

    real(8) cos_lat

    cos_lat = cos(point%lat)
    point%x = cos_lat * cos(point%lon)
    point%y = cos_lat * sin(point%lon)
    point%z = sin(point%lat)

  end subroutine cartesian_transform_2

  subroutine inverse_cartesian_transform_1(lon, lat, x, y, z)

    real(8), intent(out) :: lon, lat
    real(8), intent(in)  :: x, y, z

    lon = atan2(y, x)
    lat = asin(z)

    if (lon < 0.0d0) lon = lon + pi2

  end subroutine inverse_cartesian_transform_1

  subroutine inverse_cartesian_transform_2(point)

    class(point_type), intent(inout) :: point

    point%lon = atan2(point%y, point%x)
    point%lat = asin(point%z)

    if (point%lon < 0.0d0) point%lon = point%lon + pi2

  end subroutine inverse_cartesian_transform_2

  ! ************************************************************************** !
  ! Rotation transform                                                         !
  ! Purpose:                                                                   !
  !   Calculate the rotating transformation and its inverse of the original    !
  !   coordinate system (lon_o,lat_o) to the rotated one (lon_r, lat_r) with   !
  !   the north pole (lon_p,lat_p) defined at the original coordinate system.  !
  ! ************************************************************************** !

  subroutine rotation_transform(lon_p, lat_p, lon_o, lat_o, lon_r, lat_r)

    real(8), intent(in) :: lon_p, lat_p ! Rotated pole coordinate
    real(8), intent(in) :: lon_o, lat_o ! Original coordinate
    real(8), intent(out), optional :: lon_r, lat_r ! Rotated coordinate

    real(8) tmp1, tmp2, tmp3, dlon

    dlon = lon_o - lon_p
    if (present(lon_r)) then
        tmp1 = cos(lat_o) * sin(dlon)
        tmp2 = cos(lat_o) * sin(lat_p) * cos(dlon) - cos(lat_p) * sin(lat_o)
        lon_r = atan2(tmp1, tmp2)
        if (lon_r < 0.0d0) lon_r = PI2 + lon_r
    end if
    if (present(lat_r)) then
        tmp1 = sin(lat_o) * sin(lat_p)
        tmp2 = cos(lat_o) * cos(lat_p) * cos(dlon)
        tmp3 = tmp1 + tmp2
        tmp3 = min(max(tmp3, -1.0d0), 1.0d0)
        lat_r = asin(tmp3)
    end if

  end subroutine rotation_transform

  subroutine inverse_rotation_transform(lon_p, lat_p, lon_o, lat_o, lon_r, lat_r)

      real(8), intent(in)  :: lon_p, lat_p ! Rotated pole coordinate
      real(8), intent(out) :: lon_o, lat_o ! Original coordinate
      real(8), intent(in)  :: lon_r, lat_r ! Rotated coordinate

      real(8) sin_lon_r, cos_lon_r, sin_lat_r, cos_lat_r, sin_lat_p, cos_lat_p
      real(8) tmp1, tmp2, tmp3

      sin_lon_r = sin(lon_r)
      cos_lon_r = cos(lon_r)
      sin_lat_r = sin(lat_r)
      cos_lat_r = cos(lat_r)
      sin_lat_p = sin(lat_p)
      cos_lat_p = cos(lat_p)

      tmp1 = cos_lat_r * sin_lon_r
      tmp2 = sin_lat_r * cos_lat_p + cos_lat_r * cos_lon_r * sin_lat_p
      ! This trick is due to the inaccuracy of trigonometry calculation.
      if (abs(tmp2) < eps) tmp2 = 0.0d0
      lon_o = atan2(tmp1, tmp2)
      lon_o = lon_p + lon_o
      if (lon_o > PI2) lon_o = lon_o - PI2
      tmp1 = sin_lat_r * sin_lat_p
      tmp2 = cos_lat_r * cos_lat_p * cos_lon_r
      tmp3 = tmp1 - tmp2
      tmp3 = min(max(tmp3, -1.0d0), 1.0d0)
      lat_o = asin(tmp3)

  end subroutine inverse_rotation_transform

  real(8) function calc_distance(lon1, lat1, lon2, lat2) result(res)

    real(8), intent(in) :: lon1
    real(8), intent(in) :: lat1
    real(8), intent(in) :: lon2
    real(8), intent(in) :: lat2

    res = radius * acos(min(1.0d0, max(-1.0d0, sin(lat1) * sin(lat2) + cos(lat1) * cos(lat2) * cos(lon1 - lon2))))

  end function calc_distance

  function norm_vector(x) result(res)

    real(8), intent(in) :: x(:)
    real(8) res(size(x))

    real(8) n

    n = sqrt(sum(x * x))
    if (n /= 0) then
      res = x / n
    else
      res = x
    end if

  end function norm_vector

  ! Calculate the dihedra angle between plane AB and plane AC.

  real(8) function calc_sphere_angle(a, b, c) result(res)

    real(8), intent(in) :: a(3)
    real(8), intent(in) :: b(3)
    real(8), intent(in) :: c(3)

    real(8) nab(3) ! Normal vector of plane AB
    real(8) nac(3) ! Normal vector of plane AC

    nab = norm_vector(cross_product(a, b))
    nac = norm_vector(cross_product(a, c))
    res = acos(max(min(dot_product(nab, nac), 1.0d0), -1.0d0))

    ! Judge the cyclic direction with respect to point A to handle obtuse angle.
    if (dot_product(cross_product(b - a, c - a), a) < 0.0) res = -res

  end function calc_sphere_angle

  ! Calculate the angle between vector AB and AC on a plane.

  real(8) function calc_plane_angle(a, b, c, n) result(res)

    real(8), intent(in) :: a(:) ! Reference point
    real(8), intent(in) :: b(:) ! Point 1
    real(8), intent(in) :: c(:) ! Point 2
    real(8), intent(in) :: n(:) ! Normal vector

    real(8) ab(size(a)), ac(size(a))

    ab = norm_vector(b - a)
    ac = norm_vector(c - a)

    res = acos(max(min(dot_product(ab, ac), 1.0d0), -1.0d0))

    ! Judge the cyclic direction with respect to given normal vector to handle obtuse angle.
    if (dot_product(cross_product(ab, ac), n) < 0.0d0) res = -res

  end function calc_plane_angle

  ! Calculate the great circle arc length from A to B by assuming A and B are on the unit sphere surface.

  real(8) function calc_arc_length(a, b) result(res)

    real(8), intent(in) :: a(3)
    real(8), intent(in) :: b(3)

    res = acos(max(min(dot_product(a, b), 1.0d0), -1.0d0))

  end function calc_arc_length

  subroutine point_copy_coord(a, b)

    class(point_type), intent(inout) :: a
    class(point_type), intent(in)    :: b

    a%lon = b%lon
    a%lat = b%lat
    a%x   = b%x
    a%y   = b%y
    a%z   = b%z

  end subroutine point_copy_coord

end module sphere_geometry_mod
