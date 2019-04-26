module sphere_geometry_mod

  use const_mod
  use math_mod
  use log_mod

  implicit none

  private

  public euler_formula
  public cartesian_transform
  public inverse_cartesian_transform
  public rotation_transform
  public inverse_rotation_transform
  public calc_distance
  public calc_area
  public norm_vector
  public calc_sphere_angle
  public calc_arc_length
  public intersect
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

  interface calc_sphere_angle
    module procedure calc_sphere_angle_1
    module procedure calc_sphere_angle_2
  end interface calc_sphere_angle

  interface calc_arc_length
    module procedure calc_arc_length_1
    module procedure calc_arc_length_2
    module procedure calc_arc_length_3
  end interface calc_arc_length

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
        if (lon_r < 0.0d0) lon_r = pi2 + lon_r
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
      if (lon_o > pi2) lon_o = lon_o - pi2
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

  real(8) function calc_area(x, y, z) result(res)

    real(8), intent(in) :: x(:)
    real(8), intent(in) :: y(:)
    real(8), intent(in) :: z(:)

    integer n, im1, i, ip1
    real(8) angle

    n = size(x)
#ifndef NDEBUG
    if (n < 3) then
      call log_error('Spherical polygon number is less than 3!', __FILE__, __LINE__)
    end if
#endif
    res = 0.0
    do i = 1, n
      im1 = merge(i - 1, n, i /= 1)
      ip1 = merge(i + 1, 1, i /= n)
      angle = calc_sphere_angle([x(im1),y(im1),z(im1)], [x(i),y(i),z(i)], [x(ip1),y(ip1),z(ip1)])
      res = res + angle
    end do
    res = res - (n - 2) * pi
#ifndef NDEBUG
    if (abs(res) > 1.0 .or. res <= 0.0) then
      call log_error('Encounter bad spherical polygon to calculate area!', __FILE__, __LINE__)
    end if
#endif

  end function calc_area

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

  ! Calculate the dihedra angle between plane AB and plane BC.

  real(8) function calc_sphere_angle_1(a, b, c) result(res)

    real(8), intent(in) :: a(3)
    real(8), intent(in) :: b(3)
    real(8), intent(in) :: c(3)

    real(8) nab(3) ! Normal vector of plane AB
    real(8) nbc(3) ! Normal vector of plane BC

    nab = norm_vector(cross_product(a, b))
    nbc = norm_vector(cross_product(b, c))
    res = acos(- max(min(dot_product(nab, nbc), 1.0d0), -1.0d0))

    ! Judge the cyclic direction with respect to point A to handle obtuse angle.
    if (dot_product(cross_product(nab, nbc), a) < 0.0) res = pi2 - res

  end function calc_sphere_angle_1

  real(8) function calc_sphere_angle_2(a, b, c) result(res)

    class(point_type), intent(in) :: a
    class(point_type), intent(in) :: b
    class(point_type), intent(in) :: c

    real(8) xa(3), xb(3), xc(3)
    real(8) nab(3) ! Normal vector of plane AB
    real(8) nbc(3) ! Normal vector of plane BC

    xa = [a%x,a%y,a%z]
    xb = [b%x,b%y,b%z]
    xc = [c%x,c%y,c%z]
    nab = norm_vector(cross_product(xa, xb))
    nbc = norm_vector(cross_product(xb, xc))
    res = acos(- max(min(dot_product(nab, nbc), 1.0d0), -1.0d0))

    ! Judge the cyclic direction with respect to point A to handle obtuse angle.
    if (dot_product(cross_product(nab, nbc), xa) < 0.0) res = pi2 - res

  end function calc_sphere_angle_2

  ! Calculate the great circle arc length from A to B by assuming A and B are on the unit sphere surface.

  real(8) function calc_arc_length_1(a, b) result(res)

    real(8), intent(in) :: a(3)
    real(8), intent(in) :: b(3)

    res = acos(max(min(dot_product(a, b), 1.0d0), -1.0d0))

  end function calc_arc_length_1

  real(8) function calc_arc_length_2(a, b) result(res)

    class(point_type), intent(in) :: a
    class(point_type), intent(in) :: b

    res = acos(max(min(dot_product([a%x,a%y,a%z], [b%x,b%y,b%z]), 1.0d0), -1.0d0))

  end function calc_arc_length_2

  real(8) function calc_arc_length_3(a, b) result(res)

    class(point_type), intent(in) :: a
    real(8), intent(in) :: b(3)

    res = acos(max(min(dot_product([a%x,a%y,a%z], b), 1.0d0), -1.0d0))

  end function calc_arc_length_3

  logical function intersect(a, b, c, d, e) result(res)

    class(point_type), intent(in) :: a
    class(point_type), intent(in) :: b
    class(point_type), intent(in) :: c
    class(point_type), intent(in) :: d
    class(point_type), intent(inout) :: e

    real(8) n1(3), n2(3), v(3), r
    real(8) lon1, lon2, lat1, lat2

    n1 = cross_product([a%x,a%y,a%z], [b%x,b%y,b%z])
    n2 = cross_product([c%x,c%y,c%z], [d%x,d%y,d%z])
    v  = cross_product(n1, n2)

    r = sqrt(sum(v * v))

    if (r > eps) then
      v = v / r
    else
      res = .false.
      return
    end if
    res = .true.

    lat1 = asin(v(3))
    lat2 = -lat1
    lon1 = atan2(v(2), v(1))
    lon2 = lon1 - pi

    if (lon1 < 0.0) lon1 = lon1 + pi2
    if (lon1 > pi2) lon1 = lon1 - pi2
    if (lon2 < 0.0) lon2 = lon2 + pi2
    if (lon2 > pi2) lon2 = lon2 - pi2

    e%lon = lon1
    e%lat = lat1
    call cartesian_transform(e)

    ! Check if e is the real intersection. If not use the other one.
    if (dot_product(cross_product([a%x,a%y,a%z], [e%x,e%y,e%z]), cross_product([b%x,b%y,b%z], [e%x,e%y,e%z])) >= 0 .or. &
        dot_product(cross_product([c%x,c%y,c%z], [e%x,e%y,e%z]), cross_product([d%x,d%y,d%z], [e%x,e%y,e%z])) >= 0) then
      e%lon = lon2
      e%lat = lat2
      call cartesian_transform(e)
    end if

  end function intersect

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
