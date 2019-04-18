module const_mod

  implicit none

  integer, parameter :: max_file_path_len = 256
  integer, parameter :: max_name_len = 30

  real(8), parameter :: radius = 6371220d0
  real(8), parameter :: g      = 9.80616d0
  real(8), parameter :: omega  = 7.292d-5
  real(8), parameter :: pi     = 4.0d0 * atan(1.0d0)
  real(8), parameter :: pi2    = pi * 2.0d0
  real(8), parameter :: pi05   = pi * 0.5d0
  real(8), parameter :: rad    = pi / 180.0d0
  real(8), parameter :: deg    = 180.0d0 / pi
  real(8), parameter :: eps    = epsilon(1.0d0)
  real(8), parameter :: missing_value = -1e20

end module const_mod
