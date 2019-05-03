cmake_minimum_required (VERSION 3.0)

project (fishman LANGUAGES Fortran)

if (CMAKE_Fortran_COMPILER_ID MATCHES Intel)
  set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -Ofast")
  if (CMAKE_BUILD_TYPE MATCHES Debug)
    set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -g -gdwarf -Og -traceback -check all")
  endif ()
elseif (CMAKE_Fortran_COMPILER_ID MATCHES GNU)
  if (APPLE)
    set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -ffree-line-length-none")
  else ()
    set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -ffree-line-length-none")
  endif ()
  if (CMAKE_BUILD_TYPE MATCHES Debug)
    set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -g -gdwarf -Og -fbacktrace -fcheck=all")
  endif ()
endif ()

if (DEFINED ENV{NETCDF} AND (NOT DEFINED ENV{NETCDF_ROOT}))
  set(ENV{NETCDF_ROOT} $ENV{NETCDF})
endif ()
if (DEFINED ENV{NETCDF_ROOT})
  include_directories("$ENV{NETCDF_ROOT}/include")
  link_directories("$ENV{NETCDF_ROOT}/lib")
else ()
  find_package(PkgConfig)
  if (PKG_CONFIG_FOUND)
    pkg_search_module(NETCDF REQUIRED netcdf)
    if (NETCDF_FOUND)
      include_directories(${NETCDF_INCLUDEDIR})
      link_directories(${NETCDF_LIBDIR})
    else ()
      message(FATAL_ERROR "Pkg-config could not find netcdf library!")
    endif ()
  else ()
    message(FATAL_ERROR "Unable to find pkg-config library!")
  endif ()
endif ()

add_subdirectory(lib/container)
add_subdirectory(lib/fiona)

set(sources
  src/const_mod.F90
  src/log_mod.F90
  src/string_mod.F90
  src/random_number_mod.F90
  src/math_mod.F90
  src/sphere_geometry_mod.F90
  src/delaunay_voronoi_mod.F90
  src/mpas_mesh_mod.F90
)

add_library(fishman ${sources})
target_link_libraries(fishman netcdf netcdff fortran_container fiona)

add_executable(gen_mesh.exe src/gen_mesh.F90)
target_link_libraries(gen_mesh.exe fishman)
