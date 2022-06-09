#LIBRARIES

#EIGEN3
IF(DEFINED ENV{EIGEN3_INCLUDE_DIR})
  MESSAGE(STATUS "Found Custom EIGEN3 @ $ENV{EIGEN3_INCLUDE_DIR}")
  INCLUDE_DIRECTORIES($ENV{EIGEN3_INCLUDE_DIR})
ELSE()
  FIND_PACKAGE(Eigen3 QUIET)
  IF(EIGEN3_FOUND)
    INCLUDE_DIRECTORIES(${EIGEN3_INCLUDE_DIR})
    MESSAGE(STATUS "Found EIGEN3 @ ${EIGEN3_INCLUDE_DIR}")
  ELSE(EIGEN3_FOUND)
    MESSAGE(WARNING "Cannot find EIGEN3, using local version!")
    SET(EIGEN3_INCLUDE_DIR ${PROJECT_SOURCE_DIR}/lib/eigen3)
    INCLUDE_DIRECTORIES(${EIGEN3_INCLUDE_DIR})
  ENDIF(EIGEN3_FOUND)
ENDIF()

#LIBIGL
option(LIBIGL_WITH_OPENGL            "Use OpenGL"         ON)
option(LIBIGL_WITH_OPENGL_GLFW       "Use GLFW"           ON)
option(LIBIGL_WITH_PNG               "Use PNG"            ON)

if(LIBIGL_WITH_PNG)
  # png/ module is anomalous because it also depends on opengl it really should
  # be moved into the opengl/ directory and namespace ...
  if(TARGET igl_opengl)
    if(NOT TARGET stb_image)
      igl_download_stb()
      add_subdirectory(${LIBIGL_EXTERNAL}/stb stb_image)
    endif()
    compile_igl_module("png" "")
    target_link_libraries(igl_png ${IGL_SCOPE} igl_stb_image igl_opengl)
  endif()
endif()

find_package(LIBIGL REQUIRED QUIET)

#OMPL
find_package(ompl)
INCLUDE_DIRECTORIES(${OMPL_INCLUDE_DIRS})
LIST(APPEND ALL_LIBRARIES ${OMPL_LIBRARIES})

