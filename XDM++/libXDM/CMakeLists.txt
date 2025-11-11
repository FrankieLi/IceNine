# XDM++ Library
set( FileList 3dMath.cpp )
set( FileList ${FileList} BBox.cpp )
set( FileList ${FileList} CrystalStructure.cpp )
set( FileList ${FileList} PeakFile.cpp )
set( FileList ${FileList} Debug.cpp )
set( FileList ${FileList} Geometry.cpp )
set( FileList ${FileList} MicGrid.cpp )
set( FileList ${FileList} Quaternion.cpp )
set( FileList ${FileList} Sampling.cpp )
set( FileList ${FileList} UFF.cpp )
set( FileList ${FileList} Voxel.cpp )
set( FileList ${FileList} XDMVoxel.cpp )
set( FileList ${FileList} Symmetry.cpp )
set( FileList ${FileList} Parser.cpp )
set( FileList ${FileList} RandomMatrix.cpp )

add_library( XDM++ STATIC ${FileList} )

# Include directories for XDM++ library
target_include_directories(XDM++ PUBLIC
    ${CMAKE_CURRENT_SOURCE_DIR}
    ${Boost_INCLUDE_DIRS}
    ${MPI_CXX_INCLUDE_DIRS}
)

# Link Eigen3 if available
if(TARGET Eigen3::Eigen)
    target_link_libraries(XDM++ PUBLIC Eigen3::Eigen)
endif()
