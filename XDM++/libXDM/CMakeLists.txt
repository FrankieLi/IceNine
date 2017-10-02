#SET(BOOST_INCLUDEDIR "/Users/sfli/Libraries/boost-trunk/") 

#FIND_PACKAGE( Boost REQUIRED )
#include_directories(${BOOST_INCLUDEDIR})
include_directories(/opt/local/include/boost)
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



add_library( XDM++ ${FileList} )
