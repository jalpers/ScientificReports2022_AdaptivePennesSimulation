cmake_minimum_required (VERSION 3.1)
#option (CMAKE_GNUtoMS ON)
project (live_data_heat_simulation LANGUAGES CXX CUDA)
#set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -march=native -Wall -std=c++0x -LC:\\users\\maxro\\Documents\\Respos\\LiveHeatSimulation\\cmake\\cmake-build-release -llibcmaes" )

#set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fopenmp -std=gnu++11 -I/path/to/eigen3/ -I/path/to/include/libcmaes -L/path/to/lib -o sample_code sample-code.cc -lcmaes")

set(CMAKE_CXX_STANDARD 14)

#option(USE_CUDA "Use CUDA" ON)

set(CMAKE_AUTOMOC ON)
set(CMAKE_AUTOUIC ON)
set(CMAKE_INCLUDE_CURRENT_DIR ON)


set(Qt5_DIR C:/Qt/Qt5.12.1/5.12.1/msvc2017_64/lib/cmake/Qt5)
set(VTK_DIR C:/VTK-8.2.0/out/build/x64-Release)
#set(ALGLIB_DIR C:/alglib-cmake/out/install/x64-Debug/lib/cmake/ALGLIB)


#set(CMAKE_INCLUDE_CURRENT_DIR ON)

#set(VTK_INCLUDE_DIRS C:/VTK/cmake-build-release/bin)
#include_directories(C:/VTK/cmake-build-release/bin)
#include_directories(${VTK_INCLUDE_DIRS})
#set(QT_INCLUDE_DIRS C:/QT5_MinGW/5.15.2/mingw81_64/bin/)
#include_directories(${QT_INCLUDE_DIRS})
#include_directories(C:/QT5_MinGW/5.15.2/mingw81_64/bin/)



include_directories(
	"${CMAKE_CURRENT_SOURCE_DIR}/../include"
)

# Set Local Debugger Settings (Command Arguments and Environment Variables)
#set( COMMAND_ARGUMENTS "commands" )
#set( ENVIRONMENT_VARIABLES "C:\\Qt\\5.12.1\\msvc2017_64\\bin;C:\\Users\\maxro\\Documents\\Respos\\libcmaes\\cmake-build-release\\src;C:\\VTK-8.2.0\\build\\bin\\Release";)
#set( ENVIRONMENT_VARIABLES "C:/QT5_MinGW/5.15.2/mingw81_64/bin/ ;C:/Users/maxro/Documents/Respos/libcmaes/cmake-build-release/src;C:/VTK/cmake-build-release/bin";)
#configure_file( template.vcxproj.user.in ${CMAKE_BINARY_DIR}/live_data_heat_simulation.vcxproj.user @ONLY )

set (HEADER ${HEADER}
	../include/LiveDataHeatSimulation.h
	../include/CoordinatesConverter.h
	../include/DicomHandler.h
	../include/DataVolume.h
	../include/PennesHeatTransfer.h
	../include/Matrix3D.h
	../include/ParameterOptimization.h
	../include/SimulationRessources.h
	../include/Pathfinder.h
	../include/Debugger.h
	../include/VesselMap.h
	../include/Test_cuda.h
	../include/cuda_util.h
	../include/PennesHeatTransfer_gpu.h
	../include/TestClass.h
	../include/ThermalDose.h
	../include/ExtractROI.h
	../include/NecrosisMapComputation.h
)

set (SOURCE ${SOURCE}
	../sources/main.cpp
	../sources/LiveDataHeatSimulation.cpp
	../sources/CoordinatesConverter.cpp
	../sources/DicomHandler.cpp
	../sources/DataVolume.cpp
	../sources/PennesHeatTransfer.cpp
	../sources/Matrix3D.cpp
	../sources/ParameterOptimization.cpp
	../sources/SimulationRessources.cpp
	../sources/VesselMap.cpp
	../sources/PennesHeatTransfer_gpu.cu
	../sources/Test_cuda.cu
	../sources/ThermalDose.cpp
	../sources/ExtractROI.cpp
	../sources/NecrosisMapComputation.cpp
	
)

set (UI ${UI}
	../sources/LiveDataHeatSimulation.ui
)

source_group("Source" FILES ${SOURCE})
source_group("Header" FILES ${HEADER})
source_group("UI" FILES ${UI})

#add QT5 dependencies 	
find_package(Qt5 COMPONENTS Core Gui Network WebSockets Widgets REQUIRED)

find_package(VTK REQUIRED)
find_package( OpenCV REQUIRED )
#find_package(ALGLIB REQUIRED)
find_package(OpenMP REQUIRED)

#if(USE_CUDA)
#	find_package(CUDA REQUIRED)
#endif()
find_package(CUDAToolkit REQUIRED)
#include_directories(C:/Program Files/NVIDIA GPU Computing Toolkit/CUDA/v11.6/include)
#find_package(CUDA REQUIRED)
#include_directories("${CUDA_INCLUDE_DIRS}")
#include ("${CMAKE_CURRENT_LIST_DIR}/FindEigen3.cmake")
#find_package(Eigen3 REQUIRED)
include(${VTK_USE_FILE})
include_directories(C:/opencv/include/opencv2;C:/opencv/out/build/x64-Release/opencv2 )
include_directories(C:/alglib-cmake/out/install/x64-Release/include)
#include_directories(C:/alglib-cmake/out/build/x64-Release/src/cpp/src ) 
#include_directories(C:/alglib-cmake/out/build/x64-Release/src/cpp/src ) 
#include_directories(C:/Users/maxro/Documents/Respos/libcmaes/cmake-build-release/include;C:/Users/maxro/Documents/Respos/libcmaes/include)
#find_package(OpenGL REQUIRED COMPONENTS OpenGL)
#export ALGLIB_DIR=C:/Users/maxro/Documents/Respos/alglib-cmake/out/install/Mingw64-Release/lib/cmake/ALGLIB
#link_directories(libcmaes.dll)
add_executable(live_data_heat_simulation ${SOURCE} ${HEADER} ${UI}
"ALGLIB/alglibinternal.cpp" "ALGLIB/alglibmisc.cpp" "ALGLIB/ap.cpp" "ALGLIB/dataanalysis.cpp" "ALGLIB/diffequations.cpp" "ALGLIB/fasttransforms.cpp" "ALGLIB/integration.cpp" "ALGLIB/interpolation.cpp" "ALGLIB/linalg.cpp" "ALGLIB/optimization.cpp" "ALGLIB/solvers.cpp" "ALGLIB/specialfunctions.cpp" "ALGLIB/statistics.cpp"
)
include_directories(${EIGEN_INCLUDE_DIR})
target_link_libraries(live_data_heat_simulation
	Qt5::Widgets 
	Qt5::Core
	Qt5::Network
	Qt5::WebSockets
	Qt5::Gui
	${VTK_LIBRARIES}
	${OpenCV_LIBS}
	#ALGLIB
	OpenMP::OpenMP_CXX
	#CUDA
	#${CUDA_LIBRARIES}
	CUDA::cudart
)