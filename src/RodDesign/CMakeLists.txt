cmake_minimum_required(VERSION 3.1.0)
project(RodDesignGUI)

IF(NOT CMAKE_BUILD_TYPE)
  SET(CMAKE_BUILD_TYPE Release)
ENDIF()

MESSAGE("Build type: " ${CMAKE_BUILD_TYPE})


set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS}  -Wall  -O3 -march=native ")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall   -O3 -march=native")

# Check C++11 or C++0x support
set (CMAKE_CXX_STANDARD 11)

set(CMAKE_INCLUDE_CURRENT_DIR ON)
set(CMAKE_AUTOMOC ON)
set(CMAKE_AUTORCC ON)
set(CMAKE_AUTOUIC ON)

find_package(Eigen3 3.1.0 REQUIRED)
find_package(OpenGL REQUIRED)
find_package(GLUT REQUIRED)
find_package(Qt5 COMPONENTS Core REQUIRED)
find_package(Qt5 COMPONENTS Widgets REQUIRED)
find_package(Qt5 COMPONENTS OpenGL REQUIRED)

include_directories(
  ${PROJECT_SOURCE_DIR}/include
  ${EIGEN3_INCLUDE_DIRS}
  ${OPENGL_INCLUDE_DIRS}
  ${GLUT_INCLUDE_DIRS}
)

add_executable(${PROJECT_NAME}
    main.cpp
    src/AABB_data_structure.cpp
    src/DisplayFunctions.cpp
    src/EditTools.cpp
    src/MaterialSettingDialog.cpp
    src/MomentReductionMap.cpp
    src/RodNetwork.cpp
    src/Setting.cpp
    src/SimuRodDisplay.cpp
    src/basic_geometric_calculation.cpp
    src/curve_design_display.cpp
    src/curve_network.cpp
    src/geometric_data_structure.cpp
    src/loop_subdivision.cpp
    src/mainwindow.cpp
    include/EditTools.h
    include/MaterialSettingDialog.h
    include/SimuRodDisplay.h
    include/curve_design_display.h
    include/mainwindow.h
)

set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${PROJECT_SOURCE_DIR}/)
target_link_libraries(${PROJECT_NAME} 
                        ${EIGEN3_LIBS}
                        ${OPENGL_LIBRARIES}
                        ${GLUT_LIBRARIES}
                        Qt5::Core
                        Qt5::Widgets
                        Qt5::OpenGL)
