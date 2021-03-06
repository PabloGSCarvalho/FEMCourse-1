# @file FemSC/CMakeList.txt  -- First file to read for CMake 

cmake_minimum_required (VERSION 2.8.0)

option(USING_PZ "Whether the neoPZ package will be linked with the library" OFF)

if(USING_PZ)

IF (WIN32)
# Prevents timespec redefinition problem with Visual Studio 2015
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DHAVE_STRUCT_TIMESPEC")
ELSE (WIN32)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fPIC")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wno-narrowing")#this flag is for preventing errors when casting from double to float
ENDIF()


set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

IF (WIN32)
find_package( PZ PATHS "C:/Arquivos de Programas/PZ" REQUIRED NO_DEFAULT_PATH)
ELSE()
find_package( PZ REQUIRED )
ENDIF (WIN32)


include_directories(${PZ_INCLUDE_DIRS})
link_libraries(pz)

endif(USING_PZ)

#define se biblioteca Armadillo deve ser usada
option(USING_ARMADILLO "Biblioteca Armadillo deve ser usada neste projeto" OFF)

if(USING_ARMADILLO)
find_package( Armadillo REQUIRED NO_MODULE)
include_directories(${ARMADILLO_INCLUDE_DIRS})
endif(USING_ARMADILLO)


project (FemCourse)

add_subdirectory(headers)
add_subdirectory(headers2018)
add_subdirectory("src vazio")
add_subdirectory(src_2018)
add_subdirectory(tinyxml)


#Specify where the headers files must be found
include_directories(
	"${PROJECT_SOURCE_DIR}/headers2018"
	"${PROJECT_SOURCE_DIR}/headers"
	"${PROJECT_SOURCE_DIR}/tinyxml"
)

## inserting main definitions
#specifying the path to FemSC source code
set(FemSCSOURCE ${FemSC_SRCS} ${FemSC_tinyxml} ${FemSC_HEAD})

set(FemSOURCE2018 ${BASIC_TYPES} ${FemSC_tinyxml} ${HEADERS2018})

add_subdirectory(mainprograms)

## Some definitions for windows
if(WIN32)
	#add a define saying it is a VC compiler
	set(VC "#Will add a define VC on config.h")
	add_definitions(-DVC)
	#define _USE_MATH_DEFINES for example: M_PI constant definitions.
	add_definitions(-D_USE_MATH_DEFINES)
	#disabling VC warnings
	add_definitions(-D_SCL_SECURE_NO_WARNINGS)
	add_definitions(-D_CRT_SECURE_NO_WARNINGS)
	#disabling All warnings
	set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} /W0")
	set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} /W0")
endif()


option(USING_FAD "Whether the FAD package will be linked with the library" OFF)

if(USING_FAD)
	add_subdirectory("Fad")
endif(USING_FAD)
