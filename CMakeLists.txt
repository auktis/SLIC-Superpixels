
PROJECT(SLICSuperpixelsDGtal)

CMAKE_MINIMUM_REQUIRED(VERSION 2.6)

SET(FILES slicSuperpixels)


FIND_PACKAGE(DGtal REQUIRED)
INCLUDE_DIRECTORIES(${DGTAL_INCLUDE_DIRS})
LINK_DIRECTORIES(${DGTAL_LIBRARY_DIRS})

FIND_PACKAGE(Boost 1.46.0 REQUIRED COMPONENTS program_options)
if ( Boost_FOUND )
  message(STATUS "Boost and boost_program_options found.")
  include_directories( ${Boost_INCLUDE_DIRS} )
      SET(DGtalToolsLibDependencies ${DGtalToolsLibDependencies} 
     ${Boost_LIBRAIRIES}  
     ${Boost_PROGRAM_OPTIONS_LIBRARY})
   SET(DGtalLibInc ${Boost_INCLUDE_DIRS})
endif( Boost_FOUND) 

if ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU")
  SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -Wextra")
endif()

FOREACH(FILE ${FILES})
   ADD_EXECUTABLE(${FILE} ${FILE}  slic )
   TARGET_LINK_LIBRARIES(${FILE}  ${DGTAL_LIBRARIES}  ${DGtalToolsLibDependencies} )
ENDFOREACH(FILE)




