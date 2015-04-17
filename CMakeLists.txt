cmake_minimum_required( VERSION 2.8 )

project( lana )

set( Tinyxml_Hdr
  src/externlib/tinyxml/tinystr.h
  src/externlib/tinyxml/tinyxml.h
)

set( Tinyxml_Src
  src/externlib/tinyxml/tinystr.cpp
  src/externlib/tinyxml/tinyxml.cpp
  src/externlib/tinyxml/tinyxmlerror.cpp
  src/externlib/tinyxml/tinyxmlparser.cpp
)

set( Gml-parser_Hdr
  src/externlib/gml-parser/gml_parser.h
  src/externlib/gml-parser/gml_scanner.h
)

set( Gml-parser_Src
  src/externlib/gml-parser/gml_parser.c
  src/externlib/gml-parser/gml_scanner.c
)

set( Lana_Hdr
  src/lana.h
  src/lna.h
  src/verbose.h
  src/input/stringparser.h
  src/input/parser.h
  src/input/matchinggraph.h
  src/input/lgfparser.h
  src/input/ledaparser.h
  src/input/graphmlparser.h
  src/input/gmlparser.h
  src/input/csvparser.h
  src/input/bpparser.h
  src/input/bpcandlistparser.h
  src/input/bpblastparser.h
  src/input/bplgfparser.h
  src/input/identityparser.h
  src/input/bpidentityparser.h
  src/input/edgelistparser.h
  src/input/matchinggraphshuffle.h

  src/molecule.h
  src/bronkerboschconnected.h
  src/bronkerbosch.h
  src/product.h

  src/output/output.h
  src/output/outputsif.h
#  src/output/outputcsv.h
#  src/output/outputdot.h
#  src/output/outputgml.h
#  src/output/outputlgf.h
#  src/output/outputneato.h
#  src/output/outputsif.h
#  src/output/outputeda.h
#  src/output/outputnoa.h
#  src/output/outputparameters.h
#  src/output/outputparametersamc.h
#  src/output/outputstatamc.h
#  src/output/outputcsvamc.h

#  src/lagrangegna.h
)


set( Lana_Src
  src/lana.cpp
  src/verbose.cpp
 )



set( Common_Src
  src/verbose.cpp
)

set( GitCommand1 "symbolic-ref HEAD 2> /dev/null | cut -b 12-" )
execute_process( COMMAND git symbolic-ref HEAD  WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}
                 COMMAND tail -c -7 COMMAND tr -d "\\n" OUTPUT_VARIABLE GitBranch )
set( GitCommand2 "log --pretty=format:\\\"%h\\\" -1" )
execute_process( COMMAND git log --pretty=format:%h WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}
                 COMMAND head -n 1  COMMAND tr -d "\\n" OUTPUT_VARIABLE GitHashtag )
set( LANA_VERSION \"${GitBranch}-${GitHashtag}\" )

set( LIBLEMON_ROOT "$ENV{HOME}/lemon" CACHE PATH "Additional search directory for lemon library" )


if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release CACHE STRING "Choose the type of build, options are: Debug Release Profile" FORCE)
endif(NOT CMAKE_BUILD_TYPE)

set(CMAKE_CXX_FLAGS_DEBUG "-Wall -Wextra -Wno-long-long -Wno-unused-parameter -Wno-unknown-pragmas -g -ggdb -O0 -DDEBUG -DVERBOSE")
set(CMAKE_CXX_FLAGS_RELEASE "-Wall -Wextra -Wno-long-long -Wno-unused-parameter -Wno-unknown-pragmas -O3 -DNDEBUG")
set(CMAKE_CXX_FLAGS_PROFILE "-Wall -Wextra -Wno-long-long -Wno-unused-parameter -Wno-unknown-pragmas -O3 -g -ggdb -DNDEBUG")

set( CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wno-write-strings " )
add_definitions( -DIL_STD )
configure_file( src/config.h.in ${PROJECT_SOURCE_DIR}/src/config.h )

if( ${CMAKE_SYSTEM_NAME} MATCHES "Darwin" )
  set( CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -static-libstdc++ " )
endif()

include_directories( ${CPLEX_INC_DIR} ${CONCERT_INC_DIR} "${LIBLEMON_ROOT}/include" "src" "src/externlib" )
link_directories( ${CPLEX_LIB_DIR} ${CONCERT_LIB_DIR} "${LIBLEMON_ROOT}/lib" )

set( CommonLibs
	emon
)

add_executable( lana ${Lana_Hdr} ${Lana_Src} ${Tinyxml_Hdr} ${Tinyxml_Src} ${Gml-parser_Hdr} ${Gml-parser_Src} ${CommonObjs} )
target_link_libraries( lana ${CommonLibs} )


