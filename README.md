lana
=======

Local pairwise network alignment

Dependencies
------------

* LEMON 1.3

Compiling
---------

Get lana from github:

    git clone <HTTPS clone URL (see on the right side of this page)>


First, LEMON 1.3 needs to be installed:

    wget http://lemon.cs.elte.hu/pub/sources/lemon-1.3.tar.gz
    tar xvzf lemon-1.3.tar.gz
    cd lemon-1.3
    cmake -DCMAKE_INSTALL_PREFIX=~/lemon
    make install
    
Note: On Mac OS 10.9, comment out the following two lines and add the code below at line 162 in `CMakeLists.txt` before `make install`


    #ADD_SUBDIRECTORY(demo) 
    #ADD_SUBDIRECTORY(tools)
    
    if( ${CMAKE_SYSTEM_NAME} MATCHES "Darwin" )
      set( CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -stdlib=libstdc++ " )
    endif()

You can remove the LEMON sources now, i.e., `rm -rf lemon-1.3`. Next, natalie can be compiled:

    mkdir build
    cd build
    cmake ..
    make

In case auto-detection of LEMON or CPLEX fails, do

    cmake -DLIBLEMON_ROOT=~/lemon ..

Running lana
---------------

To run lana:

    ./lana -if1 0 -if2 0 -ifm 0 -g1 ../data/rno.gml -g2 ../data/hsa.gml -gm ../data/rno_hsa.seqSim

For usage instructions specify `-h`.
