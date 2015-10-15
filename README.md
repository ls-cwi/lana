lana
=======

Local pairwise network alignment

Authors
-------
* Jelmer Mulder
* Mohammed El-Kebir
* Gunnar W. Klau

Dependencies
------------

* LEMON 1.3
* Boost

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

You can remove the LEMON sources now, i.e., `rm -rf lemon-1.3`. Next, lana can be compiled:

    mkdir build
    cd build
    cmake ..
    make

* Boost

Compiling
------------




In case auto-detection of LEMON or CPLEX fails, do

    cmake -DLIBLEMON_ROOT=~/lemon ..

If desired, the flags `-DCMAKE_C_COMPILER` and `-DCMAKE_CXX_COMPILER` can be used with the `cmake` command to specify a different compiler.

Running lana
---------------

For complete run instructions, please run

	 ./lana -h

An example invocation:

    ./lana -if1 0 -if2 0 -ifm 0 -g1 ../data/mmu.gml -g2 ../data/cel.gml -gm ../data/mmu_cel.seqSim -o lana_mmu_cel.out
