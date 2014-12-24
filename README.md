# mc_cpp : A C++ Monte Carlo simulations software.
[![Build Status](https://travis-ci.org/FHedin/mc_cpp.svg?branch=experimental)](https://travis-ci.org/FHedin/mc_cpp)
----------------------------------------------
## TESTED PLATFORMS AND COMPILERS:
----------------------------------------------
This software **requires a C++ compiler conforming to the C++2011 standard** (see http://en.wikipedia.org/wiki/C%2B%2B11) 

Successfully tested platforms and compilers:
* Linux x86_64 :
  * GCC-G++ 4.8 or newer
  * LLVM Clang++ 3.4 or newer 

* MS Windows x86_64 : 
  * MS Visual C++ 2013 (with Visual Studio 2013 Ultimate Edition)
  * MinGW GCC-G++ 4.8

----------------------------------------------
## COMPILE & INSTALL:
----------------------------------------------
* Be sure to have CMAKE installed (http://www.cmake.org/), available on most repositories.

* Create a build directory and move to that directory: 
        mkdir build && cd build

* For buildind a debug or release version : 
        cmake -DCMAKE_BUILD_TYPE=Debug ..
or
        cmake -DCMAKE_BUILD_TYPE=Release ..

* Debug builds are slower but useful when debugging with gdb or valgrind.

* Then, once cmake built a Makefile, just execute :
        make

* For a verbose make, use : 
        make VERBOSE=1

* Please never edit this autogenerated Makefile, edit the CMakeLists.txt instead.

----------------------------------------------
## NOTE:
----------------------------------------------
All files excepted the ones under ./rapidxml-1.13 (MIT license) are licensed under the GNU GPL3 license.

----------------------------------------------
## LICENSING:
----------------------------------------------
Copyright (C) 2013.2014  Florent Hedin

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>.

----------------------------------------------
## FEATURES:
----------------------------------------------
For the moment, limited to : 
* Simulation of any molecular system in the NVT ensemble, using cubic or orthorombic PBCs.
* The forcefield has to be in the MDBAS format (https://github.com/pcazade/MDBas.git) ; see test directory, converters may be provided in the future.

MC algorithms available are :
* Metropolis [1,2]

----------------------------------------------
## REFERENCES:
----------------------------------------------
[1] Metropolis, N.; Rosenbluth, A.W.; Rosenbluth, M.N.; Teller, A.H.; Teller, E. (1953). "Equations of State Calculations by Fast Computing Machines". Journal of Chemical Physics 21 (6): 1087–1092. Bibcode:1953JChPh..21.1087M. doi:10.1063/1.1699114.

[2] Hastings, W.K. (1970). "Monte Carlo Sampling Methods Using Markov Chains and Their Applications".

