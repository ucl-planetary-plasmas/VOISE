
        Voronoi Image SEgementation (VOISE) 

        Patrick Guio (2009)


Introduction
------------

 This is the MATLAB® implementation of the VOISE image segmentation
 described in the paper by Guio and Achilleos (2009).

Installation
------------

 In order to install VOISE from repository (CVS or github), you will
 need the GNU autotools (autoconf, automake and libtool) to generate the
 configuration command configure.
 In the root directory of VOISE, run the command

 % autoreconf -vif

 Then to configure VOISE run the configure command. Some information is
 available in the INSTALL file or get some help using the command

 % configure -h

 Running

 % configure

 should detect the location of MATLAB® but you can indicate its location
 setting the environment variable MATLABROOT.
 
 If you have installed the MATLAB® wrapper for NAIF SPICE available on
 github you can set the environment variable MICE_PATH to the matlab
 directory in order to integrate the functionalities of SPICE within VOISE.

 To start VOISE, go to the matlab directory and run matlab

 % matlab

 You should get the following output

 Setting up VOISE -- version 2.0.3
 Mice path:  /home/patrick/research/codes/spice/matlab/
 Version  : Mice 1.6.0 05-NOV-2021 (EDW) (JDR) (MCS) compiled Jan  2 2022 02:12:27

 or similar.

 You can generate some faster C++ functionalities using the MATLAB® mex
 compiler. In order to do so, you need a functioning mex compiler for C++
 and the boost C++ library https://www.boost.org/

 From matlab, run the command

 >> mexCompile

 If succesfull, you can run the compiled version of some of the
 Voronoi-based algorithms

 Some test functions are available (testVOISE.m and testVOISE1.m)
 

License
-------

 This program is free software; you can redistribute it and/or modify it
 under the terms of the GNU General Public License as published by 
 the Free Software Foundation; either version 2.  of the License, or 
 (at your option) any later version.

 This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.




Content
-------

 In the matlab directory, you will find three subdirectories automatically
 detected and added to the matlab path

 * visUtils: contains 2D/3D visualisation utilities (some requiring SPICE).

 * limbFit : contains limb fitting functionality with examples.

 * imageUtils: contains image processing tools dedicated to VOISE.


Bibliography
------------


 Guio, P. and Achilleos, N. The VOISE Algorithm: a Versatile Tool for Automatic Segmentation of Astronomical Images. Mon. Not. R. Astron. Soc., 398(3):1254–1262, September 2009.

