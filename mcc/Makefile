#
# $Id: Makefile,v 1.7 2010/10/07 15:54:59 patrick Exp $
#
# Copyright (c) 2010
# Patrick Guio <p.guio@ucl.ac.uk>
#
# All Rights Reserved.
#
# This program is free software; you can redistribute it and/or modify it
# under the terms of the GNU General Public License as published by the
# Free Software Foundation; either version 2.  of the License, or (at your
# option) any later version.
#
# This program is distributed in the hope that it will be useful, but
# WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General
# Public License for more details.
#

# You can invoke this Makefile by just typing
# make
# or 
# make MATLABROOT=[directory where MATLAB is installed]
# default MATLABROOT is /usr/local/matlabR2010b
# You can also invole the demo/check targets
# make demo1
# make demo2
# make check
# The target check makes demo1 and demo2

default: webVOISE VOISE verifyImage

# Patrick's laptop
ifeq ($(HOST),jazz2)
MATLABROOT = /usr/local/matlabR2010b
endif

# Dugan's laptop
ifeq ($(OSTYPE),darwin10.0)
MATLABROOT = /Applications/MATLAB/MATLAB_Compiler_Runtime/v714
endif

MCC = mcc
MCCFLAGS = -m -v
MCCFLAGSNODISPLAY = -R '-nodisplay' -m -v
MCCFLAGSNOJVM = -R '-nojvm,-nodisplay' -m -v

matlab_dir = ../matlab

WEBVOISE_SRC = \
$(matlab_dir)/webVOISE.m \
$(matlab_dir)/startup.m \
$(matlab_dir)/addSeedsToVR.m \
$(matlab_dir)/addSeedToVD.m \
$(matlab_dir)/mergeVD.m \
$(matlab_dir)/benchVD.m \
$(matlab_dir)/parseArgs.m \
$(matlab_dir)/computeHCThreshold.m \
$(matlab_dir)/plotHistHC.m \
$(matlab_dir)/computeVDFast.m \
$(matlab_dir)/plotVOISE.m \
$(matlab_dir)/computeVD.m \
$(matlab_dir)/prctile.m \
$(matlab_dir)/divideVD.m \
$(matlab_dir)/printFigure.m \
$(matlab_dir)/printSeeds.m \
$(matlab_dir)/getCentroidSeed.m \
$(matlab_dir)/printSet.m \
$(matlab_dir)/getCentroidVD.m \
$(matlab_dir)/printVD.m \
$(matlab_dir)/getDefaultVOISEParams.m \
$(matlab_dir)/readVOISEconf.m \
$(matlab_dir)/printVOISEsetup.m \
$(matlab_dir)/getEquidistantPoint.m \
$(matlab_dir)/randomSeeds.m \
$(matlab_dir)/getVDOp.m \
$(matlab_dir)/removeSeedFromVD.m \
$(matlab_dir)/getVRclosure.m \
$(matlab_dir)/seedDist.m \
$(matlab_dir)/getVRvertices.m \
$(matlab_dir)/homogeneousCriteria.m \
$(matlab_dir)/updateNeighboursListsFromVoronoin.m \
$(matlab_dir)/initVD.m \
$(matlab_dir)/VOISE.m \
$(matlab_dir)/isXinW.m \
$(matlab_dir)/loadImage.m \
$(matlab_dir)/fitswrite.m

VERIFYIMAGE_SRC = \
$(matlab_dir)/verifyImage.m \
$(matlab_dir)/loadImage.m \
$(matlab_dir)/startup.m

VOISE: $(WEBVOISE_SRC)
	$(MCC) $(MCCFLAGS) $^ -o $@
	-rm -f *.o*

webVOISE: $(WEBVOISE_SRC)
	$(MCC) $(MCCFLAGSNODISPLAY) $^ -o $@
	-rm -f *.o*

verifyImage: $(VERIFYIMAGE_SRC)
	$(MCC) $(MCCFLAGSNOJVM) $^ -o $@
	-rm -f *.o*

webdemo1: webVOISE
	@echo "*********** testing webVOISE demo1 ***********"
	./run_webVOISE.sh $(MATLABROOT) ../share/webVOISEdemo1.dat

webdemo2: webVOISE
	@echo "*********** testing webVOISE demo2 ***********"
	./run_webVOISE.sh $(MATLABROOT) ../share/webVOISEdemo2.dat

demo1: VOISE
	@echo "*********** testing webVOISE demo1 ***********"
	./run_VOISE.sh $(MATLABROOT) ../share/webVOISEdemo1.dat

demo2: VOISE
	@echo "*********** testing webVOISE demo1 ***********"
	./run_VOISE.sh $(MATLABROOT) ../share/webVOISEdemo1.dat

check: demo1 demo2
checkweb: webdemo1 webdemo2

clean:
	-rm -f *webVOISE* *verifyImage* mccExcludedFiles.log readme.txt
