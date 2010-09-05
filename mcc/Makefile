#
# $Id: Makefile,v 1.1 2010/09/05 10:21:33 patrick Exp $
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

# You can invoke this Makefile by just 'make' or using
# make MATLABROOT=[directory where MATLAB is installed]
# if the MATLABROOT is different from /usr/local/matlab/R2010a

default: webVOISE

MATLABROOT = /usr/local/matlab/R2010a

MCC = $(MATLABROOT)/bin/mcc
MCCFLAGS = -R '-nodisplay' -m -v

matlab_dir = ../matlab

SRC = \
$(matlab_dir)/webVOISE.m \
$(matlab_dir)/startup.m \
$(matlab_dir)/start_VOISE.m \
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
$(matlab_dir)/exportfig.m \
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
$(matlab_dir)/loadImage.m

webVOISE: $(SRC)
	$(MCC) $(MCCFLAGS) $^ -o $@
	-rm -f *.o*

check: webVOISE
	@echo "*********** testing webVOISE ***********"
	run_webVOISE.sh $(MATLABROOT) ../share/webVOISEdemo2.dat

clean:
	-rm -f *webVOISE* mccExcludedFiles.log readme.txt
