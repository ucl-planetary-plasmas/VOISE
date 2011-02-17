
if HAVE_MCC

mccVOISE: webVOISE VOISE verifyImage

MCC = mcc
MCCFLAGS = -m -v
MCCFLAGSNODISPLAY = -R '-nodisplay' -m -v
MCCFLAGSNOJVM = -R '-nojvm,-nodisplay' -m -v

matlabdir = $(top_srcdir)/matlab

WEBVOISE_SRC = \
$(matlabdir)/webVOISE.m \
$(matlabdir)/startup.m \
$(matlabdir)/addSeedsToVR.m \
$(matlabdir)/addSeedToVD.m \
$(matlabdir)/mergeVD.m \
$(matlabdir)/benchVD.m \
$(matlabdir)/parseArgs.m \
$(matlabdir)/computeHCThreshold.m \
$(matlabdir)/plotHistHC.m \
$(matlabdir)/computeVDFast.m \
$(matlabdir)/plotVOISE.m \
$(matlabdir)/computeVD.m \
$(matlabdir)/prctile.m \
$(matlabdir)/divideVD.m \
$(matlabdir)/printFigure.m \
$(matlabdir)/printSeeds.m \
$(matlabdir)/getCentroidSeed.m \
$(matlabdir)/printSet.m \
$(matlabdir)/getCentroidVD.m \
$(matlabdir)/printVD.m \
$(matlabdir)/getDefaultVOISEParams.m \
$(matlabdir)/readVOISEconf.m \
$(matlabdir)/printVOISEsetup.m \
$(matlabdir)/getEquidistantPoint.m \
$(matlabdir)/randomSeeds.m \
$(matlabdir)/getVDOp.m \
$(matlabdir)/removeSeedFromVD.m \
$(matlabdir)/getVRclosure.m \
$(matlabdir)/seedDist.m \
$(matlabdir)/getVRvertices.m \
$(matlabdir)/homogeneousCriteria.m \
$(matlabdir)/updateNeighboursListsFromVoronoin.m \
$(matlabdir)/initVD.m \
$(matlabdir)/VOISE.m \
$(matlabdir)/isXinW.m \
$(matlabdir)/loadImage.m \
$(matlabdir)/fitswrite.m

VERIFYIMAGE_SRC = \
$(matlabdir)/verifyImage.m \
$(matlabdir)/loadImage.m \
$(matlabdir)/startup.m

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
	@echo "*********** testing VOISE demo1 ***********"
	./run_VOISE.sh $(MATLABROOT) ../share/webVOISEdemo1.dat

demo2: VOISE
	@echo "*********** testing VOISE demo2 ***********"
	./run_VOISE.sh $(MATLABROOT) ../share/webVOISEdemo2.dat

check: demo1 demo2
checkweb: webdemo1 webdemo2

endif

clean-local:
	-rm -f *VOISE* *verifyImage* *startup* *.log readme.txt

