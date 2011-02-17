

SUBDIRS = matlab mcc share

EXTRA_DIST = mcc/Makefile 

update-cl:
	update-cl.sh

if HAVE_MCC

ROOT = mcc@PACKAGE_NAME@-@VERSION@-@target@

mccVOISE:
	@(cd mcc && $(MAKE) $(AM_MAKEFLAGS) $@) || exit 1;

webdemo1 webdemo2 demo1 demo2 check checkweb: mccVOISE
	@(cd mcc && $(MAKE) $(AM_MAKEFLAGS) $@) || exit 1;

dist-mccVOISE: mccVOISE
	@rm -rf $(ROOT)
	@mkdir -p $(ROOT) $(ROOT)/mcc $(ROOT)/share $(ROOT)/share/input $(ROOT)/share/output
	@cp -pr AUTHORS ChangeLog NEWS README $(ROOT)
	@cp -pr mcc/Makefile mcc/*VOISE* mcc/*verifyImage* mcc/*.log mcc/readme.txt $(ROOT)/mcc
	@cp -p share/VOISEtiming.mat share/webVOISEdemo*.dat $(ROOT)/share
	@cp -p share/input/north_proj.mat share/input/sampleint.fits $(ROOT)/share/input
	@tar -cvf $(ROOT).tar $(ROOT)
	@bzip2 -9 $(ROOT).tar
	@rm -rf $(ROOT)

endif
