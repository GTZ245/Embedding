AUTOMAKE_OPTIONS = foreign

INCLUDES = -I$(includedir) -I$(OFFLINE_MAIN)/include -I$(ROOTSYS)/include

mydatadir = $(datadir)/$(PACKAGE_INCLUDE_DIR)
dist_mydata_DATA = \
  field_map.root

lib_LTLIBRARIES = \
  libRun14AuAuLeptonEvent.la \
  libRun14AuAuLeptonConvReco.la \
  libsvxcentana.la

noinst_HEADERS = \
  MyEventLinkDef.h \
  MyEvent.h \
  MyEventConstants.h \
  ReconstructionLinkDef.h 

include_HEADERS = \
  Run14AuAuLeptonCombyConstants.h \
  embedana.h \
  Reconstruction.h

libRun14AuAuLeptonEvent_la_SOURCES = \
  MyEvent.C \
  MyEvent_Dict.C

libRun14AuAuLeptonConvReco_la_SOURCES = \
  Reconstruction.cpp \
  Reconstruction_Dict.C

libsvxcentana_la_SOURCES = \
  embedana.cc \
  QACuts.C \
  embedana_Dict.C

libRun14AuAuLeptonEvent_la_LIBADD = \
$(ROOT_LIBS)  

libsvxcentana_la_LIBADD = \
  -L$(libdir) \
  -L$(OFFLINE_MAIN)/lib \
  -lTOAD \
  -lphool \
  -lfun4all \
  -luti \
  -lCNT \
  -lsvx \
  -lreactionplane \
  -lsimreco \
  -lcompactCNT \
  -lrecal

# Rule for generating CINT dictionaries from class headers.
MyEvent_Dict.C: \
  MyEvent.h \
  MyEventConstants.h \
  MyEventLinkDef.h 
	  rootcint -f $@ -c $(DEFAULT_INCLUDES) $(INCLUDES) $^

Reconstruction_Dict.C: \
  Reconstruction.h \
  ReconstructionLinkDef.h
	rootcint -f $@ -c $(DEFAULT_INCLUDES) $(AM_CPPFLAGS) $^


embedana_Dict.C: \
  embedana.h \
  Run14AuAuLeptonCombyConstants.h \
  embedanaLinkDef.h 
	  rootcint -f $@ -c $(DEFAULT_INCLUDES) $(INCLUDES) $^


clean-local:
	rm -f *Dict*

