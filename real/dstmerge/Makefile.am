
AUTOMAKE_OPTIONS = foreign

INCLUDES = -I. -I$(srcdir) -I$(includedir) -I$(OFFLINE_MAIN)/include -I$(ROOTSYS)/include

lib_LTLIBRARIES = libDstMerger.la

noinst_HEADERS = AnalysisLinkDef.h 


include_HEADERS = DstMerger.h SkipEvents.h


libDstMerger_la_SOURCES = DstMerger.cc SkipEvents.cc AnalysisDict.C


libDstMerger_la_LIBADD = -L$(OFFLINE_MAIN)/lib -L$(libdir) -lphool -lfun4all -luti

AnalysisDict.C: \
  DstMerger.h \
  SkipEvents.h \
  AnalysisLinkDef.h
	rootcint -f $@ -c $(INCLUDES) $^

# # Rule for generating CINT dictionaries from class headers.
# %Dict.C %Dict.h: %.hh
# 	rootcint -f $*Dict.C -c $(INCLUDES) $*.hh

