AUTOMAKE_OPTIONS = foreign

BUILT_SOURCES = \
  testexternals.C

INCLUDES = -I$(includedir) -I$(OFFLINE_MAIN)/include -I$(ROOTSYS)/include

lib_LTLIBRARIES = \
  libembedreco.la 

noinst_HEADERS = \
  BbcEmbedreco.h \
  DchEmbedreco.h \
  PadEmbedreco.h \
  TecEmbedreco.h \
  TofEmbedreco.h \
  CrkEmbedreco.h \
  EmcEmbedreco.h \
  AccEmbedreco.h \
  TofwEmbedreco.h \
  MixEmbedreco.h\
  EvaEmbedreco.h \
  EmbedVertexSelect.h\
  SvxEmbedSimhit.h \
  SvxEmbedEval.h \
  EmbedrecoLinkDef.h 

#do not install the header for reco modules because no classes are allowed to inherent from them. 
include_HEADERS = 

libembedreco_la_SOURCES = \
  BbcEmbedreco.C \
  DchEmbedreco.C \
  PadEmbedreco.C \
  TecEmbedreco.C \
  TofEmbedreco.C \
  CrkEmbedreco.C \
  EmcEmbedreco.C \
  AccEmbedreco.C \
  TofwEmbedreco.C \
  MixEmbedreco.C\
  EvaEmbedreco.C \
  EmbedVertexSelect.C\
  SvxEmbedSimhit.C \
  SvxEmbedEval.C \
  Embedreco_Dict.C

AM_LDFLAGS = -L$(libdir) -L$(OFFLINE_MAIN)/lib

libembedreco_la_LIBADD = \
  -lSubsysReco\
  -lcteval \
  -ldch \
  -lemc \
  -lembed \
  -lphool \
  -ltof   \
  -ltofw  \
  -lsvx

noinst_PROGRAMS = \
  testexternals

testexternals_LDADD = \
  libembedreco.la

testexternals.C:
	echo "//*** this is a generated file. Do not commit, do not edit" > $@
	echo "int main()" >> $@
	echo "{" >> $@
	echo "  return 0;" >> $@
	echo "}" >> $@

# Rule for generating CINT dictionaries from class headers.
Embedreco_Dict.C: \
  BbcEmbedreco.h \
  DchEmbedreco.h \
  PadEmbedreco.h \
  TecEmbedreco.h \
  TofEmbedreco.h \
  TofwEmbedreco.h \
  CrkEmbedreco.h \
  EmcEmbedreco.h \
  AccEmbedreco.h \
  MixEmbedreco.h\
  EvaEmbedreco.h \
  EmbedVertexSelect.h\
  SvxEmbedSimhit.h \
  SvxEmbedEval.h \
  EmbedrecoLinkDef.h
	rootcint -f $@ -c $(DEFAULT_INCLUDES) $(INCLUDES) $^

clean-local:
	rm -f *Dict*

