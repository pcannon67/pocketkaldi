ACLOCAL_AMFLAGS = -I m4
AUTOMAKE_OPTIONS = gnu subdir-objects foreign
AM_CXXFLAGS = -I$(top_srcdir)/../kaldi/src -I$(top_srcdir)/../kaldi/tools/openfst/src/include/ -I$(top_srcdir)/../kaldi/tools/ATLAS/include/ -I$(top_srcdir)/src -g -fPIC 
AM_LDFLAGS = -pthread -ldl

bin_PROGRAMS = main
main_SOURCES = src/main.cc 

lib_LIBRARIES = libpocketkaldi.a
libpocketkaldi_a_SOURCES = src/hashlist.c \
                           src/util.c \
                           src/fst.c \
                           src/matrix.c \
                           src/pcm_reader.c \
                           src/decoder.cc \
                           src/srfft.c \
                           src/fbank.c

main_LDADD = $(KALDI_ROOT)/src/decoder/kaldi-decoder.a \
              $(KALDI_ROOT)/src/lat/kaldi-lat.a \
              $(KALDI_ROOT)/src/fstext/kaldi-fstext.a \
              $(KALDI_ROOT)/src/hmm/kaldi-hmm.a \
              $(KALDI_ROOT)/src/feat/kaldi-feat.a \
              $(KALDI_ROOT)/src/transform/kaldi-transform.a \
              $(KALDI_ROOT)/src/gmm/kaldi-gmm.a \
              $(KALDI_ROOT)/src/tree/kaldi-tree.a \
              $(KALDI_ROOT)/src/util/kaldi-util.a \
              $(KALDI_ROOT)/src/thread/kaldi-thread.a \
              $(KALDI_ROOT)/src/matrix/kaldi-matrix.a \
              $(KALDI_ROOT)/src/base/kaldi-base.a \
              $(ATLASLIBS) \
              $(KALDI_ROOT)/tools/openfst/lib/libfst.a \
              libpocketkaldi.a \
              -ldl

if ENABLE_TOOLS
    bin_PROGRAMS += compute_fbank
    compute_fbank_SOURCES = tool/compute_fbank.cc
    compute_fbank_LDADD = $(KALDI_ROOT)/src/util/kaldi-util.a \
                          $(KALDI_ROOT)/src/thread/kaldi-thread.a \
                          $(KALDI_ROOT)/src/feat/kaldi-feat.a \
                          $(KALDI_ROOT)/src/matrix/kaldi-matrix.a \
                          $(KALDI_ROOT)/src/base/kaldi-base.a \
                          $(ATLASLIBS) \
                          libpocketkaldi.a
endif

TESTS = hashlist_test list_test fst_test srfft_test fbank_test
check_PROGRAMS = hashlist_test list_test fst_test srfft_test fbank_test

hashlist_test_SOURCES = test/hashlist_test.cc
hashlist_test_CXXFLAGS = $(AM_CXXFLAGS) -std=c++11 -I$(top_srcdir)/src
hashlist_test_LDADD = libpocketkaldi.a

list_test_SOURCES = test/list_test.cc
list_test_CXXFLAGS = $(AM_CXXFLAGS) -std=c++11 -I$(top_srcdir)/src
list_test_LDADD = libpocketkaldi.a

fst_test_SOURCES = test/fst_test.cc
fst_test_CXXFLAGS = $(AM_CXXFLAGS) -std=c++11 -I$(top_srcdir)/src -DTESTDIR=\"$(top_srcdir)/test/\"
fst_test_LDADD = libpocketkaldi.a

srfft_test_SOURCES = test/srfft_test.cc
srfft_test_CXXFLAGS = $(AM_CXXFLAGS) -std=c++11 -I$(top_srcdir)/src
srfft_test_LDADD = libpocketkaldi.a

fbank_test_SOURCES = test/fbank_test.cc
fbank_test_CXXFLAGS = $(AM_CXXFLAGS) -std=c++11 -I$(top_srcdir)/src -DTESTDIR=\"$(top_srcdir)/test/\"
fbank_test_LDADD = libpocketkaldi.a

if ENABLE_TOOLS
    TESTS_ENVIRONMENT = export testdir=$(top_srcdir)/test && export kaldiroot=$(KALDI_ROOT) &&
    TESTS += test/test_compute_fbank.sh
endif