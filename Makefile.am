ACLOCAL_AMFLAGS = -I m4
AUTOMAKE_OPTIONS = gnu subdir-objects foreign
AM_CFLAGS = -std=c99
AM_CXXFLAGS = -I$(top_srcdir)/../kaldi/src -I$(top_srcdir)/../kaldi/tools/openfst/src/include/ -I$(top_srcdir)/../kaldi/tools/ATLAS/include/ -I$(top_srcdir)/src -g -fPIC -std=c++11
AM_LDFLAGS = -pthread -L/usr/local/cuda-8.0/lib64/
LIBS = -ldl -lm -lblas

bin_PROGRAMS = pocketkaldi
pocketkaldi_SOURCES = src/main.cc

lib_LIBRARIES = libpocketkaldi.a
libpocketkaldi_a_SOURCES = src/util.cc \
                           src/fst.cc \
                           src/matrix.cc \
                           src/pcm_reader.cc \
                           src/decoder.cc \
                           src/srfft.cc \
                           src/fbank.cc \
                           src/strlcpy.cc \
                           src/cmvn.cc \
                           src/transition.cc \
                           src/nnet.cc \
                           src/am.cc \
                           src/vector.cc \
                           src/decodable.cc \
                           src/symbol_table.cc \
                           src/pocketkaldi.cc \
                           src/introspect.cc  \
                           src/hashtable.cc \
                           src/configuration.cc

pocketkaldi_LDADD = libpocketkaldi.a -lstdc++

if ENABLE_TOOLS
    bin_PROGRAMS += compute_fbank extract_id2pdf
    compute_fbank_SOURCES = tool/compute_fbank.cc
    compute_fbank_LDADD = $(KALDI_ROOT)/src/util/kaldi-util.a \
                          $(KALDI_ROOT)/src/thread/kaldi-thread.a \
                          $(KALDI_ROOT)/src/feat/kaldi-feat.a \
                          $(KALDI_ROOT)/src/matrix/kaldi-matrix.a \
                          $(KALDI_ROOT)/src/base/kaldi-base.a \
                          $(ATLASLIBS) \
                          libpocketkaldi.a

    extract_id2pdf_SOURCES = tool/extract_id2pdf.cc
    extract_id2pdf_LDADD = $(KALDI_ROOT)/src/util/kaldi-util.a \
                           $(KALDI_ROOT)/src/hmm/kaldi-hmm.a \
                           $(KALDI_ROOT)/src/matrix/kaldi-matrix.a \
                           $(KALDI_ROOT)/src/base/kaldi-base.a \
                           $(ATLASLIBS) \
                           libpocketkaldi.a
endif

TESTS = list_test \
        fst_test \
        srfft_test \
        fbank_test \
        cmvn_test \
        nnet_test \
        symbol_table_test \
        hashtable_test \
        configuration_test \
        pool_test

check_PROGRAMS = list_test \
                 fst_test \
                 srfft_test \
                 fbank_test \
                 cmvn_test \
                 nnet_test \
                 symbol_table_test \
                 hashtable_test \
                 configuration_test \
                 pool_test

configuration_test_SOURCES = test/configuration_test.cc
configuration_test_CXXFLAGS = $(AM_CXXFLAGS) -std=c++11 -I$(top_srcdir)/src -DTESTDIR=\"$(top_srcdir)/test/\"
configuration_test_LDADD = libpocketkaldi.a

hashtable_test_SOURCES = test/hashtable_test.cc
hashtable_test_CXXFLAGS = $(AM_CXXFLAGS) -std=c++11 -I$(top_srcdir)/src
hashtable_test_LDADD = libpocketkaldi.a

list_test_SOURCES = test/list_test.cc
list_test_CXXFLAGS = $(AM_CXXFLAGS) -std=c++11 -I$(top_srcdir)/src
list_test_LDADD = libpocketkaldi.a

fst_test_SOURCES = test/fst_test.cc
fst_test_CXXFLAGS = $(AM_CXXFLAGS) -std=c++11 -I$(top_srcdir)/src -DTESTDIR=\"$(top_srcdir)/test/\"
fst_test_LDADD = libpocketkaldi.a

symbol_table_test_SOURCES = test/symbol_table_test.cc
symbol_table_test_CXXFLAGS = $(AM_CXXFLAGS) -std=c++11 -I$(top_srcdir)/src -DTESTDIR=\"$(top_srcdir)/test/\"
symbol_table_test_LDADD = libpocketkaldi.a

srfft_test_SOURCES = test/srfft_test.cc
srfft_test_CXXFLAGS = $(AM_CXXFLAGS) -std=c++11 -I$(top_srcdir)/src
srfft_test_LDADD = libpocketkaldi.a

fbank_test_SOURCES = test/fbank_test.cc
fbank_test_CXXFLAGS = $(AM_CXXFLAGS) -std=c++11 -I$(top_srcdir)/src -DTESTDIR=\"$(top_srcdir)/test/\"
fbank_test_LDADD = libpocketkaldi.a

cmvn_test_SOURCES = test/cmvn_test.cc
cmvn_test_CXXFLAGS = $(AM_CXXFLAGS) -std=c++11 -I$(top_srcdir)/src -DTESTDIR=\"$(top_srcdir)/test/\"
cmvn_test_LDADD = libpocketkaldi.a

nnet_test_SOURCES = test/nnet_test.cc
nnet_test_CXXFLAGS = $(AM_CXXFLAGS) -std=c++11 -I$(top_srcdir)/src -DTESTDIR=\"$(top_srcdir)/test/\"
nnet_test_LDADD = libpocketkaldi.a

pool_test_SOURCES = test/pool_test.cc
pool_test_CXXFLAGS = $(AM_CXXFLAGS) -std=c++11 -I$(top_srcdir)/src
pool_test_LDADD = libpocketkaldi.a

if ENABLE_TOOLS
    TESTS_ENVIRONMENT = export testdir=$(top_srcdir)/test && export kaldiroot=$(KALDI_ROOT) &&
    TESTS += test/test_compute_fbank.sh
endif