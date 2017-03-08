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
                           src/srfft.c

main_LDADD = ../../kaldi/src/decoder/kaldi-decoder.a \
              ../../kaldi/src/lat/kaldi-lat.a \
              ../../kaldi/src/fstext/kaldi-fstext.a \
              ../../kaldi/src/hmm/kaldi-hmm.a \
              ../../kaldi/src/feat/kaldi-feat.a \
              ../../kaldi/src/transform/kaldi-transform.a \
              ../../kaldi/src/gmm/kaldi-gmm.a \
              ../../kaldi/src/tree/kaldi-tree.a \
              ../../kaldi/src/util/kaldi-util.a \
              ../../kaldi/src/thread/kaldi-thread.a \
              ../../kaldi/src/matrix/kaldi-matrix.a \
              ../../kaldi/src/base/kaldi-base.a \
              /usr/lib/libatlas.so.3 \
              /usr/lib/libf77blas.so.3 \
              /usr/lib/libcblas.so.3 \
              /usr/lib/liblapack_atlas.so.3 \
              ../../kaldi/tools/openfst/lib/libfst.a \
              libpocketkaldi.a \
              -ldl

TESTS = hashlist_test list_test fst_test srfft_test pcm_reader_test
check_PROGRAMS = hashlist_test list_test fst_test srfft_test pcm_reader_test

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

pcm_reader_test_SOURCES = test/pcm_reader_test.cc
pcm_reader_test_CXXFLAGS = $(AM_CXXFLAGS) -std=c++11 -I$(top_srcdir)/src -DTESTDIR=\"$(top_srcdir)/test/\"
pcm_reader_test_LDADD = libpocketkaldi.a

