ACLOCAL_AMFLAGS = -I m4
AUTOMAKE_OPTIONS = gnu subdir-objects foreign
AM_CXXFLAGS = -I$(top_srcdir)/../kaldi/src -I$(top_srcdir)/../kaldi/tools/openfst/src/include/ -I$(top_srcdir)/../kaldi/tools/ATLAS/include/ -I$(top_srcdir)/src -g -fPIC 
AM_LDFLAGS = -pthread -ldl

bin_PROGRAMS = main
main_SOURCES = src/main.cc \
               src/pk-simple-decoder.cc 

lib_LIBRARIES = libpocketkaldi.a
libpocketkaldi_a_SOURCES = src/alloc.c \
                           src/hashlist.c

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

TESTS = hashlist_test
check_PROGRAMS = hashlist_test
hashlist_test_SOURCES = src/alloc.c \
                        src/hashlist.c \
                        test/hashlist_test.cc
hashlist_test_CXXFLAGS = $(AM_CXXFLAGS) -std=c++11 -I$(top_srcdir)/src
