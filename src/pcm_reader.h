// Changed from feat/wave-reader.cc

// Copyright 2009-2011  Karel Vesely;  Petr Motlicek
//                2013  Florent Masson
//                2013  Johns Hopkins University (author: Daniel Povey)

// See ../../COPYING for clarification regarding multiple authors
//
// Licensed under the Apache License, Version 2.0 (the "License");
// you may not use this file except in compliance with the License.
// You may obtain a copy of the License at
//
//  http://www.apache.org/licenses/LICENSE-2.0
//
// THIS CODE IS PROVIDED *AS IS* BASIS, WITHOUT WARRANTIES OR CONDITIONS OF ANY
// KIND, EITHER EXPRESS OR IMPLIED, INCLUDING WITHOUT LIMITATION ANY IMPLIED
// WARRANTIES OR CONDITIONS OF TITLE, FITNESS FOR A PARTICULAR PURPOSE,
// MERCHANTABLITY OR NON-INFRINGEMENT.
// See the Apache 2 License for the specific language governing permissions and
// limitations under the License.

// 2017-1-3

#ifndef POCKETKALDI_PCM_READER_H_
#define POCKETKALDI_PCM_READER_H_

#include <stdint.h>
#include "alloc.h"
#include "util.h"

typedef struct pk_16kpcm_t {
  float *data;
  int num_samples;
} pk_16kpcm_t;

// Reads 16k sampling rate, mono-channel, PCM formatted wave file, and stores
// the data into data. If any error occured, set status to failed
POCKETKALDI_EXPORT
void pk_16kpcm_read(
    const char *filename,
    pk_16kpcm_t *pcm_data,
    pk_alloc_t *alloc,
    pk_status_t *status);

#endif

