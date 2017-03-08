// 2017-01-24

#include <stdio.h>
#include <assert.h>
#include "pcm_reader.h"
#include "util.h"

void TestPcmReader() {
  pk_status_t status;
  pk_16kpcm_t pcm_data;

  pk_status_init(&status);

  const char *pcm_filename = TESTDIR "en-us-hello.wav";
  pk_16kpcm_read(pcm_filename, &pcm_data, &status);
  if (!status.ok) {
    puts(status.message);
  }
  assert(status.ok);
  assert(pcm_data.num_samples == 7802);
}

int main() {
  TestPcmReader();
  return 0;
}