// Created at 2017-03-29

#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "pocketkaldi.h"
#include "util.h"

void check_and_report_error(const pk_status_t *status) {
  if (!status->ok) {
    printf("pocketkaldi: %s\n", status->message);
    exit(1);
  }
}

// Process one utterance
void process_audio(pk_t *recognizer, const char *filename) {
  pk_status_t status;
  pk_status_init(&status);

  pk_utterance_t utt;
  pk_utterance_init(&utt);
  pk_read_audio(&utt, filename, &status);
  check_and_report_error(&status);

  pk_process(recognizer, &utt);
  printf("%s\t%s\t%f\n", filename, utt.hyp, utt.loglikelihood_per_frame);

  pk_utterance_destroy(&utt);
}

// Process a list of utterances
void process_scp(pk_t *recognizer, const char *filename) {
  pk_status_t status;
  pk_status_init(&status);
  pk_readable_t *fd = pk_readable_open(filename, &status);
  check_and_report_error(&status);

  // Read each line in scp file
  char audio_file[2048];
  while (pk_readable_readline(fd, audio_file, sizeof(audio_file), &status)) {
    process_audio(recognizer, audio_file);
  }
  check_and_report_error(&status);
}

// Print the usage of this program and exit
void print_usage() {
  puts("Usage: pocketkaldi <model-file> <input-file>");
  puts("  Input-file:");
  puts("    *.wav: decode this file.");
  puts("    *.scp: decode audios listed in it.");
  exit(1);
}

int main(int argc, char **argv) {
  if (argc != 3) print_usage();

  const char *model_file = argv[1];
  const char *input_file = argv[2];
  if (strlen(input_file) < 4) print_usage();

  pk_t recognizer;
  pk_status_t status;
  pk_status_init(&status);
  pk_init(&recognizer);
  pk_load(&recognizer, model_file, &status);
  check_and_report_error(&status);

  const char *suffix = input_file + strlen(input_file) - 4;
  if (strcmp(suffix, ".wav") == 0) {
    process_audio(&recognizer, input_file);
  } else {
    process_scp(&recognizer, input_file);
  }

  pk_destroy(&recognizer);
  return 0;
}
