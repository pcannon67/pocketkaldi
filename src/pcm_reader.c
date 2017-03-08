// 2017-01-23

#include "pcm_reader.h"

#include <assert.h>
#include <stdbool.h>
#include <stdio.h>
#include <string.h>
#include <stdint.h>
#include "alloc.h"
#include "util.h"

// Reads an int32 value from ptr and move it forward
static int32_t read_int32(char **ptr) {
  int32_t val = *((int32_t *)*ptr);
  *ptr += sizeof(int32_t);
  return val;
}

// Reads an int16 value from ptr and move it forward
static int16_t read_int16(char **ptr) {
  int16_t val = *((int16_t *)*ptr);
  *ptr += sizeof(int16_t);
  return val;
}

// Reads an int8 value from ptr and move it forward
static int8_t read_int8(char **ptr) {
  int8_t val = *((int8_t *)*ptr);
  *ptr += sizeof(int8_t);
  return val;
}

// Reads a string with length strlen(expected) from ptr, and compared with
// expected. If the same returns true, otherwise returns false. Then move ptr
// forward strlen(expected) bytes
static bool check_tag(char **ptr, const char *expected) {
  int len = strlen(expected);
  bool is_same = strncmp(*ptr, expected, len) == 0;
  *ptr += len;
  return is_same;
}

// Reads 16k sampling rate, mono-channel, PCM formatted wave file, and stores
// the data into data. If any error occured, set status to failed
void pk_16kpcm_read(
    const char *filename,
    pk_16kpcm_t *pcm_data,
    pk_alloc_t *alloc,
    pk_status_t *status) {
  FILE *fd = fopen(filename, "rb");
  if (fd == NULL) {
    pk_status_fail(status, PK_STATUS_IOERROR, "failed to open: %s", filename);
  }

  // Gets file length
  int64_t file_size = 0;
  if (status->ok) {
    fseek(fd, 0, SEEK_END);
    file_size = ftell(fd);
    fseek(fd, 0, SEEK_SET);

    // A PCM file should have at least 44bytes
    if (file_size < 44) {
      pk_status_fail(status, PK_STATUS_IOERROR, "file is too small: %s", filename);
    }
  }

  // Read file content into pcm_buffer
  char *pcm_buffer = NULL;
  char *current_ptr = NULL;
  if (status->ok) {
    pcm_buffer = (char *)pk_alloc(alloc, file_size);
    current_ptr = pcm_buffer;
    int bytes_read = fread(pcm_buffer, 1, file_size, fd);
    if (bytes_read != file_size) {
      pk_status_fail(
          status,
          PK_STATUS_IOERROR,
          "failed to read, %lld bytes expected, but %lld got: %s",
          file_size,
          bytes_read,
          filename);
    }
  }

  // RIFF chunk
  if (status->ok && check_tag(&current_ptr, "RIFF") == false) {
    pk_status_fail(status, PK_STATUS_CORRUPTED,"chunk_name == 'RIFF' expected: %s", filename);
  }

  // Chunk size
  if (status->ok) {
    int chunk_size = read_int32(&current_ptr);
    if (chunk_size != file_size - 8) {
      pk_status_fail(
          status,
          PK_STATUS_CORRUPTED,
          "chunk_size == %lld expected, but %lld found: %s",
          file_size - 8,
          chunk_size,
          filename);
    }
  }

  // Format == "WAVE"
  if (status->ok && check_tag(&current_ptr, "WAVE") == false) {
    pk_status_fail(
        status,
        PK_STATUS_CORRUPTED,
        "Format == 'WAVE' expected: %s",
        filename);
  }

  // subchunk1 is "fmt "
  if (status->ok && check_tag(&current_ptr, "fmt ") == false) {
    pk_status_fail(
        status,
        PK_STATUS_CORRUPTED,
        "subchunk1 == 'fmt ' expected: %s",
        filename);
  }  

  // subchunk1_size
  if (status->ok) {
    int subchunk1_size = read_int32(&current_ptr);
    if (subchunk1_size != 16) {
      pk_status_fail(
          status,
          PK_STATUS_CORRUPTED,
          "subchunk1_size == 16 expected, but %d found: %s",
          subchunk1_size,
          filename);
    }
  }

  // audio_format
  if (status->ok) {
    int audio_format = read_int16(&current_ptr);
    if (audio_format != 1) {
      pk_status_fail(
          status,
          PK_STATUS_CORRUPTED,
          "audio_format == 1 (PCM) expected, but %d found: %s",
          audio_format,
          filename);
    }
  }

  // num_channels
  if (status->ok) {
    int num_channels = read_int16(&current_ptr);
    if (num_channels != 1) {
      pk_status_fail(
          status,
          PK_STATUS_CORRUPTED,
          "num_channels == 1 (mono) expected, but %d found: %s",
          num_channels,
          filename);
    }
  }

  // sample_rate
  int sample_rate = 0;
  if (status->ok) {
    sample_rate = read_int32(&current_ptr);
    if (sample_rate != 16000) {
      pk_status_fail(
          status,
          PK_STATUS_CORRUPTED,
          "sample_rate == 1 (mono) expected, but %d found: %s",
          sample_rate,
          filename);
    }
  }

  // bytes_rate, block_align, bits_per_sample
  int bits_per_sample = 0;
  if (status->ok) {
    int bytes_rate = read_int32(&current_ptr);
    int block_align = read_int16(&current_ptr);
    bits_per_sample = read_int16(&current_ptr);

    if (bytes_rate != sample_rate * bits_per_sample / 8) {
      pk_status_fail(
          status,
          PK_STATUS_CORRUPTED,
          "bytes_rate == %d expected, but %d found: %s",
          sample_rate * bits_per_sample / 8,
          bytes_rate,
          filename);
    }

    if (block_align != bits_per_sample / 8) {
      pk_status_fail(
          status,
          PK_STATUS_CORRUPTED,
          "block_align == %d expected, but %d found: %s",
          bits_per_sample / 8,
          block_align,
          filename); 
    }
  }

  // subchunk2 "data"
  if (status->ok && check_tag(&current_ptr, "data") == false) {
    pk_status_fail(
        status,
        PK_STATUS_CORRUPTED,
        "subchunk2 == 'data' expected: %s",
        filename);
  }

  // subchunk2_size
  int64_t subchunk2_size = 0;
  if (status->ok) {
    subchunk2_size = read_int32(&current_ptr);
    if (subchunk2_size != file_size - 44) {
      pk_status_fail(
          status,
          PK_STATUS_CORRUPTED,
          "subchunk2_size == %lld expected, but %lld found: %s",
          file_size - 44,
          subchunk2_size,
          filename);
    }
  }

  // Read data
  if (status->ok) {
    int num_samples = (int)subchunk2_size / (bits_per_sample / 8);
    pcm_data->data = (float *)pk_alloc(alloc, num_samples * sizeof(float));
    pcm_data->num_samples = num_samples;
    for (int i = 0; i < num_samples && status->ok; ++i) {
      switch (bits_per_sample) {
        case 8:
          pcm_data->data[i] = read_int8(&current_ptr);
          break;
        case 16:
          pcm_data->data[i] = read_int16(&current_ptr);
          break;
        case 32:
          pcm_data->data[i] = read_int32(&current_ptr);
          break;
        default:
          pk_status_fail(
              status,
              PK_STATUS_CORRUPTED,
              "bits_per_sample == 8, 16 or 32 expected, but %d found: %s",
              bits_per_sample,
              filename);
      }
    }
  }

  if (status->ok) {
    assert(current_ptr == pcm_buffer + file_size);
  }

  // Free resources
  if (fd != NULL) fclose(fd);
  if (pcm_buffer != NULL) pk_free(alloc, pcm_buffer);
}
