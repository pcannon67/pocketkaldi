// 2017-03-23


#ifndef POCKETKALDI_VECTOR_H_
#define POCKETKALDI_VECTOR_H_

#include "matrix.h"

typedef struct pk_matrix_t pk_matrix_t;

// pk_vector_t is a simple column vector
typedef struct pk_vector_t {
  int dim;
  float *data;
} pk_vector_t;

// Initialize the vector and fill with a number. If fill_with == NAN, doesn't
// fill anything
POCKETKALDI_EXPORT
void pk_vector_init(pk_vector_t *self, int dim, float fill_with);

// Read vector from file. If failed, set status to fail
POCKETKALDI_EXPORT
void pk_vector_read(pk_vector_t *self, pk_readable_t *fd, pk_status_t *status);

// Resize the vector
POCKETKALDI_EXPORT
void pk_vector_resize(pk_vector_t *self, int dim);

// Fill the vector with value
POCKETKALDI_EXPORT
void pk_vector_fill(pk_vector_t *self, float value);

// Copy n values from source to the vector. n should be less than self->dim
POCKETKALDI_EXPORT
void pk_vector_copyfrom(pk_vector_t *self, const float *source, int n);

// Copy vector from src to dest
POCKETKALDI_EXPORT
void pk_vector_copy(pk_vector_t *dest, const pk_vector_t *src);

// Apply floor to the vector
POCKETKALDI_EXPORT
void pk_vector_floor(pk_vector_t *self, float floor);

// Apply log to the vector
POCKETKALDI_EXPORT
void pk_vector_log(pk_vector_t *self);

// Apply scale to the vector
POCKETKALDI_EXPORT
void pk_vector_scale(pk_vector_t *self, float scale);

// Calculates the dot product of self and vec
POCKETKALDI_EXPORT
float pk_vector_dot(const pk_vector_t *self, const pk_vector_t *vec);

// Calculate (self^T dot W)^T, snd stores the output into out
POCKETKALDI_EXPORT
float pk_vector_dotmat(
    const pk_vector_t *self, 
    const pk_matrix_t *W,
    pk_vector_t *out);

// Add vec to self
POCKETKALDI_EXPORT
void pk_vector_add(const pk_vector_t *self, const pk_vector_t *vec);

// Add scale * vec to self
POCKETKALDI_EXPORT
void pk_vector_addscale(
    const pk_vector_t *self,
    float scale,
    const pk_vector_t *vec);

// Add scalar * vec to self
POCKETKALDI_EXPORT
void pk_vector_scalaradd(
    const pk_vector_t *self,
    float scalar,
    const pk_vector_t *vec);

// Subtract vec from self
POCKETKALDI_EXPORT
void pk_vector_subtract(const pk_vector_t *self, const pk_vector_t *vec);

// Add element-wise squaring of vec to self
POCKETKALDI_EXPORT
void pk_vector_add2(const pk_vector_t *self, const pk_vector_t *vec);

// Subtract element-wise squaring of vec from self
POCKETKALDI_EXPORT
void pk_vector_subtract2(const pk_vector_t *self, const pk_vector_t *vec);

// Creates a subvector according to self. It just borrows the memory of self.
// So to avoid issues it couldn't be destroyed and modified.
POCKETKALDI_EXPORT const pk_vector_t 
pk_vector_subvector(const pk_vector_t *self, int start, int dim);

// Destroy a vector
POCKETKALDI_EXPORT
void pk_vector_destroy(pk_vector_t *self);

#endif  // POCKETKALDI_VECTOR_H_
