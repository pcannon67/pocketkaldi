// 2017-03-23 vector.h
// 2017-05-25 matrix/kaldi-vector.h

// Copyright 2009-2012   Ondrej Glembek;  Microsoft Corporation;  Lukas Burget;
//                       Saarland University (Author: Arnab Ghoshal);
//                       Ariya Rastrow;  Petr Schwarz;  Yanmin Qian;
//                       Karel Vesely;  Go Vivace Inc.;  Arnab Ghoshal
//                       Wei Shi;
//                2015   Guoguo Chen

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


#ifndef POCKETKALDI_VECTOR_H_
#define POCKETKALDI_VECTOR_H_

#include "matrix.h"
#include "status.h"
#include "util.h"

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

namespace pocketkaldi {

template<typename Real>
class SubVector;

///  Provides a vector abstraction class.
///  This class provides a way to work with vectors in kaldi.
///  It encapsulates basic operations and memory optimizations.
template<typename Real>
class VectorBase {
 public:
  /// Set vector to all zeros.
  void SetZero();

  /// Returns true if matrix is all zeros.
  bool IsZero(Real cutoff = 1.0e-06) const;     // replace magic number

  /// Set all members of a vector to a specified value.
  void Set(Real f);

  /// Set vector to random normally-distributed noise.
  void SetRandn();

  /// Returns the  dimension of the vector.
  inline int Dim() const { return dim_; }

  /// Returns a pointer to the start of the vector's data.
  inline Real* Data() { return data_; }

  /// Returns a pointer to the start of the vector's data (const).
  inline const Real* Data() const { return data_; }

  /// Indexing  operator (const).
  inline Real operator() (int i) const {
    assert(i < dim_);
    return *(data_ + i);
  }

  /// Indexing operator (non-const).
  inline Real & operator() (int i) {
    assert(i < dim_);
    return *(data_ + i);
  }

  /** @brief Returns a sub-vector of a vector (a range of elements).
   *  @param o [in] Origin, 0 < o < Dim()
   *  @param l [in] Length 0 < l < Dim()-o
   *  @return A SubVector object that aliases the data of the Vector object.
   *  See @c SubVector class for details   */
  SubVector<Real> Range(const int o, const int l) {
    return SubVector<Real>(*this, o, l);
  }

  /** @brief Returns a const sub-vector of a vector (a range of elements).
   *  @param o [in] Origin, 0 < o < Dim()
   *  @param l [in] Length 0 < l < Dim()-o
   *  @return A SubVector object that aliases the data of the Vector object.
   *  See @c SubVector class for details   */
  const SubVector<Real> Range(const int o, const int l) const {
    return SubVector<Real>(*this, o, l);
  }

  /// Copy data from another vector (must match own size).
  void CopyFromVec(const VectorBase<Real> &v);

  /// Copy data from another vector of different type (double vs. float)
  template<typename OtherReal>
  void CopyFromVec(const VectorBase<OtherReal> &v);

  // Returns dot product between this and v.
  Real VecVec(const VectorBase<Real> &v) const;

  /// Apply natural log to all elements.  Throw if any element of
  /// the vector is negative (but doesn't complain about zero; the
  /// log will be -infinity
  void ApplyLog();

  /// Applies floor to all elements. Returns number of elements floored.
  int ApplyFloor(Real floor_val);

  /// Add vector : *this = *this + alpha * rv (with casting between floats and
  /// doubles)
  template<typename OtherReal>
  void AddVec(const Real alpha, const VectorBase<OtherReal> &v);

  /// Multiplies all elements by this constant.
  void Scale(Real alpha);

  friend class VectorBase<double>;
  friend class VectorBase<float>;
 protected:
  /// Destructor;  does not deallocate memory, this is handled by child classes.
  /// This destructor is protected so this object so this object can only be
  /// deleted via a child.
  ~VectorBase() {}

  /// Empty initializer, corresponds to vector of zero size.
  explicit VectorBase(): data_(NULL), dim_(0) {
  }

// Took this out since it is not currently used, and it is possible to create
// objects where the allocated memory is not the same size as dim_ : Arnab
//  /// Initializer from a pointer and a size; keeps the pointer internally
//  /// (ownership or non-ownership depends on the child class).
//  explicit VectorBase(Real* data, MatrixIndexT dim)
//      : data_(data), dim_(dim) {}

  // Arnab : made this protected since it is unsafe too.
  /// Load data into the vector: sz must match own size.
  void CopyFromPtr(const Real* Data, int sz);

  /// data memory area
  Real* data_;
  /// dimension of vector
  int dim_;
}; // class VectorBase

/** @brief A class representing a vector.
 *
 *  This class provides a way to work with vectors in kaldi.
 *  It encapsulates basic operations and memory optimizations.  */
template<typename Real>
class Vector: public VectorBase<Real> {
 public:
  enum {
    kSetZero,
    kUndefined,
    kCopyData
  };

  /// Constructor that takes no arguments.  Initializes to empty.
  Vector(): VectorBase<Real>() {}

  /// Constructor with specific size.  Sets to all-zero by default
  /// if set_zero == false, memory contents are undefined.
  explicit Vector(const int s, int resize_type = kSetZero)
      : VectorBase<Real>() {  Resize(s, resize_type);  }

  // Read vector from ReadableFile
  Status Read(util::ReadableFile *fd);

  /// Swaps the contents of *this and *other.  Shallow swap.
  void Swap(Vector<Real> *other);

  /// Set vector to a specified size (can be zero).
  /// The value of the new data depends on resize_type:
  ///   -if kSetZero, the new data will be zero
  ///   -if kUndefined, the new data will be undefined
  ///   -if kCopyData, the new data will be the same as the old data in any
  ///      shared positions, and zero elsewhere.
  /// This function takes time proportional to the number of data elements.
  void Resize(int length, int resize_type = kSetZero);

  /// Destructor.  Deallocates memory.
  ~Vector() { Destroy(); }

 private:
  /// Init assumes the current contents of the class are invalid (i.e. junk or
  /// has already been freed), and it sets the vector to newly allocated memory
  /// with the specified dimension.  dim == 0 is acceptable.  The memory contents
  /// pointed to by data_ will be undefined.
  void Init(const int dim);

  /// Destroy function, called internally.
  void Destroy();

};


/// Represents a non-allocating general vector which can be defined
/// as a sub-vector of higher-level vector [or as the row of a matrix].
template<typename Real>
class SubVector : public VectorBase<Real> {
 public:
  /// Constructor from a Vector or SubVector.
  /// SubVectors are not const-safe and it's very hard to make them
  /// so for now we just give up.  This function contains const_cast.
  SubVector(const VectorBase<Real> &t, int origin, int length) {
    // following assert equiv to origin>=0 && length>=0 &&
    // origin+length <= rt.dim_
    assert(origin + length <= t.Dim());
    VectorBase<Real>::data_ = const_cast<Real*> (t.Data()+origin);
    VectorBase<Real>::dim_ = length;
  }

  /// Constructor from a pointer to memory and a length.  Keeps a pointer
  /// to the data but does not take ownership (will never delete).
  SubVector(Real *data, int length) : VectorBase<Real> () {
    VectorBase<Real>::data_ = data;
    VectorBase<Real>::dim_ = length;
  }

  ~SubVector() {}  ///< Destructor (does nothing; no pointers are owned here).

 private:
  /// Disallow assignment operator.
  SubVector & operator = (const SubVector &other) {}
};

}  // namespace pocketkaldi

#endif  // POCKETKALDI_VECTOR_H_
