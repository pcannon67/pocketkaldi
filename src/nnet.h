// Created at 2017-03-13

#ifndef POCKETKALDI_NNET_H_
#define POCKETKALDI_NNET_H_

#include "matrix.h"
#include "util.h"

#define PK_NNET_SECTION "NNT0"
#define PK_NNET_LAYER_SECTION "LAY0"

// Layer types
#define PK_NNET_LINEAR_LAYER 0
#define PK_NNET_RELU_LAYER 1
#define PK_NNET_NORMALIZE_LAYER 2
#define PK_NNET_SOFTMAX_LAYER 3
#define PK_NNET_ADD_LAYER 4
#define PK_NNET_MUL_LAYER 5


namespace pocketkaldi {

// The base class for different type of layers
class Layer {
 public:
  // Kinds of linear types
  enum {
    kLinear = 0,
    kReLU = 1,
    kNormalize = 2,
    kSoftmax = 3
  };

  // Propogate a batch of input vectors through this layer. And the batch of
  // output vectors are in `out`
  virtual void Propagate(const pk_matrix_t *in, pk_matrix_t *out) const = 0;
  virtual ~Layer() {}
};

// Linear layer: x^T dot W + b
class LinearLayer : public Layer {
 public:
  // Initialize the linear layer with parameter W and b. It just copies the
  // values from W and b.
  LinearLayer(const pk_matrix_t *W, const pk_vector_t *b);
  ~LinearLayer();

  void Propagate(const pk_matrix_t *in, pk_matrix_t *out) const override;

 private:
  pk_matrix_t W_;
  pk_vector_t b_;
};

// Softmax layer
class SoftmaxLayer : public Layer {
 public:
  void Propagate(const pk_matrix_t *in, pk_matrix_t *out) const override;
};

// ReLU layer
class ReLULayer : public Layer {
 public:
  void Propagate(const pk_matrix_t *in, pk_matrix_t *out) const override;
};

// Normalize layer
class NormalizeLayer : public Layer {
 public:
  void Propagate(const pk_matrix_t *in, pk_matrix_t *out) const override;
};

// The neural network class. It have a stack of different kinds of `Layer`
// instances. And the batch matrix could be propogate through this neural
// network using `Propagate` method
class Nnet {
 public:
  Nnet();

  // Read the nnet from file
  Status Read(pk_readable_t *fd);

  // Propogate batch matrix through this neural network
  void Propagate(const pk_matrix_t *in, pk_matrix_t *out) const;

 private:
  std::vector<std::unique_ptr<Layer>> layers_;

  // Read a layer from `fd` and store into layers_
  Status ReadLayer(pk_readable_t *fd);
};

}  // namespace pocketkaldi

#endif
