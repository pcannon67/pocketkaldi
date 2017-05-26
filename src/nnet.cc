// Created at 2017-03-13

#include "nnet.h"

#include <assert.h>
#include <math.h>

namespace pocketkaldi {

LinearLayer::LinearLayer(
    const MatrixBase<float> &W,
    const VectorBase<float> &b) {
  assert(b.Dim() == W.NumRows() && 
         "linear layer: dimension mismatch in W and b");
  W_.Resize(W.NumCols(), W.NumRows());
  W_.CopyFromMat(W, MatrixBase<float>::kTrans);
  b_.Resize(b.Dim());
  b_.CopyFromVec(b);
}

void LinearLayer::Propagate(
    const MatrixBase<float> &in,
    Matrix<float> *out) const {
  out->Resize(in.NumRows(), W_.NumCols());

  // xW
  SimpleMatMat(in, W_, out);

  // + b
  for (int row_idx = 0; row_idx < out->NumRows(); ++row_idx) {
    SubVector<float> row = out->Row(row_idx);
    row.AddVec(1.0f, b_);
  }
}

void SoftmaxLayer::Propagate(
    const MatrixBase<float> &in,
    Matrix<float> *out) const {
  out->Resize(in.NumRows(), in.NumCols());
  out->CopyFromMat(in);
  for (int row_idx = 0; row_idx < out->NumRows(); ++row_idx) {
    SubVector<float> row = out->Row(row_idx);
    row.ApplySoftMax();
  }
}

void ReLULayer::Propagate(
    const MatrixBase<float> &in,
    Matrix<float> *out) const {
  out->Resize(in.NumRows(), in.NumCols());
  out->CopyFromMat(in);
  for (int row_idx = 0; row_idx < out->NumRows(); ++row_idx) {
    SubVector<float> row = out->Row(row_idx);
    for (int col_idx = 0; col_idx < row.Dim(); ++ col_idx) {
      if (row(col_idx) < 0.0f) row(col_idx) = 0.0f;
    }
  }
}

void NormalizeLayer::Propagate(
    const MatrixBase<float> &in,
    Matrix<float> *out) const {
  float D = in.NumCols();
  out->Resize(in.NumRows(), in.NumCols());
  out->CopyFromMat(in);
  for (int row_idx = 0; row_idx < out->NumRows(); ++row_idx) {
    SubVector<float> row = out->Row(row_idx);

    double squared_sum = row.VecVec(row);
    float scale = static_cast<float>(sqrt(D / squared_sum));
    row.Scale(scale);
  }
}

Nnet::Nnet() {
}

Status Nnet::ReadLayer(util::ReadableFile *fd) {
  Matrix<float> W;
  Vector<float> b;

  // Read section name
  int32_t section_size;
  PK_CHECK_STATUS(fd->ReadAndVerifyString(PK_NNET_LAYER_SECTION));
  PK_CHECK_STATUS(fd->ReadValue<int32_t>(&section_size));

  // Read layer type
  int layer_type;
  PK_CHECK_STATUS(fd->ReadValue<int32_t>(&layer_type));

  // Check size of this section
  int expected_size = 4;
  if (expected_size != section_size) {
    return Status::Corruption(util::Format(
        "read_layer: section_size == {} expected, but {} found ({})",
        expected_size,
        section_size,
        fd->filename()));
  }

  // Read additional parameters and initialize layer
  float scale = 0.0f;
  LinearLayer *layer;
  switch (layer_type) {
  case Layer::kLinear:
    PK_CHECK_STATUS(W.Read(fd));
    PK_CHECK_STATUS(b.Read(fd));

    layers_.emplace_back(new LinearLayer(W, b));
    break;
  case Layer::kReLU:
    layers_.emplace_back(new ReLULayer());
    break;
  case Layer::kNormalize:
    layers_.emplace_back(new NormalizeLayer());
    break;
  case Layer::kSoftmax:
    layers_.emplace_back(new SoftmaxLayer());
    break;
  default:
    return Status::Corruption(util::Format(
        "read_layer: unexpected layer type: {} ({})",
        layer_type,
        fd->filename()));
  }

  return Status::OK();
}

Status Nnet::Read(util::ReadableFile *fd) {
  // Read section name
  int32_t section_size;
  PK_CHECK_STATUS(fd->ReadAndVerifyString(PK_NNET_SECTION));
  PK_CHECK_STATUS(fd->ReadValue<int32_t>(&section_size));

  int num_layers;
  PK_CHECK_STATUS(fd->ReadValue<int32_t>(&num_layers));

  // Read each layers
  for (int layer_idx = 0; layer_idx < num_layers; ++layer_idx) {
    PK_CHECK_STATUS(ReadLayer(fd));
  }

  return Status::OK();
}

void Nnet::Propagate(const pk_matrix_t *in, pk_matrix_t *out) const {
  SubMatrix<float> input(in->data, in->ncol, in->nrow, in->nrow);
  Matrix<float> layer_input, layer_output;
  
  layer_input.Resize(input.NumRows(), input.NumCols());
  layer_input.CopyFromMat(input);
  for (const std::unique_ptr<Layer> &layer : layers_) {
    layer->Propagate(layer_input, &layer_output);
    layer_input.Swap(&layer_output);
  }

  pk_matrix_resize(out, layer_input.NumCols(), layer_input.NumRows());
  SubMatrix<float> output(out->data, out->ncol, out->nrow, out->nrow);
  output.CopyFromMat(layer_input);
}

}  // namespace pocketkaldi
