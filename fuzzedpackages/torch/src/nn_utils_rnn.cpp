#include "torch_types.h"
#include "utils.h"

// [[Rcpp::export]]
Rcpp::XPtr<XPtrTorchPackedSequence> cpp_nn_utils_rnn_pack_padded_sequence (
    Rcpp::XPtr<XPtrTorchTensor> input, Rcpp::XPtr<XPtrTorchTensor> lengths,
    bool batch_first, bool enforce_sorted)
{
  XPtrTorchPackedSequence out = lantern_nn_utils_rnn_pack_padded_sequence(input->get(), lengths->get(), batch_first,
                                            enforce_sorted);
  return make_xptr<XPtrTorchPackedSequence>(out);
}

// [[Rcpp::export]]
Rcpp::XPtr<XPtrTorchPackedSequence> cpp_nn_utils_pack_sequence (
    Rcpp::XPtr<XPtrTorchTensorList> sequence, bool enforce_sorted) {
  XPtrTorchPackedSequence out = lantern_nn_utils_rnn_pack_sequence(
    sequence->get(), enforce_sorted
  );
  return make_xptr<XPtrTorchPackedSequence>(out);
}

// [[Rcpp::export]]
Rcpp::XPtr<XPtrTorchTensorList> cpp_nn_utils_pad_packed_sequence (
  Rcpp::XPtr<XPtrTorchPackedSequence> sequence,
  bool batch_first,
  double padding_value,
  Rcpp::XPtr<XPtrTorchoptional_int64_t> total_length
) {
  XPtrTorchTensorList out = lantern_nn_utils_rnn_pad_packed_sequence(
    sequence->get(),
    batch_first, padding_value,
    total_length->get()
  );  
  return make_xptr<XPtrTorchTensorList>(out);
}

// [[Rcpp::export]]
Rcpp::XPtr<XPtrTorchTensor> cpp_nn_utils_pad_sequence (
  Rcpp::XPtr<XPtrTorchTensorList> sequence,
  bool batch_first,
  double padding_value
) {
  XPtrTorchTensor out = lantern_nn_utils_rnn_pad_sequence(sequence->get(), batch_first,
                                                           padding_value);
  return make_xptr<XPtrTorchTensor>(out);
}

// [[Rcpp::export]]
Rcpp::XPtr<XPtrTorchPackedSequence> cpp_nn_utils_PackedSequence_new(Rcpp::XPtr<XPtrTorchTensor> data,
                                                                    Rcpp::XPtr<XPtrTorchTensor> batch_sizes,
                                                                    Rcpp::XPtr<XPtrTorchTensor> sorted_indices,
                                                                    Rcpp::XPtr<XPtrTorchTensor> unsorted_indices)
{
  XPtrTorchPackedSequence out = lantern_nn_utils_rnn_PackedSequence_new( 
    data->get(),
    batch_sizes->get(),
    sorted_indices->get(),
    unsorted_indices->get()
  );
  return make_xptr<XPtrTorchPackedSequence>(out);
}

// [[Rcpp::export]]
Rcpp::XPtr<XPtrTorchTensor> cpp_nn_utils_PackedSequence_data (Rcpp::XPtr<XPtrTorchPackedSequence> x)
{
  XPtrTorchTensor out = lantern_nn_utils_PackedSequence_data(x->get());
  return make_xptr<XPtrTorchTensor>(out);
}

// [[Rcpp::export]]
Rcpp::XPtr<XPtrTorchTensor> cpp_nn_utils_PackedSequence_batch_sizes (Rcpp::XPtr<XPtrTorchPackedSequence> x)
{
  XPtrTorchTensor out = lantern_nn_utils_PackedSequence_batch_sizes(x->get());
  return make_xptr<XPtrTorchTensor>(out);
}

// [[Rcpp::export]]
Rcpp::XPtr<XPtrTorchTensor> cpp_nn_utils_PackedSequence_sorted_indices (Rcpp::XPtr<XPtrTorchPackedSequence> x)
{
  XPtrTorchTensor out = lantern_nn_utils_PackedSequence_sorted_indices(x->get());
  return make_xptr<XPtrTorchTensor>(out);
}

// [[Rcpp::export]]
Rcpp::XPtr<XPtrTorchTensor> cpp_nn_utils_PackedSequence_unsorted_indices (Rcpp::XPtr<XPtrTorchPackedSequence> x)
{
  XPtrTorchTensor out = lantern_nn_utils_PackedSequence_unsorted_indices(x->get());
  return make_xptr<XPtrTorchTensor>(out);
}
