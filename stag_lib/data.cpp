/*
   This file is provided as part of the STAG library and released under the GPL
   license.
*/
#include "data.h"

//------------------------------------------------------------------------------
// Implementation of the DataPoint class.
//------------------------------------------------------------------------------
stag::DataPoint::DataPoint(DenseMat& all_data, StagInt row_index) {
  dimension = all_data.cols();
  coordinates = all_data.row(row_index).data();
}

stag::DataPoint::DataPoint(std::vector<StagReal>& point_vector) {
  dimension = point_vector.size();
  coordinates = &point_vector[0];
}

//------------------------------------------------------------------------------
// Converting between data formats.
//------------------------------------------------------------------------------
std::vector<stag::DataPoint> stag::matrix_to_datapoints(DenseMat* data) {
  std::vector<stag::DataPoint> res;

  StagInt n = data->rows();
  StagInt d = data->cols();
  for (auto i = 0; i < n; i++) {
    res.emplace_back(d, data->row(i).data());
  }
  return res;
}
