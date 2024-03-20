/*
   This file is provided as part of the STAG library and released under the GPL
   license.
*/

#include "kde.h"

/**
 * Compute the squared distance between the given two data points.
 *
 * @param u the first data point
 * @param v the second data point
 * @return the squared distance between \f$u\f$ and \f$v\f$.
 */
StagReal squared_distance(const stag::DataPoint& u, const stag::DataPoint& v) {
  assert(u.dimension == v.dimension);
  StagReal result = 0;
  for (StagUInt i = 0; i < u.dimension; i++) {
    result += SQR(u.coordinates[i] - v.coordinates[i]);
  }
  return result;
}

StagReal stag::gaussian_kernel(StagReal a, StagReal c) {
  return exp(-(a * c));
}

StagReal stag::gaussian_kernel(StagReal a, const stag::DataPoint& u,
                         const stag::DataPoint& v) {
  return stag::gaussian_kernel(a, squared_distance(u, v));
}
