/*
   Copyright (c) 2004-2005 Massachusetts Institute of Technology.
   Authors: Alexandr Andoni (andoni@mit.edu), Piotr Indyk (indyk@mit.edu)

   Modified by Peter Macgregor, 2024.

   This file is provided as part of the STAG library and released under the GPL
   license.
*/

#include <algorithm>

#include <definitions.h>
#include <lsh.h>
#include <random.h>

#define TWO_ROOT_TWOPI 5.0132565
#define TWO_ROOT_TWO 2.828427124

// Returns TRUE iff |p1-p2|_2^2 <= threshold
inline bool isDistanceSqrLeq(StagUInt dimension, const stag::DataPoint& p1,
                             const stag::DataPoint& p2, StagReal threshold){
  StagReal result = 0;

  for (StagUInt i = 0; i < dimension; i++){
    StagReal temp = p1.coordinates[i] - p2.coordinates[i];
    result += SQR(temp);
    if (result > threshold){
      return false;
    }
  }

  return true;
}

// Generate a random real distributed uniformly in [rangeStart,
// rangeEnd]. Input must satisfy: rangeStart <= rangeEnd. The
// granularity of generated random reals is given by RAND_MAX.
StagReal genUniformRandom(StagReal rangeStart, StagReal rangeEnd){
  assert(rangeStart <= rangeEnd);

  std::uniform_real_distribution<StagReal> dist(rangeStart, rangeEnd);
  StagReal r = dist(*stag::get_global_rng());

  assert(r >= rangeStart && r <= rangeEnd);
  return r;
}

// Generate a random real from normal distribution N(0,1).
StagReal genGaussianRandom(){
  std::normal_distribution<StagReal> dist(0, 1);
  StagReal z = dist(*stag::get_global_rng());
  return z;
}

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
// Implementation of the LSHFunction class.
//------------------------------------------------------------------------------
stag::LSHFunction::LSHFunction(StagUInt dimension) {
  dim = dimension;

  // The variable is a random Gaussian vector.
  a.reserve(dim);
  for(StagUInt d = 0; d < dim; d++){
    a.emplace_back(genGaussianRandom());
  }

  // The variable b is a random offset on the random vector.
  b = genUniformRandom(0, LSH_PARAMETER_W);
}

StagUInt stag::LSHFunction::apply(const DataPoint& point) {
  assert(point.dimension == dim);

  StagReal value = 0;
  for(StagUInt d = 0; d < dim; d++){
    value += point.coordinates[d] * a[d];
  }

  return (StagUInt) floor((value + b) / LSH_PARAMETER_W);
}

StagReal stag::LSHFunction::collision_probability(StagReal c) {
  StagReal eight_over_c_squared = 8 / SQR(c);
  return - ((1 / TWO_ROOT_TWOPI)
            * c
            * exp(-eight_over_c_squared)
            * (exp(eight_over_c_squared) - 1))
         + erf(TWO_ROOT_TWO / c);
}

//------------------------------------------------------------------------------
// Implementation of the E2LSH class.
//------------------------------------------------------------------------------
stag::E2LSH::E2LSH(RNNParametersT algParameters, StagUInt nPoints,
                   std::vector<DataPoint>& dataSet){
  initialise_fields_from_parameters(algParameters, nPoints);

  // Set the fields <nPoints> and <points>.
  points = dataSet;

  // Given the number of points, let's set the hash table size to be 1/100 the
  // number of points to be hashed.
  auto hashTableSize = (StagUInt) std::max((StagUInt) 1, nPoints / 100);

  // Initialise the empty hash tables
  for (StagUInt l = 0; l < parameterL; l++) {
    hashTables.emplace_back(hashTableSize, parameterK);
  }

  // Add the points to the hash tables
  for(StagUInt i = 0; i < nPoints; i++){
    for(StagUInt l = 0; l < parameterL; l++){
      BucketHashingIndexT bucket_index = hashTables[0].compute_bucket_index(
          compute_lsh(l, dataSet[i]));
      hashTables[l].add_bucket_entry(bucket_index, i);
    }
  }
}

void stag::E2LSH::initialise_hash_functions() {
  lshFunctions.reserve(parameterL);
  for(StagUInt i = 0; i < parameterL; i++){
    lshFunctions.emplace_back();
    lshFunctions[i].reserve(dimension);
    for(StagUInt j = 0; j < parameterK; j++){
      lshFunctions[i].emplace_back(dimension);
    }
  }
}

void stag::E2LSH::initialise_fields_from_parameters(RNNParametersT algParameters, StagUInt nPointsEstimate) {
  parameterR2 = algParameters.parameterR2;
  checkDistanceWhenReturning = algParameters.checkDistance;
  parameterK = algParameters.parameterK;
  parameterL = algParameters.parameterL;
  dimension = algParameters.dimension;

  // create the hash functions
  initialise_hash_functions();

  // init fields that are used only in operations ("temporary" variables for operations).

  // init the vector and the vector
  // <precomputedHashesOfULSHs>
  sizeMarkedPoints = nPointsEstimate;
  markedPoints.resize(nPointsEstimate, false);
  markedPointsIndices.resize(sizeMarkedPoints);
}

std::vector<StagUInt> stag::E2LSH::compute_lsh(StagUInt gNumber, const DataPoint& point) {
  std::vector<StagUInt> result(parameterK);
  for(StagUInt i = 0; i < parameterK; i++){
    result[i] = lshFunctions[gNumber][i].apply(point);
  }
  return result;
}

std::vector<stag::DataPoint> stag::E2LSH::get_near_neighbors(const DataPoint& query) {
  std::vector<DataPoint> near_points;

  StagUInt nMarkedPoints = 0;// the number of marked points
  for(StagUInt l = 0; l < parameterL; l++){
    BucketHashingIndexT bucket_index = hashTables[0].compute_bucket_index(
        compute_lsh(l, query));
    LSHBucket* bucket = hashTables[l].get_bucket(bucket_index);

    // circle through the bucket and add to <result> the points that are near.
    if (bucket != nullptr) {
      for (StagInt candidatePIndex : bucket->points) {
        DataPoint candidatePoint = points[candidatePIndex];
        if (!checkDistanceWhenReturning ||
            isDistanceSqrLeq(dimension, query, candidatePoint, parameterR2)){
          if (!markedPoints[candidatePIndex]) {
            near_points.push_back(candidatePoint);
            markedPointsIndices[nMarkedPoints] = candidatePIndex;
            markedPoints[candidatePIndex] = true; // do not include more points with the same index
            nMarkedPoints++;
          }
        }
      }
    }
  }

  // we need to clear the array nnStruct->nearPoints for the next query.
  for(StagUInt i = 0; i < nMarkedPoints; i++){
    assert(markedPoints[markedPointsIndices[i]]);
    markedPoints[markedPointsIndices[i]] = false;
  }

  return near_points;
}
