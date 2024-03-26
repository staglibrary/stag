/*
   Copyright (c) 2004-2005 Massachusetts Institute of Technology.
   Authors: Alexandr Andoni (andoni@mit.edu), Piotr Indyk (indyk@mit.edu)

   Modified by Peter Macgregor, 2024.

   This file is provided as part of the STAG library and released under the GPL
   license.
*/

#include <algorithm>
#include <unordered_set>

#include "definitions.h"
#include "lsh.h"
#include "Graph/random.h"

#define TWO_ROOT_TWOPI 5.0132565
#define TWO_ROOT_TWO 2.828427124

// 4294967291 = 2^32-5
#define UH_PRIME_DEFAULT 4294967291U

// Returns TRUE iff |p1-p2|_2^2 <= threshold
// TODO: delete this method
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
// Implementation of the LSHFunction class.
//------------------------------------------------------------------------------
// TODO: combine multiple LSH functions (the K iterations) into a single matrix?
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

StagInt stag::LSHFunction::apply(const DataPoint& point) {
  assert(point.dimension == dim);

  StagReal value = 0;
  for(StagUInt d = 0; d < dim; d++){
    value += point.coordinates[d] * a[d];
  }

  return (StagInt) floor((value + b) / LSH_PARAMETER_W);
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
// Generate a random 32-bits unsigned (Uns32T) in the range
// [rangeStart, rangeEnd]. Inputs must satisfy: rangeStart <=
// rangeEnd.
StagInt genRandomInt(StagInt rangeStart, StagInt rangeEnd){
  assert(rangeStart <= rangeEnd);

  std::uniform_int_distribution<StagInt> dist(rangeStart, rangeEnd);
  StagInt r = dist(*stag::get_global_rng());

  assert(r >= rangeStart && r <= rangeEnd);
  return r;
}

stag::E2LSH::E2LSH(StagUInt K,
                   StagUInt L,
                   std::vector<DataPoint>& dataSet){
  if (dataSet.size() > 0) {
    dimension = dataSet[0].dimension;
  } else {
    dimension = 1;
  }

  parameterK = K;
  parameterL = L;
  StagUInt nPoints = dataSet.size();

  // create the hash functions
  initialise_hash_functions();

  // TODO: don't copy the dataset?
  points = dataSet;

  // Initialise the empty hash tables
  hashTables.resize(parameterL);

  // Add the points to the hash tables
  for(StagUInt i = 0; i < nPoints; i++){
    for(StagUInt l = 0; l < parameterL; l++){
      StagInt this_lsh = compute_lsh(l, dataSet[i]);
      if (!hashTables[l].contains(this_lsh)) {
        hashTables[l][this_lsh] = std::vector<StagUInt>();
      }
      hashTables[l][this_lsh].push_back(i);
    }
  }
}

void stag::E2LSH::initialise_hash_functions() {
  rnd_vec.resize(parameterK);
  for(StagUInt i = 0; i < parameterK; i++){
    rnd_vec[i] = genRandomInt(1, MAX_HASH_RND);
  }

  lshFunctions.reserve(parameterL);
  for(StagUInt i = 0; i < parameterL; i++){
    lshFunctions.emplace_back();
    lshFunctions[i].reserve(dimension);
    for(StagUInt j = 0; j < parameterK; j++){
      lshFunctions[i].emplace_back(dimension);
    }
  }
}

StagUInt stag::E2LSH::compute_lsh(StagUInt gNumber, const DataPoint& point) {
  StagInt h = 0;
  for(StagUInt i = 0; i < parameterK; i++){
    h += rnd_vec[i] * lshFunctions[gNumber][i].apply(point);

    if (h < 0) h += UH_PRIME_DEFAULT;
    if (h >= UH_PRIME_DEFAULT) h -= UH_PRIME_DEFAULT;
  }

  assert(h >= 0);
  return h;
}

std::vector<stag::DataPoint> stag::E2LSH::get_near_neighbors(const DataPoint& query) {
  std::vector<DataPoint> near_points;
  std::unordered_set<StagUInt> near_indices;

  for(StagUInt l = 0; l < parameterL; l++){
    StagUInt this_lsh = compute_lsh(l, query);

    if (!hashTables[l].contains(this_lsh)) {
      continue;
    } else {
      for (StagUInt candidatePIndex : hashTables[l][this_lsh]) {
        DataPoint& candidatePoint = points[candidatePIndex];
        if (near_indices.find(candidatePIndex) == near_indices.end()) {
          near_points.push_back(candidatePoint);
          near_indices.insert(candidatePIndex);
        }
      }
    }
  }

  return near_points;
}

StagReal stag::E2LSH::collision_probability(StagUInt K, StagUInt L,
                                            StagReal distance) {
  StagReal pc = stag::LSHFunction::collision_probability(distance);
  return 1 - pow(1.0 - pow(pc, (double) K), (double) L);
}

StagReal stag::E2LSH::collision_probability(StagReal distance) {
  return stag::E2LSH::collision_probability(parameterK, parameterL, distance);
}
