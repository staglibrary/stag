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

// Generate a random integer in the range [rangeStart,
// rangeEnd]. Inputs must satisfy: rangeStart <= rangeEnd.
StagInt genRandomInt(StagInt rangeStart, StagInt rangeEnd){
  assert(rangeStart <= rangeEnd);

  StagInt r;
  r = rangeStart + (StagInt)((rangeEnd - rangeStart + 1.0) * random() / (RAND_MAX + 1.0));

  assert(r >= rangeStart && r <= rangeEnd);
  return r;
}

// Generate a random 32-bits unsigned (Uns32T) in the range
// [rangeStart, rangeEnd]. Inputs must satisfy: rangeStart <=
// rangeEnd.
StagUInt genRandomUInt(StagUInt rangeStart, StagUInt rangeEnd){
  assert(rangeStart <= rangeEnd);
  assert(RAND_MAX >= rangeEnd - rangeStart);

  StagUInt r;
  r = rangeStart + (StagUInt)((rangeEnd - rangeStart + 1.0) * random() / (RAND_MAX + 1.0));

  assert(r >= rangeStart && r <= rangeEnd);
  return r;
}

// Generate a random real distributed uniformly in [rangeStart,
// rangeEnd]. Input must satisfy: rangeStart <= rangeEnd. The
// granularity of generated random reals is given by RAND_MAX.
StagReal genUniformRandom(StagReal rangeStart, StagReal rangeEnd){
  assert(rangeStart <= rangeEnd);

  StagReal r;
  r = rangeStart + ((rangeEnd - rangeStart) * (StagReal)random() / (StagReal)RAND_MAX);

  assert(r >= rangeStart && r <= rangeEnd);
  return r;
}

// Generate a random real from normal distribution N(0,1).
StagReal genGaussianRandom(){
  // Use Box-Muller transform to generate a point from normal
  // distribution.
  StagReal x1, x2;
  do{
    x1 = genUniformRandom(0.0, 1.0);
  } while (x1 == 0); // cannot take log of 0.
  x2 = genUniformRandom(0.0, 1.0);
  StagReal z;
  z = sqrt(-2.0 * log(x1)) * cos(2.0 * M_PI * x2);
  return z;
}

// Generate a random real from Cauchy distribution N(0,1).
StagReal genCauchyRandom(){
  StagReal x, y;
  x = genGaussianRandom();
  y = genGaussianRandom();
  if (abs(y) < 0.0000001) {
    y = 0.0000001;
  }
  return x / y;
}

//------------------------------------------------------------------------------
// Implementation of the LSHBucket class
//------------------------------------------------------------------------------
void stag::LSHBucket::add_point(StagInt new_point) {
  points.push_back(new_point);
}

//------------------------------------------------------------------------------
// Implementation of the LSHTable class
//------------------------------------------------------------------------------

// Creates a new hash table (initializes the hash table and the hash
// functions used).
stag::LSHTable::LSHTable(StagUInt hashTableSize, StagUInt bucketVectorLength) {
  tableSize = hashTableSize;
  nHashedPoints = 0;

  prime = UH_PRIME_DEFAULT;
  hashedDataLength = bucketVectorLength;

  hashTable.resize(tableSize);

  // Initializing the main hash function.
  mainHashA.resize(hashedDataLength);
  for(StagUInt i = 0; i < hashedDataLength; i++){
    mainHashA[i] = genRandomUInt(1, MAX_HASH_RND);
  }

  // Initializing the control hash functions.
  controlHash1.resize(hashedDataLength);
  for(StagUInt i = 0; i < hashedDataLength; i++){
    controlHash1[i] = genRandomUInt(1, MAX_HASH_RND);
  }
}

// Computes (a.b)mod UH_PRIME_DEFAULT.
inline StagUInt compute_product_mod_prime(const std::vector<StagUInt>& a,
                                          const std::vector<StagUInt>& b){
  assert(a.size() == b.size());

  StagUInt h = 0;
  for(StagUInt i = 0; i < a.size(); i++){
    h = h + (StagUInt)a[i] * (StagUInt)b[i];
    h = (h & TWO_TO_32_MINUS_1) + 5 * (h >> 32);
    if (h >= UH_PRIME_DEFAULT) {
      h = h - UH_PRIME_DEFAULT;
    }
  }
  return h;
}

// Compute fuction ((rndVector . data)mod prime)mod hashTableSize
inline StagUInt compute_uhash_function(const std::vector<StagUInt>& rndVector,
                                       const std::vector<StagUInt>& data,
                                       StagUInt hashTableSize){
  assert(rndVector.size() == data.size());
  StagUInt h = compute_product_mod_prime(rndVector, data) % hashTableSize;

  assert(h < hashTableSize);
  return h;
}

// Adds the bucket entry (a point <point>) to the bucket defined by
// the given index. If no such bucket exists, then it is first created.
void stag::LSHTable::add_bucket_entry(const BucketHashingIndexT& bucketIndex,
                                      StagInt pointIndex){
  bool found = false;
  for (auto& b : hashTable[bucketIndex.main_index]) {
    if (b.controlValue1 == bucketIndex.control_index) {
      b.add_point(pointIndex);
      found = true;
      break;
    }
  }

  // If we didn't find a bucket with the correct control value, add one to
  // the end of the hash table vector.
  if (!found) {
    hashTable[bucketIndex.main_index].emplace_back(bucketIndex.control_index);
    hashTable[bucketIndex.main_index].back().add_point(pointIndex);
  }

  // Keep statistics on number of hashed points.
  nHashedPoints++;
}

// Returns the bucket defined by the given bucket index.
stag::LSHBucket* stag::LSHTable::get_bucket(const BucketHashingIndexT& bucketIndex) {
  for (StagUInt i = 0; i < hashTable[bucketIndex.main_index].size(); i++) {
    if (hashTable[bucketIndex.main_index][i].controlValue1 == bucketIndex.control_index) {
      return &hashTable[bucketIndex.main_index][i];
    }
  }

  // We have failed to find a bucket - calling code must check for this.
  return nullptr;
}

stag::BucketHashingIndexT stag::LSHTable::compute_bucket_index(
    const std::vector<StagUInt>& uVector){
  StagUInt main = compute_uhash_function(mainHashA, uVector, tableSize);
  StagUInt control = compute_product_mod_prime(controlHash1, uVector);
  return {main, control};
}

// Returns TRUE iff |p1-p2|_2^2 <= threshold
inline bool isDistanceSqrLeq(StagUInt dimension, const stag::LSHDataPointT& p1,
                             const stag::LSHDataPointT& p2, StagReal threshold){
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

StagUInt stag::LSHFunction::apply(const LSHDataPointT& point) {
  assert(point.dimension == dim);

  StagReal value = 0;
  for(StagUInt d = 0; d < dim; d++){
    value += point.coordinates[d] * a[d];
  }

  return (StagUInt) floor((value + b) / LSH_PARAMETER_W);
}


//------------------------------------------------------------------------------
// Implementation of the E2LSH class.
//------------------------------------------------------------------------------
stag::E2LSH::E2LSH(RNNParametersT algParameters, StagUInt nPoints,
                   std::vector<LSHDataPointT>& dataSet){
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

std::vector<StagUInt> stag::E2LSH::compute_lsh(StagUInt gNumber, const LSHDataPointT& point) {
  std::vector<StagUInt> result(parameterK);
  for(StagUInt i = 0; i < parameterK; i++){
    result[i] = lshFunctions[gNumber][i].apply(point);
  }
  return result;
}

std::vector<stag::LSHDataPointT> stag::E2LSH::get_near_neighbors(const LSHDataPointT& query) {
  std::vector<LSHDataPointT> near_points;

  StagUInt nMarkedPoints = 0;// the number of marked points
  for(StagUInt l = 0; l < parameterL; l++){
    BucketHashingIndexT bucket_index = hashTables[0].compute_bucket_index(
        compute_lsh(l, query));
    LSHBucket* bucket = hashTables[l].get_bucket(bucket_index);

    // circle through the bucket and add to <result> the points that are near.
    if (bucket != nullptr) {
      for (StagInt candidatePIndex : bucket->points) {
        LSHDataPointT candidatePoint = points[candidatePIndex];
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
