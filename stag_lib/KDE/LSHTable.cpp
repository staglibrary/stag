/*
   Copyright (c) 2004-2005 Massachusetts Institute of Technology.
   Authors: Alexandr Andoni (andoni@mit.edu), Piotr Indyk (indyk@mit.edu)

   Modified by Peter Macgregor, 2024.

   This file is provided as part of the STAG library and released under the GPL
   license.
*/

#include "LSHTable.h"

#include "Graph/random.h"

// Generate a random 32-bits unsigned (Uns32T) in the range
// [rangeStart, rangeEnd]. Inputs must satisfy: rangeStart <=
// rangeEnd.
StagUInt genRandomUInt(StagUInt rangeStart, StagUInt rangeEnd){
  assert(rangeStart <= rangeEnd);

  std::uniform_int_distribution<StagUInt> dist(rangeStart, rangeEnd);
  StagUInt r = dist(*stag::get_global_rng());

  assert(r >= rangeStart && r <= rangeEnd);
  return r;
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

// Adds the bucket entry (a point <point>) to the bucket defined by
// the given index. If no such bucket exists, then it is first created.
void stag::LSHTable::add_bucket_entry(
    const stag::BucketHashingIndexT& bucketIndex,
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
stag::LSHBucket* stag::LSHTable::get_bucket(
    const BucketHashingIndexT& bucketIndex) {
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
