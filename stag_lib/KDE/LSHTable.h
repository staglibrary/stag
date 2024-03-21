/*
   Copyright (c) 2004-2005 Massachusetts Institute of Technology.
   Authors: Alexandr Andoni (andoni@mit.edu), Piotr Indyk (indyk@mit.edu)

   Modified by Peter Macgregor, 2024.

   This file is provided as part of the STAG library and released under the GPL
   license.
*/
#ifndef STAG_LIBRARY_LSHTABLE_H
#define STAG_LIBRARY_LSHTABLE_H

#include "definitions.h"

/**
 * \cond
 * Do not document anything in this file - it is used internally for the LSH
 * functionality only.
 */

// 4294967291 = 2^32-5
#define UH_PRIME_DEFAULT 4294967291U

// 2^29
#define MAX_HASH_RND 536870912U

// 2^32-1
#define TWO_TO_32_MINUS_1 4294967295U

namespace stag {
  typedef struct BucketHashingIndexT {
    StagUInt main_index;
    StagUInt control_index;
    BucketHashingIndexT(StagUInt main, StagUInt control)
        : main_index(main), control_index(control) {};
    BucketHashingIndexT() : main_index(0), control_index(0) {};
  } BucketHashingIndexT;

  // The type definition for an LSH bucket. A bucket is a container for points
  // which has all been hashed to the same value for some hash function g. The
  // function g is a vector of K LSH functions.
  //
  // In the bucket array, the points are referred to only as indices in the larger
  // LSH structure.
  class LSHBucket {
  public:
    explicit LSHBucket(StagUInt cv) : controlValue1(cv) {};

    void add_point(StagInt new_point) {
      points.push_back(new_point);
    }

    // These controlValues are used instead of the full k-vector (value
    // of the hash function g) describing the bucket. With a high
    // probability all buckets will have different pairs of
    // controlValues.
    StagUInt controlValue1;

    // We store only the point indices in the hash buckets - these indices
    // correspond to the indices of the points in the overall hashtable structure.
    std::vector<StagInt> points;
  };

  // A universal hash table with collision solved by chaining. The
  // chains and the buckets are stored using either singly linked lists
  // or static arrays (depending on the value of the field <typeHT>).
  class LSHTable {
  public:
    LSHTable() {};

    // Creates a new hash table (initializes the hash table and the hash
    // functions used).
    LSHTable(StagUInt hashTableSize, StagUInt bucketVectorLength);

    // Adds the bucket entry (a point <point>) to the bucket defined by
    // the given index. If no such bucket exists, then it is first created.
    void add_bucket_entry(const BucketHashingIndexT& bucketIndex,
                          StagInt pointIndex);

    // Returns the bucket defined by the given bucket index.
    LSHBucket* get_bucket(const BucketHashingIndexT& bucketIndex);

    BucketHashingIndexT compute_bucket_index(
        const std::vector<StagUInt>& uVector);

  private:
    // The array containing the hash slots of the universal hashing.
    std::vector<std::vector<LSHBucket>> hashTable;

    // The size of hashTable.
    StagUInt tableSize;

    StagUInt nHashedPoints;

    StagUInt prime; // the prime used for the universal hash functions.
    StagUInt hashedDataLength;// the number of IntT's in an element from U (U is the set of values to hash).

    // The hash functions used for the universal hashing.

    // The main hash function (that defines the index
    // of the slot in the table).
    // The type of the hash function is: h_{a}(k) = ((a\cdot k)mod p)mod hashTableSize.
    std::vector<StagUInt> mainHashA;

    // Control hash functions: used to compute/check the <controlValue>s
    // of <GBucket>s.
    // The type of the hash function is: h_{a}(k) = (a\cdot k)mod p
    std::vector<StagUInt> controlHash1;

  };
}

/**
 * \endcond
 */

#endif //STAG_LIBRARY_LSHTABLE_H
