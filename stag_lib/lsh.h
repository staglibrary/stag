/*
   Copyright (c) 2004-2005 Massachusetts Institute of Technology.
   Authors: Alexandr Andoni (andoni@mit.edu), Piotr Indyk (indyk@mit.edu)

   Modified by Peter Macgregor, 2024.

   This file is provided as part of the STAG library and released under the GPL
   license.
*/
/**
 * @file lsh.h
 * \brief Provides an implementation of Euclidean locality-sensitive hashing.
 */
#ifndef STAG_LIBRARY_LSH_H
#define STAG_LIBRARY_LSH_H

#include <vector>

#include <definitions.h>
#include <LSHTable.h>

/**
 * \cond
 */
// The value for algorithm parameter W.
#define LSH_PARAMETER_W 4.0
/**
 * \endcond
 */

namespace stag {

  /**
   * \brief A data point in d-dimensional space.
   *
   * A point is defined by its dimension, and a pointer to a C-style array
   * of coordinates. This structure uses C arrays in order that the calling code
   * can use any structure desired to store the underlying data.
   *
   * The calling code is always responsible for the underlying data, and
   * ensuring that for the life of the DataPoint object, the underlying data
   * array does not move.
   *
   * For example, the data can be stored in C++ vectors or an Eigen matrix.
   */
  class DataPoint {
  public:
    /**
     * Initialise a data point with an explicit dimension, and a pointer to the
     * data array.
     *
     * @param d the dimension of the data point
     * @param coords a pointer to the array of coordinates
     */
    DataPoint(StagUInt d, StagReal* coords) : dimension(d), coordinates(coords) {};

    /**
     * Initialise a data point to point to a given row of a dense matrix.
     *
     * @param all_data a dense matrix containing a full data set
     * @param row_index the index of the row containing this data point
     */
    DataPoint(DenseMat& all_data, StagInt row_index);

    /**
     * Initialise a data point to point to the data array of a given vector.
     *
     * @param point_vector a vector containing one data point.
     */
    DataPoint(std::vector<StagReal>& point_vector);

    /**
     * The dimension of the data point.
     */
    StagUInt dimension;

    /**
     * A pointer to a C-style array containing the data point.
     */
    StagReal *coordinates;
  };

  /**
   * \cond
   */
  // A function drawn from the locality-sensitive family of hash functions.
  class LSHFunction {
  public:
    explicit LSHFunction(StagUInt dimension);

    StagUInt apply(const DataPoint& point);

  private:
    std::vector<StagReal> a;
    StagReal b;
    StagUInt dim;
  };

  typedef struct RNNParametersT {
    StagUInt dimension; // dimension of points.
    StagReal parameterR2; // = parameterR^2
    bool checkDistance; // whether to check the query distance is less than R
    StagUInt parameterK; // parameter K of the algorithm.
    StagUInt parameterL; // parameter L of the algorithm.
  } RNNParametersT;


  class E2LSH {
  public:
    E2LSH() {};

    E2LSH(RNNParametersT algParameters, StagUInt nPoints,
          std::vector<DataPoint>& dataSet);

    std::vector<DataPoint> get_near_neighbors(const DataPoint& query);

  private:
    void initialise_fields_from_parameters(RNNParametersT algParameters,
                                           StagUInt nPointsEstimate);
    void initialise_hash_functions();

    std::vector<StagUInt> compute_lsh(StagUInt gNumber, const DataPoint& point);

    StagUInt dimension; // dimension of points.
    StagUInt parameterK; // parameter K of the algorithm.
    StagUInt parameterL; // parameter L of the algorithm.
    StagReal parameterR2; // = parameterR^2

    // Whether to check the distance to the query point is less than R when
    // returning points from the hash.
    bool checkDistanceWhenReturning;

    // The array of pointers to the points that are contained in the
    // structure. Some types of this structure (of UHashStructureT,
    // actually) use indices in this array to refer to points (as
    // opposed to using pointers).
    std::vector<DataPoint> points;

    // This table stores the LSH functions. There are <nHFTuples> rows
    // of <hfTuplesLength> LSH functions.
    std::vector<std::vector<LSHFunction>> lshFunctions;

    // The set of non-empty buckets
    std::vector<LSHTable> hashTables;

    // ***
    // The following vectors are used only for temporary operations
    // within this R-NN structure during a query operation.
    // ***

    // This vector is used for storing marked points in a query
    // operation (for computing distances to a point at most once). If
    // markedPoints[i]=TRUE then point <i> was examined already.
    std::vector<bool> markedPoints;

    // This vector stored the indices in the vector <markedPoints> of all
    // TRUE entries.
    std::vector<StagUInt> markedPointsIndices;

    // the size of <markedPoints> and of <markedPointsIndeces>
    StagUInt sizeMarkedPoints;
  };

  /**
   * \endcond
   */
}

#endif //STAG_LIBRARY_LSH_H
