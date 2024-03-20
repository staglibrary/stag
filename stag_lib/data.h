/*
   This file is provided as part of the STAG library and released under the GPL
   license.
*/

/**
 * @file data.h
 * \brief Methods for dealing with datasets.
 *
 * This file contains methods for reading and writing data to disk, and for
 * converting data between different formats.
 */

#ifndef STAG_LIBRARY_DATA_H
#define STAG_LIBRARY_DATA_H

#include <definitions.h>

namespace stag{
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
   * Convert data in an eigen matrix to an array of data point pointers.
   *
   * @param data the eigen matrix containing the data
   * @return a vector of stag::DataPoint pointing to the rows of the matrix
   */
  std::vector<stag::DataPoint> matrix_to_datapoints(DenseMat* data);
}

#endif //STAG_LIBRARY_DATA_H
