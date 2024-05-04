/*
   This file is provided as part of the STAG library and released under the GPL
   license.
*/
/**
 * @file definitions.h
 * \brief Definitions of data types used throughout the library.
 */

#ifndef STAG_LIBRARY_DEFINITIONS_H
#define STAG_LIBRARY_DEFINITIONS_H

#include <Eigen/Sparse>
#include <cstdint>

/**
 * \cond
 */
// Define some mathematical operations.
#define MAX(a, b) ((a) > (b) ? (a) : (b))
#define MIN(a, b) ((a) < (b) ? (a) : (b))
#define SQR(a) ((a) * (a))
/**
 * \endcond
 */

/**
 * The integer type used throughout the library.
 */
typedef int64_t StagInt;

/**
 * The unsigned integer type used throughout the library.
 */
typedef size_t StagUInt;

/**
 * The real type used throughout the library.
 */
typedef double StagReal;

/**
 * The fundamental datatype used in this library is the sparse matrix.
 * We use the `Eigen::SparseMatrix` class in column-major format.
 */
typedef Eigen::SparseMatrix<StagReal, Eigen::ColMajor, StagInt> SprsMat;

/**
 *  The Dense Matrix data type used in the STAG library.
 */
typedef Eigen::Matrix<StagReal, Eigen::Dynamic,
                      Eigen::Dynamic, Eigen::RowMajor> DenseMat;

/**
 * An Eigen::Triplet representing an edge in a graph. Stores the row, column, and value
 * of an entry in a graph adjacency matrix.
 */
typedef Eigen::Triplet<StagReal, StagInt> EdgeTriplet;

/**
 * \cond
 */
// Redefine the eigen index type to be the same as StagInt
#undef EIGEN_DEFAULT_DENSE_INDEX_TYPE
#define EIGEN_DEFAULT_DENSE_INDEX_TYPE StagInt

// Define a small epsilon
#define EPSILON 0.0000000001
/**
 * \endcond
 */

#endif //STAG_LIBRARY_DEFINITIONS_H
