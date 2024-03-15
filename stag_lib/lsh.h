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
   * @brief A Euclidean locality-sensitive hash function.
   *
   * A function drawn at random from the standard family of Euclidean
   * locality-sensitive hash functions.
   * This function is defined by a random vector \f$a \in \mathbb{R}^d\f$
   * drawn from a Gaussian distribution, and a random offset
   * \f$b \in [0, 4]\f$ drawn uniformly at random.
   * Applying the function to some data point \f$x \in \mathbb{R}^d\f$ is
   * equivalent to computing
   *
   * \f[
   *    h = \left\lfloor \frac{\langle a, x \rangle + b}{4} \right\rfloor.
   * \f]
   *
   * Given two data points \f$x_1\f$ and \f$x_2\f$, the probability that they
   * are hashed to the same value is given by
   *
   * \f[
   *    p\left(\|x_1 - x_2\|_2\right) = \int_0^4 \frac{1}{\|x_1 - x_2\|_2} f\left(\frac{t}{\|x_1 - x_2\|_2}\right)\left(1 - \frac{t}{4}\right) \mathrm{dt},
   * \f]
   *
   * where \f$f(\cdot)\f$ is the probability density function of the Gaussian
   * distribution. The stag::LSHFunction::collision_probability function
   * computes this value.
   */
  class LSHFunction {
  public:
    /**
     * Initialise a random LSH function with the given dimension.
     *
     * @param dimension the dimensionality of the data
     */
    explicit LSHFunction(StagUInt dimension);

    /**
     * Apply this hash function to the given data point.
     *
     * @param point a data point to be hashed
     * @return a positive integer indicating which hash bucket this data point
     *         was hashed to
     */
    StagUInt apply(const DataPoint& point);

    /**
     * For two points at a given distance \f$c\f$, compute the probability that
     * they will collide in a random Euclidean LSH function.
     *
     * This probability is given by
     *
     * \f[
     *    p(c) = \int_0^4 \frac{1}{c} f\left(\frac{t}{c}\right)\left(1 - \frac{t}{4}\right) \mathrm{dt},
     * \f]
     *
     * where \f$f(\cdot)\f$ is the probability density function of the Gaussian
     * distribution.
     * This is equivalent to
     *
     * \f[
     *   p(c)=-\frac{1}{2\sqrt{2\pi}}\left( c e^{-\frac{8}{c^2}} \right) \left( e^{\frac{8}{c^2}} -1 \right) + \mathrm{erf}\left(\frac{2\sqrt{2}}{c}\right),
     * \f]
     *
     * where \f$\mathrm{erf}(\cdot)\f$ is the [error function](https://en.wikipedia.org/wiki/Error_function).
     *
     * @param distance the distance \f$c\f$.
     * @return the collision probability of two points at distance \f$c\f$.
     */
    static StagReal collision_probability(StagReal distance);

  private:
    std::vector<StagReal> a;
    StagReal b;
    StagUInt dim;
  };

  /**
   * \cond
   */
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
