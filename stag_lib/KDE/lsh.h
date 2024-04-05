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
 *
 * Locality sensitive hashing is a primitive used for near-neighbour search.
 * In particular, a locality sensitive hash function (see stag::LSHFunction)
 * hashes vectors into buckets such that two vectors are more likely to be hashed
 * to the same bucket if their Euclidean distance is small.
 *
 * The stag::E2LSH hash table implements a full approximate-near neighbour
 * datastructure based on basic Euclidean locality sensitive hash functions.
 *
 * \par Reference
 * Andoni, Alexandr, and Piotr Indyk. "Near-optimal hashing algorithms for approximate nearest neighbor in high
 * dimensions." Communications of the ACM 51.1 (2008): 117-122.
 */
#ifndef STAG_LIBRARY_LSH_H
#define STAG_LIBRARY_LSH_H

#include <vector>

#include <unordered_map>
#include "definitions.h"
#include "data.h"

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
   *
   * Typical STAG users will use the stag::E2LSH hash table rather than these
   * the LSHFunction class directly.
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
     * @return an integer indicating which hash bucket this data point
     *         was hashed to
     */
    StagInt apply(const DataPoint& point);

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
  class MultiLSHFunction {
  public:
    MultiLSHFunction(StagInt dimension, StagInt num_functions);
    StagInt apply(const stag::DataPoint& point);
  private:
    StagInt L;
    DenseMat rand_proj;
    Eigen::VectorXd rand_offset;
    Eigen::Matrix<StagInt, Eigen::Dynamic, 1> uhash_vector;
  };
  /**
   * \endcond
   */

  /**
   * \brief A Euclidean locality sensitive hash table.
   *
   * The E2LSH hash table is constructed with some set of data points, which are
   * hashed with several copies of the stag::LSHFunction.
   *
   * Then, for any query point, the data structure returns the points in the
   * original dataset which are close to the query point.
   * The probability that a given point \f$x\f$ in the data set is returned for
   * query \f$q\f$ is dependent on the distance between \f$q\f$ and \f$x\f$.
   *
   * The E2LSH hash table takes two parameters, K and L, which control the
   * probability that two points will collide in the hash table.
   * For query point \f$q\f$, a data point \f$x\f$ at distance \f$c \in \mathbb{R}\f$
   * from \f$q\f$ is returned with probability
   * \f[
   *    1 - (1 - p(c)^K)^L,
   * \f]
   * where \f$p(c)\f$ is the probability that a single stag::LSHFunction will
   * hash \f$q\f$ and \f$x\f$ to the same value.
   * This probability can be computed with the stag::E2LSH::collision_probability
   * method.
   *
   * Larger values of K and L will increase both the construction and query time
   * of the hash table.
   */
  class E2LSH {
  public:
    E2LSH() {};

    /**
     * Initialise the E2LSH hash table.
     *
     * @param K parameter K of the hash table
     * @param L parameter L of the hash table
     * @param dataSet a pointer to the dataSet to be hashed into the hash
     *                table. The actual data should be stored and controlled by
     *                the calling code, and this vector of data point pointers
     *                will be used by the LSH table.
     */
    E2LSH(StagUInt K,
          StagUInt L,
          std::vector<DataPoint>& dataSet);

    /**
     * Query the LSH table to find the near neighbors of a given query point.
     *
     * Each point in the dataset will be returned with some probability
     * dependent on the distance to the query point and the parameters K and L.
     *
     * @param query the data point to be queried.
     * @return
     */
    std::vector<DataPoint> get_near_neighbors(const DataPoint& query);

    /**
     * Compute the probability that a data point at a given distance from a query
     * point will be returned by this hash table.
     *
     * @param distance the distance between a query point and data point
     * @return
     */
    StagReal collision_probability(StagReal distance);

    /**
     * For a given value of K and L, return the probability that a point at a
     * given distance from a query point will be returned by an E2LSH table
     * with the given parameters.
     *
     * @param K the parameter K of the ELSH table
     * @param L the parameter L of the E2LSH table
     * @param distance the distance between a query point and data point
     * @return
     */
    static StagReal collision_probability(StagUInt K, StagUInt L,
                                          StagReal distance);

  private:
    void initialise_hash_functions();

    StagInt compute_lsh(StagUInt gNumber, const DataPoint& point);

    StagUInt dimension; // dimension of points.
    StagUInt parameterK; // parameter K of the algorithm.
    StagUInt parameterL; // parameter L of the algorithm.

    std::vector<StagInt> rnd_vec; // used for hashing vectors

    // The array of pointers to the points that are contained in the
    // structure. Some types of this structure (of UHashStructureT,
    // actually) use indices in this array to refer to points (as
    // opposed to using pointers).
    std::vector<DataPoint> points;

    // This table stores the LSH functions. There are <nHFTuples> rows
    // of <hfTuplesLength> LSH functions.
    std::vector<MultiLSHFunction> lshFunctions;

    // The set of non-empty buckets
    std::vector<std::unordered_map<StagInt,std::vector<StagUInt>>> hashTables;
  };
}

#endif //STAG_LIBRARY_LSH_H
