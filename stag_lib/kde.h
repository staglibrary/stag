/*
   This file is provided as part of the STAG library and released under the GPL
   license.
*/

/**
 * @file kde.h
 * \brief Methods for computing approximate kernel density estimation.
 *
 * Given some *kernel function*
 * \f$k: \mathbb{R}^d \times \mathbb{R}^d \rightarrow \mathbb{R}\f$,
 * a set of *data points* \f$x_1, \ldots, x_n \in \mathbb{R}^d\f$,
 * and a *query point* \f$q \in \mathbb{R}^d\f$,
 * the *kernel density* of \f$q\f$ is given by
 *
 * \f[
 *    K(q) = \frac{1}{n} \sum_{i = 1}^n k(q, x_i).
 * \f]
 *
 * A common kernel function is the Gaussian kernel function
 * \f[
 *    k(u, v) = \exp\left(- a \|u - v\|_2^2\right),
 * \f]
 * where \f$a \geq 0\f$ is a parameter controlling the 'bandwidth' of the
 * kernel.
 *
 * Computing the kernel density for a query point exactly requires computing the
 * distance from the query point to every data point, requiring
 * \f$\Omega(n d)\f$ time.
 * This motivates the study of *kernel density estimation*, in which the goal
 * is to estimate the kernel density within some error tolerance, in faster
 * time than computing it exactly.
 * Specifically, given some error parameter \f$\epsilon\f$, a kernel density
 * estimation algorithm will return \f$\mathrm{KDE}(q)\f$ for some query point
 * \f$q\f$ such that
 * \f[
 *  (1 - \epsilon) K(q) \leq \mathrm{KDE}(q) \leq (1 + \epsilon) K(q).
 * \f]
 *
 * This module provides the stag::CKNSGaussianKDE data structure which takes
 * \f$O(\epsilon^{-1} n^{1.25})\f$ time for initialisation, and can then provide
 * KDE estimates in time \f$O(\epsilon^{-2} n^{0.25})\f$ for each query.
 */

#ifndef STAG_LIBRARY_KDE_H
#define STAG_LIBRARY_KDE_H

#include <mutex>

#include "definitions.h"
#include "lsh.h"

namespace stag {

  /**
   * Compute the Gaussian kernel similarity between the points u and v.
   *
   * Given a parameter \f$a \geq 0\f$ and points \f$u, v \in \mathbb{R}^n\f$,
   * the Gaussian kernel similarity between \f$u\f$ and \f$v\f$ is given by
   *
   * \f[
   *    k(u, v) = \exp\left( - a \|u - v\|^2_2 \right).
   * \f]
   *
   * Note that the Gaussian kernel is sometimes parameterised by \f$\sigma^2\f$,
   * which is related to our parameter \f$a\f$ by
   *
   * \f[
   *    a = \frac{1}{\sigma^2}.
   * \f]
   *
   * @param a the parameter a in the Gaussian kernel.
   * @param u a data point \f$u\f$
   * @param v a data point \f$v\f$
   * @return the Gaussian kernel similarity between \f$u\f$ and \f$v\f$.
   */
  StagReal gaussian_kernel(StagReal a,
                           const stag::DataPoint& u,
                           const stag::DataPoint& v);

  /**
   * Compute the Gaussian kernel similarity for two points at a squared distance
   * \f$c\f$.
   *
   * Given a parameter \f$a \geq 0\f$, the Gaussian kernel similarity between
   * two points at distance \f$c\f$ is given by
   *
   * \f[
   *    \exp\left( - a c \right).
   * \f]
   *
   * @param a the parameter a in the Gaussian kernel.
   * @param c the squared distance between two points.
   * @return
   */
  StagReal gaussian_kernel(StagReal a, StagReal c);

  /**
   * \cond
   * A helper class for the CKNSGaussianKDE data structure. Undocumented.
   */
  class CKNSGaussianKDEHashUnit {
  public:
    CKNSGaussianKDEHashUnit(StagReal a, DenseMat* data, StagInt log_nmu,
                            StagInt j, StagReal K2_constant, StagInt prob_offset);
    StagReal query(const stag::DataPoint& q);

  private:
    StagReal query_neighbors(const stag::DataPoint& q,
                             const std::vector<stag::DataPoint>& neighbors);
    bool below_cutoff;
    bool final_shell;
    stag::E2LSH LSH_buckets;
    StagInt j;
    StagInt J;
    StagInt log_nmu;
    StagReal a;
    StagInt sampling_offset;
    StagInt n;

    // Used only if the number of data points is below the cutoff.
    std::vector<stag::DataPoint> all_data;
  };
  /**
   * \endcond
   */

  /**
   * \brief A CKNS Gaussian KDE data structure.
   *
   * This data structure implements the CKNS algorithm for kernel density
   * estimation. Given data \f$x_1, \ldots, x_n \in \mathbb{R}^d\f$, in matrix
   * format, this data structure will preprocess the data in
   * \f$O(\epsilon^{-2} n^{1.25})\f$ time, such that for any query point,
   * a \f$(1 + \epsilon)\f$-approximate kernel density estimate can be returned
   * in \f$O(\epsilon^{-2} n^{0.25})\f$ time.
   *
   * \par References
   * Charikar, Moses, et al. "Kernel density estimation through density
   * constrained near neighbor search." 2020 IEEE 61st Annual Symposium on
   * Foundations of Computer Science (FOCS). IEEE, 2020.
   */
  class CKNSGaussianKDE {
  public:
    /**
     * \cond
     * Do not document the default constructor. It should only be used as a
     * placeholder.
     */
    CKNSGaussianKDE() {}
    /**
     * \endcond
     */

    /**
     * Initialise a new KDE data structure with the given dataset.
     *
     * Data should be a pointer to a matrix \f$X \in \mathbb{R}^{n \times d}\f$
     * where each row represents a data point.
     *
     * The \f$\epsilon\f$ parameter is used to control the error guarantee of the
     * CKNS data structure. A lower value of \f$\epsilon\f$ will give a more accurate
     * estimate, at the cost of a higher processing time.
     * The data structure should produce estimates which are within a
     * \f$(1 \pm \epsilon)\f$ factor of the true kernel density.
     *
     * The initialisation time complexity of the data structure is
     * \f$O(\epsilon^{-2} n^{1.25} \log^2(n))\f$ and the query time for each
     * query point is \f$(O(\epsilon^{-2} n^{0.25} \log^2(n))\f$.
     *
     * @param data pointer to an \f$(n \times d)\f$ matrix containing the dataset.
     * @param a the parameter \f$a\f$ of the Gaussian kernel function.
     * @param eps (optional) the error parameter \f$\epsilon\f$ of the KDE data
     *            structure. Default is 0.5.
     * @param min_mu (optional) the minimum kernel density value of any query
     *               point. A smaller number will give longer preprocessing and
     *               query time complexity. If a query point has a kernel density
     *               smaller than this value, then the data structure may not
     *               return the correct result.
     *               Default is 1 / n.
     */
    CKNSGaussianKDE(DenseMat* data, StagReal a, StagReal eps, StagReal min_mu);

    /**
     * \overload
     */
    CKNSGaussianKDE(DenseMat* data, StagReal a, StagReal eps);

    /**
     * \overload
     */
    CKNSGaussianKDE(DenseMat* data, StagReal a);

    /**
     * Initialise a new KDE data structure with the given dataset.
     *
     * Data should be a pointer to a matrix \f$X \in \mathbb{R}^{n \times d}\f$
     * where each row represents a data point.
     *
     * This constructor gives fine-grained control over the constants used by
     * the data structure to control the accuracy and variance of the estimator.
     * For casual application, the simpler constructors will be sufficient, and
     * will select sensible default values for the constants.
     *
     * @param data pointer to an \f$(n \times d)\f$ matrix containing the dataset.
     * @param a the parameter \f$a\f$ of the Gaussian kernel function.
     * @param min_mu the minimum kernel density value of any query
     *               point. A smaller number will give longer preprocessing and
     *               query time complexity. If a query point has a kernel density
     *               smaller than this value, then the data structure may not
     *               return the correct result.
     * @param K1 the number of copies of the data structure to create in parallel.
     *           This parameter controls the variance of the estimator returned
     *           by the algorithm. Is is usually set to \f$\epsilon^{-2} \cdot \log(n)\f$.
     * @param K2_constant controls the collision probability of each of the E2LSH
     *                    hash tables used within the data structure. A higher value
     *                    will give more accurate estimates at the cost of higher memory
     *                    and time complexity. It is usually set to \f$0.1 \log(n)\f$.
     * @param sampling_offset the CKNS algorithm samples the dataset with various
     *                        sampling probabilities. Setting a sampling offset of
     *                        \f$k\f$ will further subsample the data by a factor
     *                        of \f$1/2^k\f$. This will speed up the algorithm
     *                        at the cost of some accuracy.
     *                        It is usually set to \f$0\f$.
     */
    CKNSGaussianKDE(DenseMat* data, StagReal a, StagReal min_mu, StagInt K1,
                    StagReal K2_constant, StagInt sampling_offset);

    /**
     * \cond
     * Hidden constructor, used to specify a minimuum and maximum id of the
     * provided matrix.
     */
    CKNSGaussianKDE(DenseMat* data, StagReal a, StagReal min_mu, StagInt K1,
                    StagReal K2_constant, StagInt sampling_offset, StagInt min_idx,
                    StagInt max_idx);
    /**
     * \nocond
     */

    /**
     * Calculate an estimate of the KDE value for each of the data points in
     * the query matrix.
     *
     * The KDE estimate for each row in the query matrix will be returned as
     * a vector.
     *
     * @param q a pointer to an \f$(m \times d)\f$ matrix containing the \f$m\f$
     *          query points.
     * @return a vector of length \f$m\f$ containing the KDE estimates
     */
    std::vector<StagReal> query(DenseMat* q);

    /**
     * Calculate the KDE estimate for the given query point.
     *
     * \note
     * If you would like to obtain a KDE query for many query points, you should
     * use the other query method to pass them all simultaneously.
     *
     * @param q the query data point
     * @return the KDE estimate for the given query point
     */
    StagReal query(const stag::DataPoint& q);

  private:
    void initialize(DenseMat* data, StagReal a, StagReal min_mu, StagInt K1,
                    StagReal K2_constant, StagInt prob_offset, StagInt min_idx,
                    StagInt max_idx);
    StagInt add_hash_unit(StagInt log_nmu_iter,
                          StagInt log_nmu,
                          StagInt iter,
                          StagInt j,
                          DenseMat* data,
                          std::mutex& units_mutex);
    std::vector<StagReal> chunk_query(DenseMat* query, StagInt chunk_start, StagInt chunk_end);

    std::vector<std::vector<std::vector<CKNSGaussianKDEHashUnit>>> hash_units;
    std::vector<DenseMat> data_copies;
    StagInt min_id;
    StagInt max_id;
    StagInt max_log_nmu;
    StagInt min_log_nmu;
    StagInt num_log_nmu_iterations;
    StagInt sampling_offset;
    StagInt n;
    StagReal a;
    StagInt k1;
    StagReal k2_constant;
  };

  /**
   * \brief A data structure for computing the exact Gauussian KDE.
   *
   * This data structure uses a brute-force algorithm to compute the kernel
   * density of each query point.
   *
   * The time complexity of initialisation with \f$n\f$ data points is \f$O(n)\f$.
   * The query time complexity is \f$O(m n d)\f$, where \f$m\f$ is the number
   * of query points, and \f$d\f$ is the dimensionality of the data.
   */
  class ExactGaussianKDE {
  public:
    /**
     * \cond
     * Do not document the default constructor. It should only be used as a
     * placeholder.
     */
    ExactGaussianKDE() {}
    /**
     * \endcond
     */

    /**
     * Initialise the data structure with the given dataset and Gaussian kernel
     * parameter \f$a\f$.
     *
     * The initialisation time for this data structure is \f$O(1)\f$.
     *
     * @param data
     * @param a
     */
    ExactGaussianKDE(DenseMat* data, StagReal a);


    /**
     * \cond
     * Compute the exact Gaussian KDE for indexes between the min and max
     * for the given data.
     *
     * Additionally, provide a threshold value, beneath which the Gaussian
     * kernel is taken to be 0.
     */
    ExactGaussianKDE(DenseMat* data, StagReal a, StagInt min_idx, StagInt max_idx);
    /**
     * \endcond
     */

    /**
     * Calculate the KDE value for each of the data points in
     * the query matrix.
     *
     * The query time complexity for each query point is \f$O(n)\f$.
     *
     * @param q
     * @return a vector of KDE values.
     */
    std::vector<StagReal> query(DenseMat* q);

    /**
     * Calculate the exact kernel density for the given query point.
     *
     * \note
     * If you would like to obtain a kernel density for many query points, you
     * should use the other query method to pass them all simultaneously.
     *
     * @param q the query data point
     * @return the kernel density for the given query point
     */
    StagReal query(const stag::DataPoint& q);

    /**
     * \cond
     * The sample neighbor method is used only by the approximate similarity
     * graph code - we do not expect end-users to use it.
     */
     std::vector<StagInt> sample_neighbors(const stag::DataPoint& q, StagReal degree,
                                           std::vector<StagReal> rs);
     /**
      * \endcond
      */

  private:
    std::vector<stag::DataPoint> all_data;
    StagReal a;
    StagInt min_id;
    StagInt max_id;
  };
}


#endif //STAG_LIBRARY_KDE_H
