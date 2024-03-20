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

#include <definitions.h>
#include <lsh.h>

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
    CKNSGaussianKDEHashUnit(StagReal a, DenseMat* data, StagInt log_nmu, StagInt j);
    StagReal query(const stag::DataPoint& q);

  private:
    bool below_cutoff;
    stag::E2LSH LSH_buckets;
    StagInt j;
    StagInt log_nmu;
    StagReal a;

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
   * For more information on the Kernel Density Estimation problem, see kde.h.
   *
   * \par References
   * Charikar, Moses, et al. "Kernel density estimation through density
   * constrained near neighbor search." 2020 IEEE 61st Annual Symposium on
   * Foundations of Computer Science (FOCS). IEEE, 2020.
   */
  class CKNSGaussianKDE {
  public:
    /**
     * Initialise a new KDE data structure with the given dataset.
     *
     * Data should be a pointer to a matrix \f$X \in \mathbb{R}^{n \times d}\f$
     * where each row represents a data point.
     *
     * \note
     * The calling code is responsible for the memory management of the data
     * matrix, and it must be available throughout the life of the CKNS data
     * structure.
     *
     * @param data pointer to an \f$(n \times d)\f$ matrix containing the dataset.
     * @param a the parameter \f$a\f$ of the Gaussian kernel function.
     * @param eps the error parameter \f$\epsilon\f$ of the KDE data structure.
     */
    CKNSGaussianKDE(DenseMat* data, StagReal a, StagReal eps);

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
    StagReal query(stag::DataPoint& q);

  private:
    StagInt add_hash_unit(StagInt log_nmu_iter,
                          StagInt log_nmu,
                          StagInt iter,
                          StagInt j,
                          DenseMat* data);

    std::vector<std::vector<std::vector<CKNSGaussianKDEHashUnit>>> hash_units;
    std::mutex hash_units_mutex;
    StagInt max_log_nmu;
    StagInt num_log_nmu_iterations;
    StagInt n;
    StagReal a;
    StagInt k1;
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
     * Initialise the data structure with the given dataset and Gaussian kernel
     * parameter \f$a\f$.
     *
     * @param data
     * @param a
     */
    ExactGaussianKDE(DenseMat* data, StagReal a);

    /**
     * Calculate the KDE value for each of the data points in
     * the query matrix.
     *
     * @param q
     * @return a vector of KDE values.
     */
    std::vector<StagReal> query(DenseMat* q);

  private:
    std::vector<stag::DataPoint> all_data;
  };
}


#endif //STAG_LIBRARY_KDE_H
