/*
   This file is provided as part of the STAG library and released under the GPL
   license.
*/

#include "kde.h"

#include "Graph/random.h"
#include "multithreading/ctpl_stl.h"
#include "data.h"
#include <iostream>

#define TWO_ROOT_TWO 2.828427124
#define TWO_ROOT_TWOPI 5.0132565
#define LOG_TWO 0.69314718056

// The CKNS algorithm has a couple of 'magic constants' which control the
// probability guarantee and variance bounds.
#define K2_DEFAULT_CONSTANT 5       // K_2 = C log(n) p^{-k_j}
#define K1_DEFAULT_CONSTANT 1       // K_1 = C log(n) / eps^2
#define EPS_DEFAULT 0.5             // K_1 = C log(n) / eps^2

// At a certain number of sampled points, we might as well brute-force the hash
// unit.
#define HASH_UNIT_CUTOFF 1000

/*
 * Used to disable compiler warning for unused variable.
 */
template<class T> void ignore_warning(const T&){}

/**
 * Compute the squared distance between the given two data points.
 *
 * @param u the first data point
 * @param v the second data point
 * @return the squared distance between \f$u\f$ and \f$v\f$.
 */
StagReal squared_distance(const stag::DataPoint& u, const stag::DataPoint& v) {
  assert(u.dimension == v.dimension);
  StagReal result = 0;
  for (StagUInt i = 0; i < u.dimension; i++) {
    result += SQR(u.coordinates[i] - v.coordinates[i]);
  }
  return result;
}

StagReal stag::gaussian_kernel(StagReal a, StagReal c) {
  return exp(-(a * c));
}

StagReal stag::gaussian_kernel(StagReal a, const stag::DataPoint& u,
                         const stag::DataPoint& v) {
  return stag::gaussian_kernel(a, squared_distance(u, v));
}

//------------------------------------------------------------------------------
// Beginning of CKNS Implementation
//------------------------------------------------------------------------------
/**
 * Compute the value of J for the CKNS algorithm, given
 *   - n, the number of data points
 *   - log2(n * mu), as an integer (at most log2(n)).
 *
 * @param n the number of data points
 * @param log_nmu the value of log2(n * mu)
 * @return
 */
StagInt ckns_J(StagInt n, StagInt log_nmu) {
  assert(log_nmu < log2(n));
  return ((StagInt) ceil(log2((StagReal) n))) - log_nmu;
}

/**
 * Compute the sampling probability for the CKNS algorithm, given a value of
 * j and log2(n * mu).
 *
 * @param j
 * @param log_nmu the value of log2(n * mu)
 * @return
 */
StagReal ckns_p_sampling(StagInt j, StagInt log_nmu) {
  return pow(2, (StagReal) -j) * pow(2, (StagReal) -log_nmu);
}

StagReal ckns_gaussian_rj_squared(StagInt j, StagReal a) {
  return (StagReal) j * LOG_TWO / a;
}

/**
 * Create the LSH parameters for a particular application of E2LSH as part of
 * the CKNS algorithm with the Gaussian kernel.
 *
 * @param J the total number of sampling levels
 * @param j the current sampling level
 * @param n the total number of data points
 * @param d the dimension of the data
 * @param a the scaling parameter of the Gaussian kernel
 * @return the K and L parameters for this E2LSH table.
 */
std::vector<StagUInt> ckns_gaussian_create_lsh_params(
    StagInt J, StagInt j, StagReal a, StagReal K2_constant) {
  StagReal r_j = sqrt((StagReal) j * log(2) / a);
  StagReal p_j = stag::LSHFunction::collision_probability(r_j);
  StagReal phi_j = ceil((((StagReal) j)/((StagReal) J)) * (StagReal) (J - j + 1));
  StagUInt k_j = MAX(1, floor(- phi_j / log2(p_j)));
  StagUInt K_2 = ceil(K2_constant * pow(2, phi_j));
  return {
      k_j, // parameter K
      K_2, // parameter L
  };
}

//------------------------------------------------------------------------------
// Implementation of the CKNS Gaussian KDE Hash Unit
//
// The CKNS KDE data structure is made up of several E2LSH hashes of different
// subsets of the dataset. We define a data structure (class since we're in C++
// land) which represents one such E2LSH hash as part of the CKNS algorithm.
//
// Each unit corresponds to a given 'guess' of the value of a query KDE value,
// and a distance level from the query point. The KDE value guess is referred to
// as mu in the analysis of the CKNS algorithm, and the level is referred to by
// the index j.
//------------------------------------------------------------------------------
stag::CKNSGaussianKDEHashUnit::CKNSGaussianKDEHashUnit(
    StagReal kern_param, DenseMat* data, StagInt lognmu, StagInt j_small,
    StagReal K2_constant) {
  StagInt n = data->rows();
  StagInt d = data->cols();
  a = kern_param;
  log_nmu = lognmu;
  j = j_small;

  // Get the J parameter from the value of n and log(n * mu)
  StagInt J = ckns_J(n, log_nmu);
  assert(j <= J);

  // Initialise the random number generator for the sampling of the
  // dataset
  std::uniform_real_distribution<StagReal> uniform_distribution(0.0, 1.1);

  // Create an array of PPointT structures which will be used to point to the
  // Eigen data matrix.
  std::vector<stag::DataPoint> lsh_data;

  // We need to sample the data set with the correct sampling probability.
  StagReal p_sampling = ckns_p_sampling(j, log_nmu);

  StagUInt num_sampled_points = 0;
  for (StagInt i = 0; i < n; i++) {
    // Sample with probability p_sampling
    if (uniform_distribution(*get_global_rng()) <= p_sampling) {
      lsh_data.emplace_back(d, data->row(i).data());
      num_sampled_points++;
    }
  }

  // If the number of sampled points is below the cutoff, don't create an LSH
  // table, just store the points and we'll search through them at query time.
  if (num_sampled_points <= HASH_UNIT_CUTOFF) {
    below_cutoff = true;
    all_data = lsh_data;
  } else {
    below_cutoff = false;

    // Create the LSH parameters
    std::vector<StagUInt> lsh_parameters = ckns_gaussian_create_lsh_params(
        J, j, a, K2_constant);

    // Construct the LSH object
    LSH_buckets = stag::E2LSH(lsh_parameters[0],
                              lsh_parameters[1],
                              lsh_data);
  }
}

StagReal stag::CKNSGaussianKDEHashUnit::query(const stag::DataPoint& q) {
  std::vector<stag::DataPoint> near_neighbours;
  if (below_cutoff) {
    near_neighbours = all_data;
  } else {
    // Recover points within L_j - their distance is less than r_j.
    near_neighbours = LSH_buckets.get_near_neighbors(q);
  }

  StagReal p_sampling = ckns_p_sampling(j, log_nmu);
  StagReal rj_squared = ckns_gaussian_rj_squared(j, a);
  StagReal rj_minus_1_squared = 0;
  if (j > 1) {
    rj_minus_1_squared = ckns_gaussian_rj_squared(j-1, a);
  }

  StagReal total = 0;

  for (const auto& neighbor : near_neighbours) {
    StagReal d_sq = squared_distance(q, neighbor);

    // We return only points that are in L_j - that is, in the annulus between
    // r_{j-1} and r_j.
    if (d_sq <= rj_squared && d_sq > rj_minus_1_squared) {
      // Include this point in the estimate
      total += gaussian_kernel(a, d_sq) / p_sampling;
    }
  }
  return total;
}

//------------------------------------------------------------------------------
// CKNS Gaussian KDE
//
// We now come to the implementation of the full CKNS KDE data structure.
//------------------------------------------------------------------------------
void stag::CKNSGaussianKDE::initialize(DenseMat* data,
                                       StagReal gaussian_param,
                                       StagReal min_mu,
                                       StagInt K1,
                                       StagReal K2_constant) {
#ifndef NDEBUG
  std::cout << "Warning: STAG in debug mode!" << std::endl;
#endif
  n = data->rows();
  a = gaussian_param;

  // We are going to create a grid of LSH data structures:
  //   log2(n * mu) ranges from 0 to floor(log2(n))
  //   i ranges from 1 to k1.
  //   j ranges from 1 to J.
  min_log_nmu = (StagInt) MAX(0, floor(log2((StagReal) n * min_mu)));
  max_log_nmu = (StagInt) ceil(log2((StagReal) n));
  num_log_nmu_iterations = ceil((StagReal) (max_log_nmu - min_log_nmu) / 2);

  k1 = K1;
  k2_constant = K2_constant;
#ifndef NDEBUG
  std::cout << "[STAG] k1: " << k1 << std::endl;
  std::cout << "[STAG] k2_constant: " << k2_constant << std::endl;
#endif

  hash_units.resize(num_log_nmu_iterations);
  for (StagInt log_nmu_iter = 0;
       log_nmu_iter < num_log_nmu_iterations;
       log_nmu_iter++){
    hash_units[log_nmu_iter].resize(k1);
  }

  // For each value of n * mu, we'll create an array of LSH data structures.
  StagInt num_threads = std::thread::hardware_concurrency();
  ctpl::thread_pool pool((int) num_threads);
  StagInt max_J = ckns_J(n, 0);
  std::vector<std::future<StagInt>> futures;
  for (StagInt j_offset = max_J - 1; j_offset >= 0; j_offset--) {
    for (StagInt log_nmu_iter = 0;
         log_nmu_iter < num_log_nmu_iterations;
         log_nmu_iter++) {
      StagInt log_nmu = min_log_nmu + (log_nmu_iter * 2);
      StagInt J = ckns_J(n, log_nmu);
      StagInt j = J - j_offset;
      if (j >= 1) {
        for (StagInt iter = 0; iter < k1; iter++) {
          futures.push_back(
              pool.push(
                  [&, log_nmu_iter, log_nmu, iter, j](int id) {
                    ignore_warning(id);
                    return add_hash_unit(log_nmu_iter, log_nmu, iter, j, data);
                  }
              )
          );
        }
      }
    }
  }

  // Join all the threads
  for (auto & future : futures) {
    future.get();
  }

  // Close the thread pool.
  pool.stop();
}

stag::CKNSGaussianKDE::CKNSGaussianKDE(DenseMat *data,
                                       StagReal a,
                                       StagReal eps,
                                       StagReal min_mu) {
  n = data->rows();
  StagInt K1 = ceil(K1_DEFAULT_CONSTANT * log((StagReal) n) / SQR(eps));
  StagReal K2_constant = K2_DEFAULT_CONSTANT * log((StagReal) n);
  initialize(data, a, min_mu, K1, K2_constant);
}

stag::CKNSGaussianKDE::CKNSGaussianKDE(DenseMat *data, StagReal a) {
  n = data->rows();
  StagInt K1 = ceil(K1_DEFAULT_CONSTANT * log((StagReal) n) / SQR(EPS_DEFAULT));
  StagReal K2_constant = K2_DEFAULT_CONSTANT * log((StagReal) n);
  StagReal min_mu = 1.0 / (StagReal) n;
  initialize(data, a, min_mu, K1, K2_constant);
}

stag::CKNSGaussianKDE::CKNSGaussianKDE(
    DenseMat *data, StagReal a, StagReal eps) {
  n = data->rows();
  StagInt K1 = ceil(K1_DEFAULT_CONSTANT * log((StagReal) n) / SQR(eps));
  StagReal K2_constant = K2_DEFAULT_CONSTANT * log((StagReal) n);
  StagReal min_mu = 1.0 / (StagReal) n;
  initialize(data, a, min_mu, K1, K2_constant);
}

stag::CKNSGaussianKDE::CKNSGaussianKDE(DenseMat *data,
                                       StagReal a,
                                       StagReal min_mu,
                                       StagInt K1,
                                       StagReal K2_constant) {
  initialize(data, a, min_mu, K1, K2_constant);
}


StagInt stag::CKNSGaussianKDE::add_hash_unit(StagInt log_nmu_iter,
                                             StagInt log_nmu,
                                             StagInt iter,
                                             StagInt j,
                                             DenseMat* data) {
  assert(log_nmu < max_log_nmu);
  assert(log_nmu >= min_log_nmu);
  CKNSGaussianKDEHashUnit new_hash_unit = CKNSGaussianKDEHashUnit(
      a, data, log_nmu, j, k2_constant);
  hash_units_mutex.lock();
  hash_units[log_nmu_iter][iter].push_back(new_hash_unit);
  hash_units_mutex.unlock();
  return 0;
}

/**
 * Compute the median value of a vector.
 */
StagReal median(std::vector<StagReal> &v)
{
  size_t n = v.size() / 2;
  std::nth_element(v.begin(), v.begin()+n, v.end());
  return v[n];
}

std::vector<StagReal> stag::CKNSGaussianKDE::query(DenseMat* query_mat) {
  std::vector<StagReal> results(query_mat->rows());

  std::vector<stag::DataPoint> query_points = matrix_to_datapoints(
      query_mat);

  StagInt num_threads = std::thread::hardware_concurrency();

  // Split the query into num_threads chunks.
  if (query_mat->rows() < num_threads) {
    for (auto i = 0; i < query_mat->rows(); i++) {
      results[i] = this->query(query_points.at(i));
    }
  } else {
    // Start the thread pool
    ctpl::thread_pool pool((int) num_threads);

    StagInt chunk_size = floor(query_mat->rows() / num_threads);

    // The query size is large enough to be worth splitting.
    std::vector<std::future<std::vector<StagReal>>> futures;
    for (auto chunk_id = 0; chunk_id < num_threads; chunk_id++) {
      futures.push_back(
          pool.push(
              [&, chunk_size, chunk_id, num_threads, query_points] (int id) {
                ignore_warning(id);
                assert(chunk_id < num_threads);
                StagInt this_chunk_start = chunk_id * chunk_size;
                StagInt this_chunk_end = this_chunk_start + chunk_size;
                if (chunk_id == num_threads - 1) {
                  this_chunk_end = query_mat->rows();
                }

                assert(this_chunk_start <= (StagInt) query_points.size());
                assert(this_chunk_end <= (StagInt) query_points.size());
                assert(this_chunk_end >= this_chunk_start);

                std::vector<StagReal> chunk_results(this_chunk_end - this_chunk_start);

                for (auto i = this_chunk_start; i < this_chunk_end; i++) {
                  chunk_results[i - this_chunk_start] = this->query(query_points.at(i));
                }

                return chunk_results;
              }
          )
      );
    }

    StagInt next_index = 0;
    assert((StagInt) futures.size() == num_threads);
    for (auto chunk_id = 0; chunk_id < num_threads; chunk_id++) {
      std::vector<StagReal> chunk_results = futures[chunk_id].get();

      for (auto res : chunk_results) {
        results.at(next_index) = res;
        next_index++;
      }
    }

    pool.stop();
  }

  return results;
}

StagReal stag::CKNSGaussianKDE::query(const stag::DataPoint &q) {
  // Iterate through possible values of mu , until we find a correct one for
  // the query.
  for (auto log_nmu_iter = num_log_nmu_iterations - 1;
       log_nmu_iter >= 0;
       log_nmu_iter--) {
    StagInt log_nmu = min_log_nmu + (log_nmu_iter * 2);
    StagInt J = ckns_J(n, log_nmu);
    StagReal this_mu_estimate;

    // Get an estimate from k1 copies of the CKNS data structure.
    // Take the median one to be the true estimate.
    std::vector<StagReal> iter_estimates(k1, 0);
    for (auto iter = 0; iter < k1; iter++) {
      // Iterate through the shells for each value of j
      for (auto j = 1; j <= J; j++) {
        iter_estimates[iter] += hash_units[log_nmu_iter][iter][j - 1].query(q);
      }
    }

    this_mu_estimate = median(iter_estimates);

    // Check whether the estimate is at least mu, in which case we
    // return it.
    if (log(this_mu_estimate) >= (StagReal) log_nmu) {
      return this_mu_estimate / (StagReal) n;
    }
  }

  // Didn't find a good answer, return mu.
  return exp((StagReal) min_log_nmu) / (StagReal) n;
}

//------------------------------------------------------------------------------
// Exact Gaussian KDE Implementation
//------------------------------------------------------------------------------
StagReal gaussian_kde_exact(StagReal a,
                            const std::vector<stag::DataPoint>& data,
                            const stag::DataPoint& query) {
  StagReal total = 0;
  for (const auto& i : data) {
    total += gaussian_kernel(a, query, i);
  }
  return total / (double) data.size();
}

stag::ExactGaussianKDE::ExactGaussianKDE(DenseMat *data, StagReal param) {
  all_data = stag::matrix_to_datapoints(data);
  a = param;
}

StagReal stag::ExactGaussianKDE::query(const stag::DataPoint& q) {
  return gaussian_kde_exact(a, all_data, q);
}

std::vector<StagReal> stag::ExactGaussianKDE::query(DenseMat* query_mat) {
  std::vector<StagReal> results(query_mat->rows());

  std::vector<stag::DataPoint> query_points = stag::matrix_to_datapoints(
      query_mat);

  StagInt num_threads = std::thread::hardware_concurrency();

  // Split the query into num_threads chunks.
  if (query_mat->rows() < num_threads) {
    for (auto i = 0; i < query_mat->rows(); i++) {
      results[i] = this->query(query_points.at(i));
    }
  } else {
    // Start the thread pool
    ctpl::thread_pool pool((int) num_threads);

    StagInt chunk_size = floor(query_mat->rows() / num_threads);

    // The query size is large enough to be worth splitting.
    std::vector<std::future<std::vector<StagReal>>> futures;
    for (auto chunk_id = 0; chunk_id < num_threads; chunk_id++) {
      futures.push_back(
          pool.push(
              [&, chunk_size, chunk_id, num_threads, query_points] (int id) {
                ignore_warning(id);
                assert(chunk_id < num_threads);
                StagInt this_chunk_start = chunk_id * chunk_size;
                StagInt this_chunk_end = this_chunk_start + chunk_size;
                if (chunk_id == num_threads - 1) {
                  this_chunk_end = query_mat->rows();
                }

                assert(this_chunk_start <= (StagInt) query_points.size());
                assert(this_chunk_end <= (StagInt) query_points.size());
                assert(this_chunk_end >= this_chunk_start);

                std::vector<StagReal> chunk_results(this_chunk_end - this_chunk_start);

                for (auto i = this_chunk_start; i < this_chunk_end; i++) {
                  chunk_results[i - this_chunk_start] = gaussian_kde_exact(
                      a, all_data, query_points[i]);
                }

                return chunk_results;
              }
          )
      );
    }

    StagInt next_index = 0;
    assert((StagInt) futures.size() == num_threads);
    for (auto chunk_id = 0; chunk_id < num_threads; chunk_id++) {
      std::vector<StagReal> chunk_results = futures[chunk_id].get();

      for (auto res : chunk_results) {
        results.at(next_index) = res;
        next_index++;
      }
    }

    pool.stop();
  }

  return results;

}
