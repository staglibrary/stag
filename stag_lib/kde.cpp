/*
   This file is provided as part of the STAG library and released under the GPL
   license.
*/

#include "kde.h"

#include "random.h"
#include "multithreading/ctpl_stl.h"
#include "data.h"
#include <iostream>
#include <random>
#include <unordered_set>
#include <set>

#define TWO_ROOT_TWO 2.828427124
#define TWO_ROOT_TWOPI 5.0132565
#define LOG_TWO 0.69314718056

// The CKNS algorithm has a couple of 'magic constants' which control the
// probability guarantee and variance bounds.
#define K2_DEFAULT_CONSTANT 0.1     // K_2 = C log(n) p^{-k_j}
#define K1_DEFAULT_CONSTANT 1       // K_1 = C log(n) / eps^2
#define EPS_DEFAULT 1               // K_1 = C log(n) / eps^2
#define CKNS_DEFAULT_OFFSET 0

// At a certain number of sampled points, we might as well brute-force the hash
// unit.
#define HASH_UNIT_CUTOFF 3000

#ifndef NDEBUG
#  define LOG_DEBUG(x) do { std::cerr << x; } while (0)
#else
#  define LOG_DEBUG(x)
#endif

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
  StagReal diff;
  for (StagUInt i = 0; i < u.dimension; i++) {
    diff = u.coordinates[i] - v.coordinates[i];
    result += diff * diff;
  }
  return result;
}

/**
 * Compute the squared distance between the given two data points, if it is at
 * most max_dist. Otherwise, return -1.
 *
 * @param u the first data point
 * @param v the second data point
 * @param max_dist the maximum dist allowed.
 * @return the squared distance between \f$u\f$ and \f$v\f$.
 */
StagReal squared_distance_at_most(const stag::DataPoint& u,
                                  const stag::DataPoint& v,
                                  StagReal max_dist) {
  assert(u.dimension == v.dimension);
  StagReal result = 0;
  StagReal diff;
  for (StagUInt i = 0; i < u.dimension; i++) {
    diff = u.coordinates[i] - v.coordinates[i];
    result += diff * diff;
    if (result > max_dist) return -1;
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
  assert((StagReal) log_nmu <= ceil(log2((StagReal) n)));
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
StagReal ckns_p_sampling(StagInt j, StagInt log_nmu, StagInt n, StagInt sampling_offset) {
  // j = 0 is the special random sampling hash unit. Sample with probability
  // 2^-offset / n.
  if (j == 0) return MIN((StagReal) 1, pow(2, -sampling_offset) / n);
  else return MIN((StagReal) 1, pow(2, (StagReal) -(j+sampling_offset)) * pow(2, (StagReal) -log_nmu));
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
    StagReal K2_constant, StagInt prob_offset) {
  n = data->rows();
  StagInt d = data->cols();
  a = kern_param;
  log_nmu = lognmu;
  j = j_small;
  sampling_offset = prob_offset;

  // Get the J parameter from the value of n and log(n * mu)
  J = ckns_J(n, log_nmu);
  assert(j <= J);

  // Create an array of DataPoint structures which will be used to point to the
  // Eigen data matrix.
  std::vector<stag::DataPoint> lsh_data;

  // We need to sample the data set with the correct sampling probability.
  StagReal p_sampling = ckns_p_sampling(j, log_nmu, n, sampling_offset);

  // Get the number of points to sample
  auto num_sampled_points = (StagInt) floor(p_sampling * n);

  // And the starting index
  StagInt starting_idx = 0;
  if (p_sampling <= (StagReal) 1/2) starting_idx = (StagInt) (1 - 2 * p_sampling) * n;
  else starting_idx = 0;
  assert(starting_idx + num_sampled_points <= data->rows());
  assert(starting_idx >= 0);

  for (StagInt i = starting_idx; i < starting_idx + num_sampled_points; i++) {
    assert(i <= data->rows());
    lsh_data.emplace_back(d, data->row(i).data());
    assert(lsh_data.back().coordinates >= data->data());
  }

  // If the number of sampled points is below the cutoff, or this is the 'outer
  // ring' don't create an LSH table, just store the points and we'll search
  // through them at query time.
  if (j == 0) {
    final_shell = true;
    below_cutoff = false;
    all_data = lsh_data;
  } else if (lsh_data.size() <= HASH_UNIT_CUTOFF) {
    final_shell = false;
    below_cutoff = true;
    all_data = lsh_data;
  } else {
    below_cutoff = false;
    final_shell = false;

    // Create the LSH parameters
    std::vector<StagUInt> lsh_parameters = ckns_gaussian_create_lsh_params(
        J, j, a, K2_constant);

    // Construct the LSH object
    LSH_buckets = stag::E2LSH(lsh_parameters[0],
                              lsh_parameters[1],
                              lsh_data);
  }
}

StagReal stag::CKNSGaussianKDEHashUnit::query_neighbors(const stag::DataPoint& q,
                                                         const std::vector<stag::DataPoint>& neighbors) {
  StagReal p_sampling = ckns_p_sampling(j, log_nmu, n, sampling_offset);
  StagReal rj_squared = ckns_gaussian_rj_squared(j, a);
  StagReal rj_minus_1_squared = 0;
  if (j > 1) {
    rj_minus_1_squared = ckns_gaussian_rj_squared(j-1, a);
  } else if (j == 0) {
    // If j = 0, then this is the 'outer' shell
    rj_minus_1_squared = ckns_gaussian_rj_squared(J, a);
  }

  // Separate computation if this is the final shell - no maximum distance.
  if (final_shell) {
    StagReal total = 0;
    for (const auto& neighbor : neighbors) {
      // We use all the points in the final shell.
      StagReal d_sq = squared_distance(q, neighbor);
      if (d_sq > rj_minus_1_squared) {
        // Include this point in the estimate
        total += gaussian_kernel(a, d_sq) / p_sampling;
      }
    }
    return total;
  } else {
    StagReal total = 0;
    for (const auto& neighbor : neighbors) {
      // We return only points that are in L_j - that is, in the annulus between
      // r_{j-1} and r_j.
      StagReal d_sq = squared_distance_at_most(q, neighbor, rj_squared);
      if (d_sq > rj_minus_1_squared) {
        // Include this point in the estimate
        total += gaussian_kernel(a, d_sq) / p_sampling;
      }
    }
    return total;
  }
}

StagReal stag::CKNSGaussianKDEHashUnit::query(const stag::DataPoint& q) {
  if (below_cutoff || final_shell) {
    return query_neighbors(q, all_data);
  } else {
    std::vector<stag::DataPoint> near_neighbours = LSH_buckets.get_near_neighbors(q);
    return query_neighbors(q, near_neighbours);
  }
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
                                       StagReal K2_constant,
                                       StagInt prob_offset,
                                       StagInt min_idx,
                                       StagInt max_idx) {
#ifndef NDEBUG
  std::cerr << "Warning: STAG compiled in debug mode. For optimal performance, compile with -DCMAKE_BUILD_TYPE=Release." << std::endl;
#endif
  assert(max_idx <= data->rows());
  assert(min_idx >= 0);
  assert(min_idx < max_idx);
  min_id = min_idx;
  max_id = max_idx;
  n = max_id - min_id;
  a = gaussian_param;
  sampling_offset = prob_offset;

  // We are going to create a grid of LSH data structures:
  //   log2(n * mu) ranges from 0 to floor(log2(n))
  //   i ranges from 1 to k1.
  //   j ranges from 1 to J.
  assert(min_mu <= 1);
  min_log_nmu = (StagInt) floor(log2((StagReal) n * min_mu));
  max_log_nmu = (StagInt) ceil(log2((StagReal) n));
  assert(min_log_nmu <= max_log_nmu);

  num_log_nmu_iterations = (StagInt) ceil((StagReal) (max_log_nmu - min_log_nmu) / 2) + 1;

  k1 = K1;
  k2_constant = K2_constant;

  hash_units.resize(num_log_nmu_iterations);
  for (StagInt log_nmu_iter = 0;
       log_nmu_iter < num_log_nmu_iterations;
       log_nmu_iter++){
    hash_units[log_nmu_iter].resize(k1);
  }

  // Create K1 permuted copies of the data submatrix we are interested in.
  data_copies.reserve(K1);
  for (StagInt iter = 0; iter < k1; iter++) {
    Eigen::PermutationMatrix<Eigen::Dynamic, Eigen::Dynamic> perm(n);
    perm.setIdentity();
    std::shuffle(perm.indices().data(), perm.indices().data() + perm.indices().size(), std::mt19937(std::random_device()()));
    DenseMat permutedMatrix = perm * data->block(min_idx, 0, n, data->cols());
    data_copies.push_back(permutedMatrix);
  }

  // For each value of n * mu, we'll create an array of LSH data structures.
  StagInt num_threads = std::thread::hardware_concurrency();
  ctpl::thread_pool pool((int) num_threads);
  std::vector<std::future<StagInt>> futures;
  std::mutex hash_units_mutex;
  for (StagInt iter = 0; iter < k1; iter++) {
    for (StagInt log_nmu_iter = 0;
         log_nmu_iter < num_log_nmu_iterations;
         log_nmu_iter++) {
      StagInt log_nmu = max_log_nmu - (log_nmu_iter * 2);
      StagInt J = ckns_J(n, log_nmu);
      assert(J >= 0);

      // Make sure everything works like we expect.
      assert(log_nmu <= max_log_nmu);
      assert(log_nmu >= min_log_nmu - 1);

      // j = 0 is the special random sampling hash unit.
      for (StagInt j = 0; j <= J; j++) {
        futures.push_back(
            pool.push(
                [&, log_nmu_iter, log_nmu, iter, j](int id) {
                  ignore_warning(id);
                  return add_hash_unit(log_nmu_iter, log_nmu, iter, j, &data_copies[iter], hash_units_mutex);
                }
            )
        );
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
  initialize(data, a, min_mu, K1, K2_constant, CKNS_DEFAULT_OFFSET, 0, n);
}

stag::CKNSGaussianKDE::CKNSGaussianKDE(DenseMat *data, StagReal a) {
  n = data->rows();
  StagInt K1 = ceil(K1_DEFAULT_CONSTANT * log((StagReal) n) / SQR(EPS_DEFAULT));
  StagReal K2_constant = K2_DEFAULT_CONSTANT * log((StagReal) n);
  StagReal min_mu = 1.0 / (StagReal) n;
  initialize(data, a, min_mu, K1, K2_constant, CKNS_DEFAULT_OFFSET, 0, n);
}

stag::CKNSGaussianKDE::CKNSGaussianKDE(
    DenseMat *data, StagReal a, StagReal eps) {
  n = data->rows();
  StagInt K1 = ceil(K1_DEFAULT_CONSTANT * log((StagReal) n) / SQR(eps));
  StagReal K2_constant = K2_DEFAULT_CONSTANT * log((StagReal) n);
  StagReal min_mu = 1.0 / (StagReal) n;
  initialize(data, a, min_mu, K1, K2_constant, CKNS_DEFAULT_OFFSET, 0, n);
}

stag::CKNSGaussianKDE::CKNSGaussianKDE(DenseMat *data,
                                       StagReal a,
                                       StagReal min_mu,
                                       StagInt K1,
                                       StagReal K2_constant,
                                       StagInt prob_offset) {
  initialize(data, a, min_mu, K1, K2_constant, prob_offset, 0, data->rows());
}

stag::CKNSGaussianKDE::CKNSGaussianKDE(DenseMat *data,
                                       StagReal a,
                                       StagReal min_mu,
                                       StagInt K1,
                                       StagReal K2_constant,
                                       StagInt prob_offset,
                                       StagInt min_idx,
                                       StagInt max_idx) {
  initialize(data, a, min_mu, K1, K2_constant, prob_offset, min_idx, max_idx);
}

StagInt stag::CKNSGaussianKDE::add_hash_unit(StagInt log_nmu_iter,
                                             StagInt log_nmu,
                                             StagInt iter,
                                             StagInt j,
                                             DenseMat* data,
                                             std::mutex& hash_units_mutex) {
  assert(log_nmu <= max_log_nmu);
  assert(log_nmu >= min_log_nmu - 1);
  CKNSGaussianKDEHashUnit new_hash_unit = CKNSGaussianKDEHashUnit(
      a, data, log_nmu, j, k2_constant, sampling_offset);
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
  StagInt num_threads = std::thread::hardware_concurrency();

  // Split the query into num_threads chunks.
  if (query_mat->rows() < num_threads) {
    return this->chunk_query(query_mat, 0, query_mat->rows());
  } else {
    // Initialise the results vector
    std::vector<StagReal> results(query_mat->rows());

    // Start the thread pool
    ctpl::thread_pool pool((int) num_threads);

    StagInt chunk_size = floor((StagReal) query_mat->rows() / (StagReal) num_threads);

    // The query size is large enough to be worth splitting.
    std::vector<std::future<std::vector<StagReal>>> futures;
    for (auto chunk_id = 0; chunk_id < num_threads; chunk_id++) {
      futures.push_back(
          pool.push(
              [&, chunk_size, chunk_id, num_threads, query_mat] (int id) {
                ignore_warning(id);
                assert(chunk_id < num_threads);
                StagInt this_chunk_start = chunk_id * chunk_size;
                StagInt this_chunk_end = this_chunk_start + chunk_size;
                if (chunk_id == num_threads - 1) {
                  this_chunk_end = query_mat->rows();
                }

                assert(this_chunk_start <= (StagInt) query_mat->rows());
                assert(this_chunk_end <= (StagInt) query_mat->rows());
                assert(this_chunk_end >= this_chunk_start);

                return this->chunk_query(query_mat, this_chunk_start, this_chunk_end);
              }
          )
      );
    }

    StagInt next_index = 0;
    assert((StagInt) futures.size() == num_threads);
    for (auto chunk_id = 0; chunk_id < num_threads; chunk_id++) {
      std::vector<StagReal> chunk_results = futures[chunk_id].get();

      for (auto res : chunk_results) {
        results[next_index] = res;
        next_index++;
      }
    }

    pool.stop();

    return results;
  }
}

std::vector<StagReal> stag::CKNSGaussianKDE::chunk_query(
    DenseMat* query, StagInt chunk_start, StagInt chunk_end) {
  // Iterate through possible values of mu , until we find a correct one for
  // each query point.
  std::vector<StagReal> results(chunk_end - chunk_start, 0);
  std::unordered_set<StagInt> unsolved_queries;
  std::unordered_map<StagInt, StagReal> last_mu_estimates;
  assert(chunk_start < chunk_end);
  for (StagInt i = chunk_start; i < chunk_end; i++) {
    assert(i < chunk_end);
    unsolved_queries.insert(i);
    last_mu_estimates.insert(std::pair<StagInt, StagReal>(i, 0));
  }
  assert((StagInt) last_mu_estimates.size() == (chunk_end - chunk_start));

  for (StagInt log_nmu_iter = 0;
       log_nmu_iter < num_log_nmu_iterations;
       log_nmu_iter++) {
    StagInt log_nmu = max_log_nmu - (log_nmu_iter * 2);
    StagInt J = ckns_J(n, log_nmu);

    // Get an estimate from k1 copies of the CKNS data structure.
    // Take the median one to be the true estimate.
    std::unordered_map<StagInt, std::vector<StagReal>> iter_estimates;
    for (StagInt i : unsolved_queries) {
      assert(i < chunk_end);
      std::vector<StagReal> temp(k1, 0);
      iter_estimates.insert(std::pair<StagInt, std::vector<StagReal>>(i, temp));
    }

    for (auto iter = 0; iter < k1; iter++) {
      // Iterate through the shells for each value of j
      // Recall that j = 0 is the special random sampling unit.
      for (auto j = 0; j <= J; j++) {
        for (auto i : unsolved_queries) {
          assert(i < chunk_end);
          iter_estimates[i][iter] += hash_units[log_nmu_iter][iter][j].query(stag::DataPoint(*query, i));
        }
      }
    }

    std::vector<StagInt> newly_solved;
    for (StagInt i : unsolved_queries) {
      assert(i < chunk_end);
      StagReal this_mu_estimate = median(iter_estimates[i]);

      // Check whether the estimate is at least mu, in which case we
      // return it.
      if (log(this_mu_estimate) >= (StagReal) 1.3 * log_nmu) {
        results[i - chunk_start] = this_mu_estimate / (StagReal) n;
        newly_solved.push_back(i);
      }

      assert(i < chunk_end);
      auto test = last_mu_estimates.end();
      ignore_warning(test);
      assert(last_mu_estimates.find(i) != last_mu_estimates.end());
      last_mu_estimates[i] = this_mu_estimate;
    }

    for (StagInt i : newly_solved) {
      assert(i < chunk_end);
      auto test = unsolved_queries.end();
      ignore_warning(test);
      assert(unsolved_queries.find(i) != unsolved_queries.end());
      unsolved_queries.erase(i);
    }
  }

  // Didn't find a good answer, return the last estimate, or 0.
  for (auto i : unsolved_queries) {
    assert(i < chunk_end);
    results[i - chunk_start] = last_mu_estimates[i] / n;
  }
  return results;
}

StagReal stag::CKNSGaussianKDE::query(const stag::DataPoint &q) {
  // Iterate through possible values of mu , until we find a correct one for
  // the query.
  StagReal last_mu_estimate = 0;

  for (auto log_nmu_iter = num_log_nmu_iterations - 1;
       log_nmu_iter >= 0;
       log_nmu_iter--) {
    StagInt log_nmu = max_log_nmu - (log_nmu_iter * 2);
    StagInt J = ckns_J(n, log_nmu);
    StagReal this_mu_estimate;

    // Get an estimate from k1 copies of the CKNS data structure.
    // Take the median one to be the true estimate.
    std::vector<StagReal> iter_estimates(k1, 0);
    for (auto iter = 0; iter < k1; iter++) {
      // Iterate through the shells for each value of j
      // Recall that j = 0 is the special random sampling hash unit.
      for (auto j = 0; j <= J; j++) {
        iter_estimates[iter] += hash_units[log_nmu_iter][iter][j].query(q);
      }
    }

    this_mu_estimate = median(iter_estimates);

    // Check whether the estimate is at least mu, in which case we
    // return it.
    if (log(this_mu_estimate) >= (StagReal) 1.3 * log_nmu) {
      return this_mu_estimate / (StagReal) n;
    }

    last_mu_estimate = this_mu_estimate;
  }

  // Didn't find a good answer, return the last estimate, or 0.
  return last_mu_estimate / n;
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
  return total / (StagReal) data.size();
}

stag::ExactGaussianKDE::ExactGaussianKDE(DenseMat *data, StagReal param) {
  min_id = 0;
  max_id = data->rows();
  all_data = stag::matrix_to_datapoints(data);
  a = param;
}

stag::ExactGaussianKDE::ExactGaussianKDE(DenseMat *data, StagReal param,
                                         StagInt min_idx, StagInt max_idx) {
  min_id = min_idx;
  max_id = max_idx;
  a = param;

  for (auto i = min_id; i < max_id; i++) {
    all_data.emplace_back(data->cols(), data->row(i).data());
  }
}

StagReal stag::ExactGaussianKDE::query(const stag::DataPoint& q) {
  return gaussian_kde_exact(a, all_data, q);
}

std::vector<StagInt> stag::ExactGaussianKDE::sample_neighbors(const stag::DataPoint &q, StagReal degree, std::vector<StagReal> rs) {
  std::vector<StagInt> samples;

  std::deque<StagReal> targets;
  for (auto r : rs) {
    targets.push_back(degree * r);
  }
  std::sort(targets.begin(), targets.end());

  StagReal total = 0;
  for (StagInt i = min_id; i < max_id; i++) {
    total += gaussian_kernel(a, q, all_data[i - min_id]);

    // Get an iterator to the first element more than the total
    auto it = std::lower_bound(targets.begin(), targets.end(), total);

    // Count the elements less than the total
    StagInt count = std::distance(targets.begin(), it);

    // Remove the satisfied targets
    targets.erase(targets.begin(), it);

    // Add the samples
    for (auto j = 0; j < count; j++) {
      samples.push_back(i);
    }

    if (targets.empty()) break;
  }

  assert((StagInt) all_data.size() == max_id - min_id);
  assert(targets.empty());
  assert(samples.size() == rs.size());

  return samples;
}

std::vector<StagReal> stag::ExactGaussianKDE::query(DenseMat* query_mat) {
  std::vector<StagReal> results(query_mat->rows());

  std::vector<stag::DataPoint> query_points = stag::matrix_to_datapoints(
      query_mat);

  StagInt num_threads = std::thread::hardware_concurrency();

  // Split the query into num_threads chunks.
  if (query_mat->rows() < num_threads) {
    for (auto i = 0; i < query_mat->rows(); i++) {
      results[i] = this->query(query_points[i]);
    }
  } else {
    // Start the thread pool
    ctpl::thread_pool pool((int) num_threads);

    StagInt chunk_size = floor((StagReal) query_mat->rows() / (StagReal) num_threads);

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
        results[next_index] = res;
        next_index++;
      }
    }

    pool.stop();
  }

  return results;

}
