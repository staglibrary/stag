//
// Methods for solving Laplacian systems of equations.
//
// This file is provided as part of the STAG library and released under the MIT
// license.
//

/**
 * @file solve.h
 * \brief Methods for solving Laplacian systems of equations.
 */

#ifndef STAG_TEST_SOLVE_H
#define STAG_TEST_SOLVE_H

#include "graph.h"

namespace stag {

  DenseVec solve_laplacian(Graph* g, DenseVec& b, double eps);

  DenseVec jacobi_iteration(const SprsMat* A, DenseVec& b, double eps);
}

#endif //STAG_TEST_SOLVE_H
