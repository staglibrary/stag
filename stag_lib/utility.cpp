/*
   This file is provided as part of the STAG library and released under the GPL
   license.
*/
#include <iterator>
#include <filesystem>
#include "utility.h"

std::vector<StagInt> stag::sprsMatInnerIndices(const SprsMat *matrix) {
  // Make sure that the given matrix is compressed
  assert(matrix->isCompressed());

  // Return the required indices vector
  const StagInt *indexPtr = matrix->innerIndexPtr();
  StagInt nonZeros = matrix->nonZeros();
  return {indexPtr, indexPtr + nonZeros};
}

std::vector<StagInt> stag::sprsMatOuterStarts(const SprsMat* matrix) {
  // Make sure that the given matrix is compressed
  assert(matrix->isCompressed());

  // Return the required indices vector
  const StagInt *indexPtr = matrix->outerIndexPtr();
  StagInt outerSize = matrix->outerSize();
  return {indexPtr, indexPtr + outerSize + 1};
}

std::vector<StagReal> stag::sprsMatValues(const SprsMat* matrix) {
  // Make sure that the given matrix is compressed
  assert(matrix->isCompressed());

  // Return the required indices vector
  const StagReal *valuePtr = matrix->valuePtr();
  StagInt nonZeros = matrix->nonZeros();
  return {valuePtr, valuePtr + nonZeros};
}

std::vector<StagReal> stag::sprsMatToVec(const SprsMat* matrix) {
  // If the number of dimensions is not given, use the dimension of the sparse
  // matrix.
  return stag::sprsMatToVec(matrix, matrix->rows());
}

std::vector<StagReal> stag::sprsMatToVec(const SprsMat* matrix, StagInt n) {
  if (n < 1) throw std::invalid_argument("Dimension n must be at least 1.");

  // Initialise the solution vector.
  std::vector<StagReal> dense_vec;

  for (StagInt i = 0; i < n; i++) {
    if (i < matrix->rows()) {
      // Get the i-th entry of the sparse matrix
      dense_vec.push_back(matrix->coeff(i, 0));
    } else {
      // If the sparse matrix is not long enough, fill the vector with 0s.
      dense_vec.push_back(0);
    }
  }
  return dense_vec;
}

SprsMat stag::sprsMatFromVectors(std::vector<StagInt>& column_starts,
                                 std::vector<StagInt>& row_indices,
                                 std::vector<StagReal>& values) {
  // The length of the row_indices and values vectors should be the same
  if (row_indices.size() != values.size()) {
    throw std::invalid_argument("Sparse matrix indices and values array length mismatch.");
  }

  // The last value in the column_starts vector should be equal to the length
  // of the data vectors.
  if (column_starts.back() != (StagInt) row_indices.size()) {
    throw std::invalid_argument("Final column starts entry should equal size of data vectors.");
  }

  SprsMat constructed_mat = Eigen::Map<SprsMat>((StagInt) column_starts.size() - 1,
                                                (StagInt) column_starts.size() - 1,
                                                (StagInt) values.size(),
                                                column_starts.data(),
                                                row_indices.data(),
                                                values.data());
  constructed_mat.makeCompressed();
  return constructed_mat;
}

bool stag::isSymmetric(const SprsMat *matrix) {
  // Iterate through the non-zero elements in the matrix
  for (int k = 0; k < matrix->outerSize(); ++k) {
    for (SprsMat::InnerIterator it(*matrix, k); it; ++it) {
      // If the value in the symmetrically opposite position is not the same,
      // then return false.
      if (it.value() != matrix->coeff(it.col(), it.row())) {
        return false;
      }
    }
  }

  // We didn't find any symmetrically opposite coefficients with different
  // values, and so this matrix is symmetric.
  return true;
}

std::istream& stag::safeGetline(std::istream& is, std::string& t)
{
    t.clear();

    // The characters in the stream are read one-by-one using a std::streambuf.
    // That is faster than reading them one-by-one using the std::istream.
    // Code that uses streambuf this way must be guarded by a sentry object.
    // The sentry object performs various tasks,
    // such as thread synchronization and updating the stream state.

    std::istream::sentry se(is, true);
    std::streambuf* sb = is.rdbuf();

    for(;;) {
        int c = sb->sbumpc();
        switch (c) {
            case '\n':
                return is;
            case '\r':
                if(sb->sgetc() == '\n')
                    sb->sbumpc();
                return is;
            case std::streambuf::traits_type::eof():
                // Also handle the case when the last line has no line ending
                if(t.empty())
                    is.setstate(std::ios::eofbit);
                return is;
            default:
                t += (char)c;
        }
    }
}

std::string random_string( size_t length )
{
  auto randchar = []() -> char
  {
    const char charset[] =
        "0123456789"
        "ABCDEFGHIJKLMNOPQRSTUVWXYZ"
        "abcdefghijklmnopqrstuvwxyz";
    const size_t max_index = (sizeof(charset) - 1);
    return charset[ std::rand() % max_index ];
  };
  std::string str(length,0);
  std::generate_n( str.begin(), length, randchar );
  return str;
}

std::string stag::getTempFilename() {
  // Get the name of the temporary directory on this filesystem
  std::filesystem::path temp_dir = std::filesystem::temp_directory_path();

  // Create the random filename
  std::filesystem::path fname("stag_temp_file." + random_string(20));

  // Create the full path
  std::filesystem::path full_path = temp_dir / fname;
  return full_path.string();
}

std::string stag::openTempFile(std::ofstream* os) {
  std::string temp_fname = stag::getTempFilename();

  // Open the ofstream
  *os = std::ofstream(temp_fname);

  // Return the name of the file
  return temp_fname;
}
