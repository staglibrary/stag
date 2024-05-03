/*
   This file is provided as part of the STAG library and released under the GPL
   license.
*/
// Standard C++ libraries
#include <iostream>
#include <random>
#include <deque>
#include <future>

// Additional libraries
#include "multithreading/ctpl_stl.h"

// STAG modules
#include "data.h"
#include "utility.h"

/*
 * Used to disable compiler warning for unused variable.
 */
template<class T> void ignore_warning(const T&){}

//------------------------------------------------------------------------------
// Implementation of the DataPoint class.
//------------------------------------------------------------------------------
stag::DataPoint::DataPoint(DenseMat& all_data, StagInt row_index) {
  dimension = all_data.cols();
  coordinates = all_data.row(row_index).data();
}

stag::DataPoint::DataPoint(std::vector<StagReal>& point_vector) {
  dimension = point_vector.size();
  coordinates = &point_vector[0];
}

//------------------------------------------------------------------------------
// Converting between data formats.
//------------------------------------------------------------------------------
std::vector<stag::DataPoint> stag::matrix_to_datapoints(DenseMat* data) {
  std::vector<stag::DataPoint> res;

  StagInt n = data->rows();
  StagInt d = data->cols();
  for (auto i = 0; i < n; i++) {
    res.emplace_back(d, data->row(i).data());
  }
  return res;
}

//------------------------------------------------------------------------------
// Loading and saving matrices to file.
//------------------------------------------------------------------------------
void stag::save_matrix(DenseMat& data, std::string& filename) {
  // Attempt to open the specified file
  std::ofstream os(filename);

  // If the file could not be opened, throw an exception
  if (!os.is_open()) {
    throw std::runtime_error(std::strerror(errno));
  }

  // Iterate through the entries in the matrix, and write the data file.
  for (auto i = 0; i < data.rows(); i++) {
    for (auto j = 0; j < data.cols(); j++) {
      os << data.coeffRef(i, j) << " ";
    }
    os << std::endl;
  }

  // Close the output file stream
  os.close();
}

/**
 * Parse a single content line of a 'csv'-style file containing numeric values.
 * This method assumes that the line is not a comment, and tries to parse either
 * by splitting on commas or whitespace.
 *
 * @param line the string to parse
 * @return a vector containing the values
 * @throw std::invalid_argument the line cannot be parsed
 */
std::vector<StagReal> parse_numeric_csv_content_line(std::string line) {
  // List the possible delimiters for the elements on the line
  std::vector<std::string> delimiters{",", " ", "\t"};

  // Split the line to extract the information
  StagUInt num_tokens_found = 0;
  std::vector<StagReal> tokens;

  // Try splitting by each delimiter in turn
  size_t split_pos = 0;
  std::string token;
  for (const std::string &delimiter: delimiters) {
    // If we have not found any delimiters of the previous types yet,
    // then try this one.
    if (num_tokens_found == 0) {
      while ((split_pos = line.find(delimiter)) != std::string::npos) {
        // Extract the portion of the line up to the delimiter
        token = line.substr(0, split_pos);
        line.erase(0, split_pos + delimiter.length());

        // If the token has length 0, then skip
        if (token.length() == 0) continue;

        // Parse the token as a real
        // throws an exception if the token cannot be parsed.
        tokens.push_back(std::stod(token));

        // Increase the counter of the number of tokens found
        num_tokens_found++;
      }
    }
  }

  // Extract the final token in the line
  try {
    // Try extracting another token from the rest of the line, but ignore any
    // errors - there might be none.
    StagReal last_token = std::stod(line);
    tokens.push_back(last_token);
    num_tokens_found++;
  } catch (std::exception &e) {
    // Ignore this
  }

  // Return the data
  return tokens;
}

DenseMat stag::load_matrix(std::string& filename) {
  // Attempt to open the provided file
  std::ifstream is(filename);

  // If the file could not be opened, throw an exception
  if (!is.is_open()) {
    throw std::runtime_error(std::strerror(errno));
  }

  // Initialise the data matrix to be 1 x 1.
  // Once we know the number of columns, we will update num_cols
  StagInt num_cols = -1;
  StagInt num_rows = 0;
  DenseMat data(1, 1);

  // Read in the file one line at a time
  std::string line;
  while (stag::safeGetline(is, line)) {
    if (line.length() > 0 && line[0] != '#' && line[0] != '/') {
      // Get all of the data from this line
      std::vector<StagReal> tokens = parse_numeric_csv_content_line(line);

      // If we don't yet know the number of columns in the data matrix, set it
      // now.
      if (num_cols < 0) {
        num_cols = (StagInt) tokens.size();
        data.conservativeResize(Eigen::NoChange_t::NoChange, num_cols);
      }

      // If the number of columns on this line does not match previous lines,
      // then throw an exception.
      if (num_cols != (StagInt) tokens.size()) {
        throw std::runtime_error("Rows of different size in matrix file.");
      }

      num_rows++;

      if (num_rows <= data.rows()) {
        for (auto i = 0; i < num_cols; i++) {
          data.coeffRef(num_rows-1, i) = tokens[i];
        }
      } else {
        // We need to resize the data matrix
        data.conservativeResize(ceil(1.1 * (StagReal) num_rows),
                                Eigen::NoChange_t::NoChange);

        for (auto i = 0; i < num_cols; i++) {
          data.coeffRef(num_rows-1, i) = tokens[i];
        }
      }
    }
  }

  // Finally, we resize the matrix to its final size.
  data.conservativeResize(num_rows, Eigen::NoChange_t::NoChange);

  return data;
}
