//
// This file is provided as part of the STAG library and released under the MIT
// license.
//
#include <iostream>
#include <fstream>
#include <stdexcept>
#include <filesystem>

#include "graph.h"
#include "utility.h"
#include "graphio.h"

/**
 * Parse a single content line of an edgelist file. This method assumes that
 * the line is not a comment, and tries to parse either by splitting on commas
 * or whitespace.
 *
 * @return a triple representing the edge (u, v, weight).
 * @throw std::invalid_argument the line cannot be parsed
 */
stag::edge parse_edgelist_content_line(std::string line) {
  // List the possible delimiters for the elements on the line
  std::vector<std::string> delimiters{",", " ", "\t"};

  // Split the line to extract the edge information
  // The weight defaults to 1 if it is not otherwise updated by the information
  // in the line.
  int num_tokens_found = 0;
  int u = -1;
  int v = -1;
  double weight = 1;

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

        // Parse the token as the appropriate data type - int or double
        // throws an exception if the token cannot be parsed.
        if (num_tokens_found == 0) u = std::stoi(token);
        if (num_tokens_found == 1) v = std::stoi(token);
        if (num_tokens_found == 2) weight = std::stod(token);

        // Increase the counter of the number of tokens found
        num_tokens_found++;
      }
    }
  }

  // Extract the final token in the line
  if (num_tokens_found == 1) {
    v = std::stoi(line);
    num_tokens_found++;
  } else if (num_tokens_found == 2) {
    try {
      // Try extracting the weight from the rest of the line, but ignore any
      // errors - the weight might not be there.
      weight = std::stod(line);
      num_tokens_found++;
    } catch (std::exception &e) {
      // Ignore this
    }
  }

  // Check that we have exactly two or three elements in the split line
  if (num_tokens_found < 2 || num_tokens_found > 3) {
    throw std::invalid_argument("Wrong number of tokens on edgelist line.");
  }

  // Make sure that the vertices u and v have been updated successfully
  if (u < 0 || v < 0) {
    throw std::invalid_argument("Parse error on edgelist line.");
  }

  // Return the triple
  return {u, v, weight};
}

stag::Graph stag::load_edgelist(std::string &filename) {
  // Attempt to open the provided file
  std::ifstream is(filename);

  // If the file could not be opened, throw an exception
  if (!is.is_open()) {
    throw std::runtime_error(std::strerror(errno));
  }

  // We will construct a vector of triples in order to construct the final
  // adjacency matrix
  std::vector<EdgeTriplet> non_zero_entries;

  // Read the file in one line at a time
  stag_int number_of_vertices = 0;
  std::string line;
  stag::edge this_edge;
  while (stag::safeGetline(is, line)) {
    if (line[0] != '#' && line[0] != '/' && line.length() > 0) {
      try {
        // This line of the input file isn't a comment, parse it.
        this_edge = parse_edgelist_content_line(line);

        // Add two edges to the adjacency matrix in order to keep it symmetric.
        non_zero_entries.emplace_back(
            EdgeTriplet(this_edge.v1, this_edge.v2, this_edge.weight));
        non_zero_entries.emplace_back(
            EdgeTriplet(this_edge.v2, this_edge.v1, this_edge.weight));

        // Update the number of vertices to be the maximum of the column and row
        // indices.
        number_of_vertices = std::max(number_of_vertices, this_edge.v1 + 1);
        number_of_vertices = std::max(number_of_vertices, this_edge.v2 + 1);
      } catch (std::invalid_argument &e) {
        // Re-throw any parsing errors
        throw(std::runtime_error(e.what()));
      }
    }
  }

  // Close the input file stream
  is.close();

  // Update the adjacency matrix from the triples constructed from the input file.
  SprsMat adj_mat(number_of_vertices, number_of_vertices);
  adj_mat.setFromTriplets(non_zero_entries.begin(), non_zero_entries.end());

  // Construct and return the graph object
  return stag::Graph(adj_mat);
}

void stag::save_edgelist(stag::Graph &graph, std::string &filename) {
  // Attempt to open the specified file
  std::ofstream os(filename);

  // If the file could not be opened, throw an exception
  if (!os.is_open()) {
    throw std::runtime_error(std::strerror(errno));
  }

  // Write header information to the output file.
  os << "# This file was automatically generated by the STAG library." << std::endl;
  os << "#   number of vertices = " << graph.number_of_vertices() << std::endl;
  os << "#   number of edges = " << graph.number_of_edges() << std::endl;

  // Iterate through the entries in the graph adjacency matrix, and write
  // the edgelist file
  const SprsMat* adj_mat = graph.adjacency();
  for (int k = 0; k < adj_mat->outerSize(); ++k) {
    for (SprsMat::InnerIterator it(*adj_mat, k); it; ++it) {
      // We only consider the 'upper triangle' of the matrix
      if (it.col() > it.row()) {
        os << it.row() << " " << it.col() << " " << it.value() << std::endl;
      }
    }
  }

  // Close the output file stream
  os.close();
}

//------------------------------------------------------------------------------
// Sorting and manipulating edgelist on disk.
//------------------------------------------------------------------------------

void stag::copy_edgelist_duplicate_edges(std::string& infile, std::string& outfile) {
  // Open the input file
  std::ifstream is(infile);
  if (!is.is_open()) {
    throw std::runtime_error(std::strerror(errno));
  }

  // Open the output file
  std::ofstream os(outfile);
  if (!os.is_open()) {
    throw std::runtime_error(std::strerror(errno));
  }

  // Read the file in one line at a time
  std::string line;
  stag::edge this_edge;
  while (stag::safeGetline(is, line)) {
    if (line[0] != '#' && line[0] != '/' && line.length() > 0) {
      try {
        // This line of the input file isn't a comment, parse it.
        this_edge = parse_edgelist_content_line(line);

        // And write both directions to the output file
        os << this_edge.v1 << " " << this_edge.v2 << " " << this_edge.weight << std::endl;
        os << this_edge.v2 << " " << this_edge.v1 << " " << this_edge.weight << std::endl;
      } catch (std::invalid_argument &e) {
        // Re-throw any parsing errors
        throw(std::runtime_error(e.what()));
      }
    } else {
      // This line is a comment - pass it verbatim to the output
      os << line << std::endl;
    }
  }

  // Close the file streams
  is.close();
  os.close();
}

/**
 * Get the number of lines and maximum ID in an edgelist file.
 *
 * This has running time O(n) where n is the number of lines in the edgelist
 * file.
 */
void get_edgelist_lines_and_nodes(std::string& filename, stag_int& lines,
                                  stag_int& max_id) {
  // Open the input file
  std::ifstream is(filename);
  if (!is.is_open()) {
    throw std::runtime_error(std::strerror(errno));
  }

  // Find the number of lines and the maximum node ID in the edgelist file
  max_id = 0;
  lines = 0;
  std::string line;
  stag::edge this_edge;
  while (stag::safeGetline(is, line)) {
    lines++;

    if (line[0] != '#' && line[0] != '/' && line.length() > 0) {
      try {
        // This line of the input file isn't a comment, parse it.
        this_edge = parse_edgelist_content_line(line);

        // Update the number of vertices to be the maximum of the column and row
        // indices.
        max_id = std::max(max_id, this_edge.v1 + 1);
        max_id = std::max(max_id, this_edge.v2 + 1);
      } catch (std::invalid_argument &e) {
        // Re-throw any parsing errors
        throw(std::runtime_error(e.what()));
      }
    }
  }

  // Close the input file - it will be re-opened in each iteration of the
  // quicksort algorithm.
  is.close();
}

/**
 * A structure representing an interval for the quicksort algorithm for edgelist
 * files.
 */
struct EdgelistSortInterval {
  stag_int start_line;
  stag_int end_line;
  stag_int min_id;
  stag_int max_id;
};

void stag::sort_edgelist(std::string &filename) {
  // Find the number of lines and the maximum node ID in the edgelist file
  stag_int max_id = 0;
  stag_int num_lines = 0;
  get_edgelist_lines_and_nodes(filename, num_lines, max_id);

  // Initialise the vector of quicksort intervals
  std::vector<EdgelistSortInterval> intervals;
  intervals.push_back({0, num_lines - 1, 0, max_id});

  // Iterate the quicksort algorithm until there are no intervals left.
  while (!intervals.empty()) {
    // Create a temporary file for this iteration.
    std::ofstream os;
    std::string temp_fname = stag::openTempFile(&os);
    if (!os.is_open()) throw std::runtime_error(std::strerror(errno));

    // Open the edgelist file as input
    std::ifstream ifs(filename);
    if (!ifs.is_open()) throw std::runtime_error(std::strerror(errno));

    // Throughout the algorithm, we must maintain a record of which line of the
    // input and output file we are pointing at.
    stag_int current_input_line = 0;
    stag_int current_output_line = 0;

    // For each pass through the input file, output any header comment lines
    // verbatim
    bool written_content = false;

    // For each interval, we will iterate over that portion of the input file
    // twice.
    std::string line;
    stag::edge this_edge;
    std::vector<EdgelistSortInterval> new_intervals;
    for (EdgelistSortInterval interval : intervals) {
      // The pivot index is half-way between the max and min.
      stag_int double_pivot = interval.min_id + interval.max_id;
      assert(double_pivot > 2 * interval.min_id);
      assert(double_pivot < 2 * interval.max_id);

      // First iteration: looking for edges with node_ids less than the pivot
      while (current_input_line < interval.start_line) {
        // Write out every line up to the start point.
        stag::safeGetline(ifs, line);
        os << line << std::endl;
        current_output_line++;
        current_input_line++;
      }
      std::streampos start_loc = ifs.tellg();
      assert(current_input_line == interval.start_line);

      stag_int new_interval_start = current_output_line;
      stag_int new_interval_min_id = interval.max_id;
      stag_int new_interval_max_id = interval.min_id;

      while (current_input_line < interval.end_line) {
        stag::safeGetline(ifs, line);
        current_input_line++;

        if (line[0] != '#' && line[0] != '/' && line.length() > 0) {
          try {
            // This line of the input file isn't a comment, parse it.
            this_edge = parse_edgelist_content_line(line);
            written_content = true;

            if (2 * this_edge.v1 < double_pivot) {
              os << line << std::endl;
              current_output_line++;

              // Update the maxs and mins
              if (this_edge.v1 < new_interval_min_id) new_interval_min_id = this_edge.v1;
              if (this_edge.v1 > new_interval_max_id) new_interval_max_id = this_edge.v1;
            }
          } catch (std::invalid_argument &e) {
            // Re-throw any parsing errors
            throw(std::runtime_error(e.what()));
          }
        } else {
          if (!written_content) {
            os << line << std::endl;
            current_output_line++;
          }
        }
      }

      assert(new_interval_max_id >= new_interval_min_id);
      stag_int new_interval_end = current_output_line;
      if (new_interval_end - new_interval_start > 1 &&
          new_interval_max_id > new_interval_min_id) {
        new_intervals.push_back({new_interval_start,
                                 new_interval_end,
                                 new_interval_min_id,
                                 new_interval_max_id});
      }

      // Return to the start of this interval
      ifs.seekg(start_loc);
      current_input_line = interval.start_line;

      new_interval_start = current_output_line;
      new_interval_min_id = interval.max_id;
      new_interval_max_id = interval.min_id;

      while (current_input_line < interval.end_line) {
        stag::safeGetline(ifs, line);
        current_input_line++;

        if (line[0] != '#' && line[0] != '/' && line.length() > 0) {
          try {
            // This line of the input file isn't a comment, parse it.
            this_edge = parse_edgelist_content_line(line);

            if (2 * this_edge.v1 >= double_pivot) {
              os << line << std::endl;
              current_output_line++;

              // Update the maxs and mins
              if (this_edge.v1 < new_interval_min_id) new_interval_min_id = this_edge.v1;
              if (this_edge.v1 > new_interval_max_id) new_interval_max_id = this_edge.v1;
            }
          } catch (std::invalid_argument &e) {
            // Re-throw any parsing errors
            throw(std::runtime_error(e.what()));
          }
        }
      }

      assert(new_interval_max_id >= new_interval_min_id);
      new_interval_end = current_output_line;
      if (new_interval_end - new_interval_start > 1 &&
          new_interval_max_id > new_interval_min_id) {
        new_intervals.push_back({new_interval_start,
                                 new_interval_end,
                                 new_interval_min_id,
                                 new_interval_max_id});
      }

    }

    // Update the intervals
    intervals = new_intervals;

    // Close the streams.
    ifs.close();
    os.close();

    // Copy the temporary file over the original edgelist.
    std::filesystem::remove(filename);
    std::filesystem::copy(temp_fname, filename);
    std::filesystem::remove(temp_fname);
  }
}

//------------------------------------------------------------------------------
// Adjacency List processing
//------------------------------------------------------------------------------

stag::edge parse_adjacencylist_edge(std::string token, stag_int source_node) {
  stag_int neighbour;
  double weight;

  // Try to split on ':' to get a weight
  size_t split_pos = token.find(':');
  if (split_pos == std::string::npos) {
    // There is no colon, so the entire token should be parseable as an integer
    // and this is an edge with weight 1.
    neighbour = std::stoi(token);
    weight = 1;
  } else {
    // There is a colon. Split the token to get the neighbour and the weight.
    neighbour = std::stoi(token.substr(0, split_pos));
    token.erase(0, split_pos + 1);
    weight = std::stod(token);
  }

  return {source_node, neighbour, weight};
}

/**
 * Comparison function for sorting edges.
 */
bool cmp_neighbors(const stag::edge& a, const stag::edge& b) {
  return a.v2 < b.v2;
}


std::vector<stag::edge> stag::parse_adjacencylist_content_line(std::string line) {
  std::vector<stag::edge> edges;

  // Begin by finding the ID of the node at the start of the line
  size_t split_pos = line.find(':');
  if (split_pos == std::string::npos) {
    throw std::invalid_argument("Couldn't extract ID on adjacencylist line.");
  }
  std::string token = line.substr(0, split_pos);
  stag_int source_node_id = std::stoi(token);
  line.erase(0, split_pos + 1);

  // Now, repeatedly split the rest of the line on spaces to find the neighbours.
  while ((split_pos = line.find(' ')) != std::string::npos) {
    token = line.substr(0, split_pos);
    line.erase(0, split_pos + 1);
    if (token.length() == 0) continue;

    // We have found a token representing some edge. Parse it
    edges.push_back(parse_adjacencylist_edge(token, source_node_id));
  }

  // Check for another neighbour at the end of the line
  try {
    edges.push_back(parse_adjacencylist_edge(line, source_node_id));
  } catch (std::exception& e) {
    // Ignore any exceptions - there may not be a neighbour to parse.
  }

  // Sort the neighbours. This is not strictly necessary, but will ensure
  // that different methods of accessing the adjacencylist graph will behave
  // as similarly as possible.
  std::stable_sort(edges.begin(), edges.end(), cmp_neighbors);

  return edges;
}

stag::Graph stag::load_adjacencylist(std::string &filename) {
  // Attempt to open the provided file
  std::ifstream is(filename);

  // If the file could not be opened, throw an exception
  if (!is.is_open()) {
    throw std::runtime_error(std::strerror(errno));
  }

  // We will construct a vector of triples in order to construct the final
  // adjacency matrix
  std::vector<EdgeTriplet> non_zero_entries;

  // Read the file in one line at a time
  stag_int number_of_vertices = 0;
  std::string line;
  std::vector<stag::edge> neighbours;
  while (stag::safeGetline(is, line)) {
    if (line[0] != '#' && line[0] != '/' && line.length() > 0) {
      try {
        // This line of the input file isn't a comment, parse it.
        neighbours = stag::parse_adjacencylist_content_line(line);

        // Add the edges to the adjacency matrix
        for (auto this_edge : neighbours) {
          non_zero_entries.emplace_back(
              this_edge.v1, this_edge.v2, this_edge.weight);

          // Update the number of vertices to be the maximum of the column and row
          // indices.
          number_of_vertices = std::max(number_of_vertices, this_edge.v1 + 1);
          number_of_vertices = std::max(number_of_vertices, this_edge.v2 + 1);
        }
      } catch (std::invalid_argument &e) {
        // Re-throw any parsing errors
        throw(std::runtime_error(e.what()));
      }
    }
  }

  // Close the input file stream
  is.close();

  // Update the adjacency matrix from the triples constructed from the input file.
  SprsMat adj_mat(number_of_vertices, number_of_vertices);
  adj_mat.setFromTriplets(non_zero_entries.begin(), non_zero_entries.end());

  // Construct and return the graph object
  return stag::Graph(adj_mat);
}

void stag::save_adjacencylist(stag::Graph &graph, std::string &filename) {
  // Attempt to open the specified file
  std::ofstream os(filename);

  // If the file could not be opened, throw an exception
  if (!os.is_open()) {
    throw std::runtime_error(std::strerror(errno));
  }

  // Write header information to the output file.
  os << "# This file was automatically generated by the STAG library." << std::endl;
  os << "#   number of vertices = " << graph.number_of_vertices() << std::endl;
  os << "#   number of edges = " << graph.number_of_edges() << std::endl;
  os << "#" << std::endl;
  os << "# This file is formatted as a STAG adjacency list. For more information," << std::endl;
  os << "# see the STAG documentation." << std::endl;
  os << "#" << std::endl;

  // Iterate through the nodes in the graph, and write
  // the adjacencylist file
  for (stag_int node = 0; node < graph.number_of_vertices(); node++) {
    os << node << ":";

    for (stag::edge e: graph.neighbors(node)) {
      os << " " << e.v2 << ":" << e.weight;
    }

    os << std::endl;
  }

  // Close the output file stream
  os.close();
}

void stag::adjacencylist_to_edgelist(std::string &adjacencylist_fname, std::string &edgelist_fname) {
  // Open input and output streams
  std::ifstream is(adjacencylist_fname);
  std::ofstream os(edgelist_fname);
  if (!is.is_open()) throw std::runtime_error(std::strerror(errno));
  if (!os.is_open()) throw std::runtime_error(std::strerror(errno));

  // We will include any comments up until the first content line.
  // This preserves 'header' information as a comment.
  bool written_content = false;

  // Read the file in one line at a time
  std::string line;
  std::vector<stag::edge> neighbours;
  while (stag::safeGetline(is, line)) {
    if (line[0] != '#' && line[0] != '/' && line.length() > 0) {
      try {
        // This line of the input file isn't a comment, parse it.
        neighbours = stag::parse_adjacencylist_content_line(line);
        written_content = true;

        // Add the edges to the edgelist file
        for (auto this_edge : neighbours) {
          // Since the adjacency list contains edges in both direction, only
          // add them to the edgelist when v1 is less than v2
          if (this_edge.v1 <= this_edge.v2) {
            os << this_edge.v1 << " " << this_edge.v2 << " " << this_edge.weight << std::endl;
          }
        }
      } catch (std::invalid_argument &e) {
        // Re-throw any parsing errors
        throw(std::runtime_error(e.what()));
      }
    } else {
      if (!written_content) {
        os << line << std::endl;
      }
    }
  }

  // Close the file streams
  is.close();
  os.close();
}

void stag::edgelist_to_adjacencylist(std::string &edgelist_fname,
                                     std::string &adjacencylist_fname) {
  // Begin by copying and sorting the edgelist file
  std::string temp_edgelist_filename = stag::getTempFilename();
  stag::copy_edgelist_duplicate_edges(edgelist_fname, temp_edgelist_filename);
  stag::sort_edgelist(temp_edgelist_filename);

  // Open the input and output streams.
  std::ifstream is(temp_edgelist_filename);
  std::ofstream os(adjacencylist_fname);

  // We will include any comments up until the first content line.
  // This preserves 'header' information as a comment.
  bool written_content = false;

  // Iterate through the edgelist file
  stag_int current_node = -1;
  std::string line;
  while (stag::safeGetline(is, line)) {
    if (line[0] != '#' && line[0] != '/' && line.length() > 0) {
      try {
        // This line of the input file isn't a comment, parse it.
        stag::edge this_edge = parse_edgelist_content_line(line);
        written_content = true;

        // If this is a larger node than we've seen so far, begin a new
        // line of the adjacency list file.
        assert(this_edge.v1 >= current_node);
        if (this_edge.v1 > current_node) {
          os << std::endl;
          os << this_edge.v1 << ":";
          current_node = this_edge.v1;
        }

        // Add the edge to the current line
        os << " " << this_edge.v2 << ":" << this_edge.weight;
      } catch (std::invalid_argument &e) {
        // Re-throw any parsing errors
        throw (std::runtime_error(e.what()));
      }
    } else {
      if (!written_content) {
        os << line << std::endl;
      }
    }
  }

  // Close the file streams
  is.close();
  os.close();

  // Delete the temporary file
  std::filesystem::remove(temp_edgelist_filename);
}
