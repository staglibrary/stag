/**
 * Tests for the data.h header. Includes methods for reading and writing
 * data to disk.
 *
 * This file is provided as part of the STAG library and released under the GPL
 * license.
*/
#include <iostream>
#include <gtest/gtest.h>
#include <filesystem>
#include "KDE/data.h"
#include "utility.h"
#include "Graph/graph.h"
#include "Graph/random.h"
#include "Cluster/cluster.h"

// Note: for many of these tests, we re-use the edgelist and adjacencylist files
// however, in general the load_matrix and save_matrix methods should be used
// for data matrices, not graph edgelists.

#define EXPECT_MATRICES_EQ(expected, actual) \
        EXPECT_EQ(expected.rows(), actual.rows());\
        EXPECT_EQ(expected.cols(), actual.cols());\
        for (auto i = 0; i < expected.rows(); i++) { \
          for (auto j = 0; j < expected.cols(); j++) { \
            EXPECT_EQ(expected.coeff(i, j), actual.coeff(i, j)) \
              << "Mismatch on row " << i << " col " << j;\
          } \
        }


TEST(DataTest, LoadMatrixSimple) {
  std::string filename = "test/data/test1.edgelist";
  DenseMat data = stag::load_matrix(filename);

  // Create the expected data for the matrix.
  DenseMat expected_data {{0, 1}, {1, 2}, {2, 0}};

  // Check that the matrices match.
  EXPECT_MATRICES_EQ(expected_data, data);
}

TEST(DataTest, LoadMatrixDoubles) {
  std::string filename = "test/data/test2.edgelist";
  DenseMat data = stag::load_matrix(filename);

  // Create the expected data for the matrix.
  DenseMat expected_data {{0, 1, 0.5}, {1, 2, 1}, {2, 0, 0.5}};

  // Check that the matrices match.
  EXPECT_MATRICES_EQ(expected_data, data);
}

TEST(DataTest, LoadMatrixDifferentColLengths) {
  // Check that when the matrix has different column lengths, the load_matrix
  // method throws an exception.
  std::string filename = "test/data/test3.edgelist";

  EXPECT_THROW(stag::load_matrix(filename),
               std::runtime_error);
}

TEST(DataTest, LoadMatrixTabs) {
  std::string filename = "test/data/test5.edgelist";
  DenseMat data = stag::load_matrix(filename);

  // Create the expected data for the matrix.
  DenseMat expected_data {{0, 1, 0.5}, {1, 2, 1}, {2, 0, 0.5}};

  // Check that the matrices match.
  EXPECT_MATRICES_EQ(expected_data, data);
}

TEST(DataTest, MatrixIOBadFilename) {
  std::string badFilename = "thisfiledoesntexist.edgelist";
  EXPECT_THROW({DenseMat data = stag::load_matrix(badFilename);},
               std::runtime_error);

  // Now try writing to a file which cannot be written to
  DenseMat data {{1, 0}, {3, 5}, {2, 6}};
  badFilename = "////cannotwritehere.edgelist";
  EXPECT_THROW(stag::save_matrix(data, badFilename),
               std::runtime_error);
}

TEST(DataTest, SaveMatrix) {
  DenseMat data {{1, 0}, {3, 5}, {2, 6}};

  // Save a matrix file
  std::string filename = "output.txt";
  stag::save_matrix(data, filename);

  // Reading the matrix back in should result in the same matrix
  DenseMat new_data = stag::load_matrix(filename);

  EXPECT_MATRICES_EQ(data, new_data);
}

TEST(DataTest, MatrixToDataPoints) {
  // Create a matrix
  DenseMat data {{1, 0}, {3, 5}, {2, 6}};

  // Convert it to an array of data point objects
  std::vector<stag::DataPoint> datapoints = stag::matrix_to_datapoints(&data);

  // Check that each point is correct
  for (auto i = 0; i < datapoints.size(); i++) {
    for (auto j = 0; j < datapoints[i].dimension; j++) {
      EXPECT_EQ(data.coeff(i, j), datapoints[i].coordinates[j]);
    }
  }
}

TEST(DataTest, ASGMNIST) {
  // Load the MNIST dataset
  std::string filename = "test/data/mnist.txt";
  DenseMat data = stag::load_matrix(filename);
  StagReal a = 0.000001;

  // Create tha approximate similarity graph from this matrix.
  stag::Graph asg = stag::approximate_similarity_graph(&data, a);
  EXPECT_EQ(asg.number_of_vertices(), data.rows());

  // Load the correct clusters
  std::string labels_filename = "test/data/mnist_labels.txt";
  DenseMat labels = stag::load_matrix(labels_filename);
  std::vector<StagInt> labels_vec(labels.rows());
  for (auto i = 0; i < labels.rows(); i++) {
    labels_vec.at(i) = (StagInt) labels(i, 0);
  }

  // Check the clustering performance with the asg roughly matches the
  // performance with the fully connected graph
  StagInt k = 10;
  std::vector<StagInt> clusters = stag::spectral_cluster(&asg, k);
  StagReal asg_ari = stag::adjusted_rand_index(clusters, labels_vec);

  stag::Graph sg = stag::similarity_graph(&data, a);
  clusters = stag::spectral_cluster(&sg, k);
  StagReal fc_ari = stag::adjusted_rand_index(clusters, labels_vec);
  EXPECT_GE(asg_ari, 0.8 * fc_ari);
}

TEST(DataTest, ASGmoons) {
  // Load the moons dataset
  std::string filename = "test/data/moons.txt";
  DenseMat data = stag::load_matrix(filename);
  StagReal a = 10;

  // Create tha approximate similarity graph from this matrix.
  stag::Graph asg = stag::approximate_similarity_graph(&data, a);
  EXPECT_EQ(asg.number_of_vertices(), data.rows());

  // Load the correct clusters
  std::string labels_filename = "test/data/moons_labels.txt";
  DenseMat labels = stag::load_matrix(labels_filename);
  std::vector<StagInt> labels_vec(labels.rows());
  for (auto i = 0; i < labels.rows(); i++) {
    labels_vec.at(i) = (StagInt) labels(i, 0);
  }

  // Check the clustering performance with the asg roughly matches the
  // performance with the fully connected graph
  StagInt k = 2;
  std::vector<StagInt> clusters = stag::spectral_cluster(&asg, k);
  StagReal asg_ari = stag::adjusted_rand_index(clusters, labels_vec);

  stag::Graph sg = stag::similarity_graph(&data, a);
  clusters = stag::spectral_cluster(&sg, k);
  StagReal fc_ari = stag::adjusted_rand_index(clusters, labels_vec);
  EXPECT_GE(fc_ari, 0.9);
  EXPECT_GE(asg_ari, 0.8 * fc_ari);
}
