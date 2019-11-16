#ifndef CLUSTER_IO_UTILS
#define CLUSTER_IO_UTILS

#include <iostream>
#include <fstream>
#include <tuple>
#include <vector>

namespace utils {
  namespace io {
    /** \brief ReadFile - Reads file provided by user and
      @par std::string file_name - Pass by reference the path to the input file
      @par std::vector<T> &points - Pass by reference a vector type T which
           represent the N points of dimension D
      @par std::vector<K> &ids - Pass by reference a vector type K which stores
           points' ids
      @par const int no_points - Number of point in file
      @par const int dim - Points' dimension
      @par ExitCode &statues - enumerated ExitCode provided from namespace utils
      return: SUCCESS or FAIL
    */
    template <typename T, typename K>
    int ReadFile(std::string &file_name, const int no_points, const int dim,
       std::vector<T> &points, std::vector<K> &ids, utils::ExitCode &status) {
      // Open file
      std::ifstream infile;
      infile.open(file_name);
      // Check if file exists
      if (!infile) {
        status = INVALID_DATASET;
        return FAIL;
      }
      // Read file
      for (size_t i = 0; i < no_points && infile; ++i) {
        infile >> ids[i];
        for (size_t j = 0; j < dim; ++j) {
          infile >> points[i * dim + j];
        }
      }
      // close the file
      infile.close();
      return SUCCESS;
    }
    /** \brief WriteFile - Output prorgam results to given output file
      @par std::string &file_name - Pass by reference the path to the output file
      @par std::vector<std::tuple<T,U,double>> &exact - Results from Brute Force
      @par std::vector<std::tuple<T,U,double>> &approx - Results from LSH
      @par std::vector<std::vector<std::pair<T,U>>> &radius_nn - Results from radius NN
      @par ExitCode &statues - enumerated ExitCode provided from namespace utils
      return: SUCCESS or FAIL
    */
    template <typename T, typename U>
    int WriteFile(std::string &file_name, std::vector<std::tuple<T,U,double>> &exact,
      std::vector<std::tuple<T,U,double>> &approx, const double R,
      std::vector<std::vector<std::pair<T,U>>> &radius_nn, utils::ExitCode &status) {

      // Open file
      std::ofstream outfile;
      outfile.open(file_name);
      // Check if file is opened
      if (outfile.is_open()) {
        /* Get number of queries executed */
        int N = exact.size();
        /* Print info for each query */
        for (size_t i = 0; i < N; ++i) {
          outfile << "Query: " << i << std::endl;
          outfile << "Nearest neighbor: " << std::get<1>(approx[i]) << std::endl;
          outfile << "distanceLSH: " <<  std::get<0>(approx[i]) << std::endl;
          outfile << "distanceTrue: " <<  std::get<0>(exact[i]) << std::endl;
          outfile << "tLSH: " <<  std::get<2>(approx[i]) << " seconds" << std::endl;
          outfile << "tTrue: " <<  std::get<2>(exact[i]) << " seconds" << std::endl;
          if (R != 0.0) {
            outfile << R << "-near neighbors: " << std::endl;
            for (int j = 0; j < radius_nn[i].size(); ++j) {
              outfile << std::get<1>(radius_nn[i][j]) << std::endl;
            }
            if (!radius_nn[i].size()) {
              outfile << "No " << R << "-near neighbors found" << std::endl;
            }
          }
          outfile << std::endl;
        }

      } else {
        status = INVALID_OUTPUT;
        return FAIL;
      }
      // close file
      outfile.close();
      return SUCCESS;
    }
    /** \brief GetDataPoints - Get the number of file data points
      @par std::string &file_name - Path to file
      @par int &no_vectors - Pass by reference the number of vectors to be returned
      @par ExitCode &statues - enumerated ExitCode provided from namespace utils
      return: SUCCESS or FAIL
    */
    int GetDataPoints(std::string &file_name, uint32_t &no_vectors, ExitCode &status);
    /** \brief GetPointsDim - Get the dimension of the points
      @par std::string &file_name - Path to file
      @par int &dim - Pass by reference the points' dimension
      @par ExitCode &statues - enumerated ExitCode provided from namespace utils
      return: SUCCESS or FAIL
    */
    int GetPointsDim(std::string &file_name, uint16_t &dim, ExitCode &status);
  }
}

#endif
