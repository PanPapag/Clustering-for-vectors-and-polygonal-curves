#ifndef CLUSTER_IO_UTILS
#define CLUSTER_IO_UTILS

#include <iostream>
#include <fstream>
#include <tuple>
#include <vector>

namespace utils {
  namespace io {
    namespace vectors {
      /** \brief ReadFile - Reads file provided by user and
        @par std::string file_name - Pass by reference the path to the input file
        @par std::vector<T> &vectors - Pass by reference a vector type T which
             represent the N vectors of dimension D
        @par std::vector<K> &ids - Pass by reference a vector type K which stores
             vectors' ids
        @par const int no_vectors - Number of point in file
        @par const int dim - Vectors' dimension
        @par ExitCode &statues - enumerated ExitCode provided from namespace utils
        return: SUCCESS or FAIL
      */
      template <typename T, typename K>
      int ReadFile(std::string &file_name, const int no_vectors, const int dim,
        std::vector<T> &vectors, std::vector<K> &ids, utils::ExitCode &status) {
        // Open file
        std::ifstream infile;
        infile.open(file_name);
        // Check if file exists
        if (!infile) {
          status = INVALID_DATASET;
          return FAIL;
        }
        // skip first line
        std::string temp_fl;
        infile >> temp_fl;
        // Read file
        for (size_t i = 0; i < no_vectors && infile; ++i) {
          infile >> ids[i];
          for (size_t j = 0; j < dim; ++j) {
            infile >> vectors[i * dim + j];
          }
        }
        // close the file
        infile.close();
        return SUCCESS;
      }
      /** \brief GetNoDatasetVectors - Get the number of vectors in the dataset
        @par file_name - Relative path to the dataset
        @par no_vectors - Total number of vectors in the dataset to be returned
        @par ExitCode &statues - enumerated ExitCode provided from namespace utils
        return: SUCCESS or FAIL
      */
      int GetNoDatasetVectors(std::string &file_name, uint32_t &no_vectors,
         ExitCode &status);
      /** \brief GetVectorsDim - Get the dimension of the vectors in the dataset
        @par std::string &file_name - Path to the dataset
        @par int &dim - Pass by reference the vectors' dimension
        @par ExitCode &statues - enumerated ExitCode provided from namespace utils
        return: SUCCESS or FAIL
      */
      int GetVectorsDim(std::string &file_name, uint16_t &dim, ExitCode &status);
    }
    /** \brief GetFirstLine - Reads first line of the input file
      @par file_name - Relative path to the dataset
      @par first_line - First line string to be returned
      @par ExitCode &statues - enumerated ExitCode provided from namespace utils
      return: SUCCESS or FAIL
    */
    int GetFirstLine(std::string &file_name, std::string &first_line,
      utils::ExitCode &status);
  }
}

#endif
