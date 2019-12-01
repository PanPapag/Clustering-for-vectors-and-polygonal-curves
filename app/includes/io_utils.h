#ifndef CLUSTER_IO_UTILS
#define CLUSTER_IO_UTILS

#include <iostream>
#include <fstream>
#include <tuple>
#include <vector>
#include <algorithm>

#include "../../core/utils/utils.h"


namespace utils {
  namespace io {
    /** \brief GetFirstLine - Reads first line of the input file
      @par file_name - Relative path to the dataset
      @par first_line - First line string to be returned
      @par ExitCode &statues - enumerated ExitCode provided from namespace utils
      return: SUCCESS or FAIL
    */
    int GetFirstLine(std::string &file_name, std::string &first_line,
      utils::ExitCode &status);
    /** \brief ReadConfig - Reads configuration file
      @par file_name - Relative path to the dataset
      @par no_clusters - Total number of clusters
      @par no_grids - Total number of grids used by grid LSh
      @par no_hf - Total number of LSH hash function
      @par no_ht - Total number of LSH hash tables
      @par ExitCode &statues - enumerated ExitCode provided from namespace utils
      return: SUCCESS or FAIL
    */
    int ReadConfig(std::string &file_name, uint8_t& no_clusters,
      uint8_t& no_grids, uint8_t& no_hf, uint8_t& no_ht, utils::ExitCode &status);
    /* Vectors I/O utils */
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
      int ReadFile(const std::string &file_name, const int no_vectors,
        const int dim, std::vector<T> &vectors, std::vector<K> &ids,
        utils::ExitCode &status) {
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
      /** \brief WriteFile - Output prorgam results to the given output file

        return: SUCCESS or FAIL
      */
      template <typename T, typename U>
      int WriteFile(const std::string &file_name, utils::ExitCode &status) {

        // Open file
        std::ofstream outfile;
        outfile.open(file_name);
        // Check if file is opened
        if (outfile.is_open()) {

        } else {
          status = INVALID_OUTPUT;
          return FAIL;
        }
        // close file
        outfile.close();
        return SUCCESS;
      }
      /** \brief GetNoDatasetVectors - Get the number of vectors in the dataset
        @par file_name - Relative path to the dataset
        @par no_vectors - Total number of vectors in the dataset to be returned
        @par ExitCode &statues - enumerated ExitCode provided from namespace utils
        return: SUCCESS or FAIL
      */
      int GetNoDataVectors(std::string &file_name, uint32_t &no_vectors,
         ExitCode &status);
      /** \brief GetVectorsDim - Get the dimension of the vectors in the dataset
        @par std::string &file_name - Path to the dataset
        @par int &dim - Pass by reference the vectors' dimension
        @par ExitCode &statues - enumerated ExitCode provided from namespace utils
        return: SUCCESS or FAIL
      */
      int GetVectorsDim(std::string &file_name, uint16_t &dim, ExitCode &status);
    }
    /* Curves I/O utils */
    namespace curves {
      /** \brief ReadFile - Reads file provided by user
        @par std::string& file_name - Pass by reference the path to the input file
        @par std::vector<T>& curves - Pass by reference a vector type T which
             represent the N curves of length m_i, i = 1..N
        @par std::vector<K>& ids - Pass by reference a vector type K which stores
             curves' ids
        @par std::vector<T>& lengths - Pass by reference a vector type int which
             stores curves' lenths
        @par std::vector<T>& offsets - Pass by reference a vector type int which
             stores offesets to the std::vector<std::pair<T,T>>& curves
        @par const int no_curves - Number of point in file
        @par ExitCode& statues - enumerated ExitCode provided from namespace utils
        return: SUCCESS or FAIL
      */
      template <typename T, typename K>
      int ReadFile(std::string& file_name, const int no_curves,
        std::vector<std::pair<T,T>>& curves, std::vector<K>& ids,
        std::vector<int>& lengths, std::vector<int>& offsets,
        utils::ExitCode& status) {

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
        std::string buffer_x;
        std::string buffer_y;
        int offset{};
        for (size_t i = 0; i < no_curves && infile; ++i) {
          infile >> ids[i];
          infile >> lengths[i];
          for (size_t j = 0; j < lengths[i]; ++j) {
            infile >> buffer_x;
            infile >> buffer_y;
            // remove noise characters such as '(', to convert successuflly to T
            buffer_x.erase(std::remove(buffer_x.begin(),buffer_x.end(),'('),
                           buffer_x.end());
            T point_x = convert_to<T>(buffer_x);
            T point_y = convert_to<T>(buffer_y);
            curves.push_back(std::make_pair(point_x,point_y));
          }
          offsets[i] = offset;
          offset += lengths[i];
        }
        // close the file
        infile.close();
        return SUCCESS;
      }
      /** \brief GetNoDataCurves - Get the number of curves in the dataset
        @par file_name - Relative path to the dataset
        @par no_vectors - Total number of curvess in the dataset to be returned
        @par ExitCode &statues - enumerated ExitCode provided from namespace utils
        return: SUCCESS or FAIL
      */
      int GetNoDataCurves(std::string& file_name, uint32_t& no_curves,
        ExitCode& status);
    }
  }
}

#endif
