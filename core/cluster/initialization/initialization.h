#ifndef INITIALIZATION
#define INITIALIZATION

#include "../../../core/metric/metric.h"
#include "../../../core/utils/utils.h"

namespace cluster {

  namespace initialization {

    namespace vectors {
      /** \brief Initilization of cluster centroids
        Details:
          (1) Calculate symmetric n × n distance matrix of all objects,
          i.e. all distances d_ij from every object i = 1,...,n to every other
          object j = 1,...,n,i != j.
          (2) For object i compute:
            v_i = sum from j = 1 to n (d_ij / sum from t = 1 to n d_jt)
            for i = 1,...,n
          (3) Return the k objects with k smallest vi values.
          Algorithm proposed in [Park-Jun’09]
      */
      template <typename T>
      std::vector<T> ParkJunInit(const std::vector<T>& dataset_vectors,
        const int& no_vectors, const int& vectors_dim, const int& no_clusters) {

        // Declare 2D d array which holds all distances between the vectors
        T d_array[no_vectors][no_vectors];
        // Calculate all distances using ManhattanDistance
        for (size_t i = 0; i < no_vectors; ++i) {
          for (size_t j = 0; j < no_vectors; ++j) {
            if (i == j) {
              d_array[i][j] = {};
            } else {
              d_array[i][j] = metric::ManhattanDistance<T>(
                std::next(dataset_vectors.begin(), i * vectors_dim),
                std::next(dataset_vectors.begin(), j * vectors_dim),
                std::next(dataset_vectors.begin(), j * vectors_dim + vectors_dim));
            }
          }
        }
        // Declare 1D v array and computes its values as the algorithm above proposes
        std::vector<T> v(no_vectors);
        for (size_t i = 0; i < no_vectors; ++i) {
          T o_sum{};
          for (size_t j = 0; j < no_vectors; ++j) {
            T i_sum{};
            for (size_t t = 0; t < no_vectors; ++t) {
              i_sum += d_array[j][t];
            }
            o_sum += d_array[i][j] / i_sum;
          }
          v[i] = o_sum;
        }
        // Find indices of no_clusters smallest v_array values
        std::vector<int> idx = utils::ArgMin<T>(v, no_clusters);
        // Initialize centroids
        std::vector<T> centroids(no_clusters * vectors_dim);
        for (size_t i = 0; i < no_clusters; ++i) {
          for (size_t j = 0; j < vectors_dim; ++j) {
            centroids[i * vectors_dim + j] = dataset_vectors[idx[i] * vectors_dim + j];
          }
        }
        // Return Initialized centroids
        return centroids;
      }
    }

    namespace curves {

    }
  }
}


#endif
