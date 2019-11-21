#ifndef ASSIGNMENT
#define ASSIGNMENT

#include <tuple>
#include <algorithm>

#include "../../../core/metric/metric.h"
#include "../../../core/utils/utils.h"

#endif

namespace cluster {
  namespace assignment {
      /* Assignment Algorithms for vectors */
      namespace vectors {  
        /** \brief Lloyd's Approach for Assignment: Each vector is assigned to 
         * its closest cluster.
         * Returns a 2D vector which holds the id of vectors assigned to each cluster.
         * We use Euclidean Distance to calculate minimum distances. 
        **/
        template <typename T>
        std::vector<std::vector<size_t>> LloydsAssignment(const std::vector<T>& dataset_vectors, 
            const std::vector<T>& centers, const int& no_vectors, const int& vectors_dim, 
            const int& no_clusters) {
          
          // Declare 2D vector which holds all vectors 
          // to their assigned cluster
          std::vector<std::vector<size_t>> d_array(no_clusters);
          // Calculate all distances using Euclidean Distance
          for (size_t i = 0; i < no_vectors; i++) {
            T min_dist{};
            size_t assigned_cluster{};
            for (size_t j = 0; j < no_clusters; j++) {
              T dist = metric::SquaredEuclidianDistance<T>(
                std::next(dataset_vectors.begin(), i * vectors_dim),
                std::next(centers.begin(), j * vectors_dim),
                std::next(centers.begin(), j * vectors_dim + vectors_dim));
              // Store minimum dist from clusters
              if (dist < min_dist || j==0) {
                min_dist = dist;
                assigned_cluster = j;
              }
            }
            // Store to closest cluster
            d_array[assigned_cluster].push_back(i);
          }
        }
      }
      namespace curves {

      }
  }
}