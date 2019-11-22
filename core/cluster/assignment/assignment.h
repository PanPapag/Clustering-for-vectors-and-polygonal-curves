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
        std::tuple<std::vector<std::vector<size_t>>,std::vector<T>> LloydsAssignment (
          const std::vector<T>& dataset_vectors, const std::vector<T>& centers, 
          const int& no_vectors, const int& vectors_dim, const int& no_clusters) {
          
          // Declare 2D vector which holds all vectors 
          // to their assigned cluster. 
          std::vector<std::vector<size_t>> d_array(no_clusters);
          std::vector<T> assign_costs(no_clusters, 0);
          // Calculate all distances using Euclidean Distance
          T total_cost{};
          for (size_t i = 0; i < no_vectors; i++) {
            T min_dist{};
            size_t assigned_cluster{};
            for (size_t j = 0; j < no_clusters; j++) {
              T dist = metric::SquaredEuclidianDistance<T>(
                std::next(dataset_vectors.begin(), i * vectors_dim),
                std::next(dataset_vectors.begin(), j * vectors_dim),
                std::next(dataset_vectors.begin(), j * vectors_dim + vectors_dim));
              // Store minimum dist from clusters
              if (dist < min_dist || j==0) {
                min_dist = dist;
                assigned_cluster = j;
              }
            }
            // Store to closest cluster
            d_array[assigned_cluster].push_back(i);
            assign_costs[assigned_cluster] += min_dist;
          }
          return std::make_tuple(d_array,assign_costs);
        }
      }
      /* Assignment Algorithms for curves */
      namespace curves {
        template<typename T>
        /** \brief Lloyd's Approach for Assignment: Each curve is assigned to 
         * its closest cluster.
         * Returns a 2D vector which holds the id of curves assigned to each cluster.
         * We use DTW to calculate minimum distances. 
        **/
        std::tuple<std::vector<std::vector<size_t>>,std::vector<T>> 
          LloydsAssignment (const std::vector<std::pair<T,T>>& dataset_curves, 
            std::tuple<std::vector<std::pair<T,T>>,std::vector<int>,std::vector<int>>& centers, 
            const std::vector<int>& dataset_curves_lengths, 
            const std::vector<int>& dataset_curves_offsets, 
            const int& no_curves, const int& no_clusters) {
            
            // Declare 2D vector which holds all curves 
            // to their assigned cluster
            std::vector<std::vector<size_t>> d_array(no_clusters);
            std::vector<T> assign_costs(no_clusters, 0);
            // Get lengths, offsets and points of curves
            // Warning: DTW does not work without const&
            const std::vector<std::pair<T,T>>& center_curves = std::get<0>(centers);
            const std::vector<int>& center_lengths = std::get<1>(centers);
            const std::vector<int>& center_offsets = std::get<2>(centers);
            // Calculate all distances using DTW Distance
            T total_cost{};
            for (size_t i = 0; i < no_curves; i++) {
              T min_dist{};
              size_t assigned_cluster{};
              for (size_t j = 0; j < no_clusters; j++) {
                T dist = metric::DTWDistance<T> (
                  std::next(dataset_curves.begin(),dataset_curves_offsets[i]),
                  std::next(dataset_curves.begin(),
                            dataset_curves_offsets[i] + dataset_curves_lengths[i]),
                  std::next(center_curves.begin(),center_offsets[j]),
                  std::next(center_curves.begin(),
                            center_offsets[j] + center_lengths[j]));
                // Store minimum dist from clusters
                if (dist < min_dist || j==0) {
                  min_dist = dist;
                  assigned_cluster = j;
                }
              }
              // Store id to closest cluster
              d_array[assigned_cluster].push_back(i);
              assign_costs[assigned_cluster] += min_dist;
            }
            return std::make_tuple(d_array,assign_costs);
          }
      }
  }
}