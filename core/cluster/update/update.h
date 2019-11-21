#ifndef UPDATE
#define UPDATE

#include <tuple>
#include <algorithm>

#include "../../../core/metric/metric.h"
#include "../../../core/utils/utils.h"

#endif

namespace cluster {
  namespace update {
    /* Update Algorithms for vectors */
    namespace vectors {  
      /** \brief Lloyd's Approach for Update:  
      **/
      template <typename T>
      std::vector<size_t> LloydsUpdate (const std::vector<T>& dataset_vectors, 
        const std::vector<size_t>& centers, const int& no_vectors, const int& vectors_dim, 
        const int& no_clusters, const std::vector<std::vector<size_t>>& clusters, 
        std::vector<T>& costs) {
          
        // Declare 1D vector which holds new best center
        std::vector<size_t> new_centers(no_clusters);
        for (size_t i = 0; i < no_clusters; i++) {
          for (auto& selected_center : clusters[i]) {
            T cost{};
            for (auto& id : clusters[i]) {
              if (id == selected_center) continue;
              T dist = metric::SquaredEuclidianDistance<T>(
                std::next(dataset_vectors.begin(), selected_center * vectors_dim),
                std::next(dataset_vectors.begin(), id * vectors_dim),
                std::next(dataset_vectors.begin(), id * vectors_dim + vectors_dim));
              cost += dist;
            }
            if (cost < costs[i]) {
              costs[i] = cost;
              new_centers[i] = selected_center;
            }
          }
        }
        return new_centers;
      }
    }
    /* Update Algorithms for curves */
    namespace curves {
      template<typename T>
      /** \brief Lloyd's Approach for Update: 
      **/
      std::vector<size_t> LloydsUpdate (
        const std::vector<std::pair<T,T>>& dataset_curves,std::vector<size_t>& centers, 
        const std::vector<int>& dataset_curves_lengths, 
        const std::vector<int>& dataset_curves_offsets, 
        const int& no_curves, const int& no_clusters, 
        const std::vector<std::vector<size_t>>& clusters, 
        std::vector<T>& costs) {
            
        // Declare 2D vector which holds all curves 
        // to their assigned cluster
        std::vector<size_t> new_centers(no_clusters);
        // Get lengths, offsets and points of curves
        // Warning: DTW does not work without const&
        // const std::vector<std::pair<T,T>>& center_curves = std::get<0>(centers);
        // const std::vector<int>& center_lengths = std::get<1>(centers);
        // const std::vector<int>& center_offsets = std::get<2>(centers);
        // Calculate all distances using DTW Distance
        T total_cost{};
        for (size_t i = 0; i < no_clusters; i++) {
          for (auto& center : clusters[i]) { 
            T cost{};
            for (auto& id : clusters[i]) {
              if (id == center) continue;
              T dist = metric::DTWDistance<T> (
              std::next(dataset_curves.begin(),dataset_curves_offsets[center]),
              std::next(dataset_curves.begin(),
                        dataset_curves_offsets[center] + dataset_curves_lengths[center]),
              std::next(dataset_curves.begin(),dataset_curves_offsets[id]),
              std::next(dataset_curves.begin(),
                        dataset_curves_offsets[id] + dataset_curves_lengths[id]));
              cost += dist;
            }
            // Store minimum dist from clusters
            if (cost < costs[i]) {
              costs[i] = cost;
              new_centers[i] = center;
            }
          }
        }
        return new_centers;
      }
    }
  }
}