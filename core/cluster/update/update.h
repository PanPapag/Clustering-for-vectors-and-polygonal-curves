#ifndef UPDATE
#define UPDATE

#include <tuple>
#include <algorithm>

#include "../../../core/metric/metric.h"
#include "../../../core/utils/utils.h"

namespace cluster {
  namespace update {
    /* Update Algorithms for vectors */
    namespace vectors {
      /** \brief PAM Algorith a'la Lloyd's
       * Details:
       * Search for potential new centers in clusters
        (1) Calculate all the distances of each vector in a cluster
            from all the other vectors of same cluster
        (2) Set as new center of cluster the vector 
            with smaller sum of distances
        @par[in] const std::vector<T>& dataset_vectors : vectors given from dataset
        @par[in] const std::vector<T>& centers : selected centers of each cluster
        @par[in] const int& no_vectors : number of vectors
        @par[in] const int& vectors_dim : vectors' dimensions (all vectors are
                                          dimensionally equal)
        @par[in] const int& no_clusters : number of clusters
        @par[in] const std::vector<T>& costs : assignment cost of each cluster's center
                                                computed in Assignment Method
      **/
      template <typename T>
      std::vector<T> PAMUpdate (const std::vector<T>& dataset_vectors,
        const std::vector<T>& centers, const int& no_vectors, const int& vectors_dim,
        const int& no_clusters, const std::vector<std::vector<size_t>>& clusters,
        std::vector<T>& costs) {

        // Declare 1D vector which holds new best centers ids
        // We initialize vector with -1 in case centers remain the same
        std::vector<ssize_t> new_centers_offsets(no_clusters,-1);
        // Declare 1D vector which holds new best centers
        std::vector<T> new_centers(no_clusters * vectors_dim);
        for (size_t i = 0; i < no_clusters; i++) {
          // iterate over each vector in cluster
          for (auto& selected_center : clusters[i]) {
            T cost{};
            // calculate dist from all other vectors in cluster
            for (auto& id : clusters[i]) {
              // skip calculation - dist would be 0
              if(id == selected_center) continue;
              T dist = metric::SquaredEuclidianDistance<T>(
                std::next(dataset_vectors.begin(), selected_center * vectors_dim),
                std::next(dataset_vectors.begin(), id * vectors_dim),
                std::next(dataset_vectors.begin(), id * vectors_dim + vectors_dim));
              cost += dist;
            }
            // update costs and new centers' ids
            if (cost < costs[i]) {
              costs[i] = cost;
              new_centers_offsets[i] = selected_center;
            }
          }
        }
        for (size_t i = 0; i < no_clusters; ++i) {
          for (size_t j = 0; j < vectors_dim; ++j) {
            if (new_centers_offsets[i] != -1) {
              new_centers[i * vectors_dim + j] = dataset_vectors[new_centers_offsets[i] * vectors_dim + j];
            } else {
              // center remained the same
              new_centers[i * vectors_dim + j] = centers[i * vectors_dim + j];
            }
          }
        }
        return new_centers;
      }

      /** \brief Lloyd's Algorithm for Update
        * Details:
        * Compute Mean and set it as new center in each cluster: 
        * m(i) = (1/T) * Sum v(i)
          @par[in] const std::vector<T>& dataset_vectors : vectors given from dataset
          @par[in] const std::vector<T>& centers : selected centers of each cluster
          @par[in] const int& no_vectors : number of vectors
          @par[in] const int& vectors_dim : vectors' dimensions (all vectors are
                                            dimensionally equal)
          @par[in] const int& no_clusters : number of clusters
        **/
        template <typename T>
        std::vector<T> LloydsUpdate (const std::vector<T>& dataset_vectors,
          const std::vector<T>& centers, const int& no_vectors, const int& vectors_dim,
          const int& no_clusters, const std::vector<std::vector<size_t>>& clusters,
          std::vector<T>& costs) {

          std::vector<T> new_centers(no_clusters * vectors_dim, 0);
          for (size_t i = 0; i < no_clusters; i++) {
            T cost{};
            for (auto& id : clusters[i]) {
              for(size_t j=0; j<vectors_dim; j++) {
                new_centers[i * vectors_dim + j] += dataset_vectors[id * vectors_dim + j];
              }
            }
            for(size_t j=0; j<vectors_dim; j++) {
              new_centers[i * vectors_dim + j] /= clusters[i].size();
            }
          }
          return new_centers;
        }
      }
    
    
    // TODO: Fix bug in offset calculation
    /* Update Algorithms for curves */
    namespace curves {
      /**
        \brief Lloyd's Approach for Update
      **/
      template <typename T>
      std::tuple<std::vector<std::pair<T,T>>,std::vector<int>,std::vector<int>>
        PAMUpdate (const std::vector<std::pair<T,T>>& dataset_curves,
          std::tuple<std::vector<std::pair<T,T>>,std::vector<int>,std::vector<int>>& centers,
          const std::vector<int>& dataset_curves_lengths,
          const std::vector<int>& dataset_curves_offsets,
          const int& no_curves, const int& no_clusters,
          const std::vector<std::vector<size_t>>& clusters,
          std::vector<T>& costs) {

          // Declare 2D vector which holds all curves
          // to their assigned cluster
          std::vector<size_t> offsets(no_clusters);
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
                offsets[i] = center;
              }
            }
          }

          int off{};
          std::vector<std::pair<T,T>> new_centers;
          std::vector<int> new_centers_offsets(no_clusters);
          std::vector<int> new_centers_lengths(no_clusters);
          for (size_t i = 0; i < no_clusters; ++i) {
            new_centers_lengths[i] = dataset_curves_lengths[offsets[i]];
            for (size_t j = 0; j < new_centers_lengths[i]; ++j) {
              new_centers.push_back(dataset_curves[dataset_curves_offsets[offsets[i]] + j]);
            }
            new_centers_offsets[i] = off;
            off += new_centers_lengths[i];
          }
          return std::make_tuple(new_centers,new_centers_lengths,new_centers_offsets);
      }
    }
  }
}

#endif
