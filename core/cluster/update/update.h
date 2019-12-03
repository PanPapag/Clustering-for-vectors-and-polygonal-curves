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
          for (const auto& selected_center : clusters[i]) {
            T cost{};
            // calculate dist from all other vectors in cluster
            for (const auto& id : clusters[i]) {
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
          const std::vector<T>& centers, const int& no_vectors,
          const int& vectors_dim, const int& no_clusters,
          const std::vector<std::vector<size_t>>& clusters) {

          std::vector<T> new_centers(no_clusters * vectors_dim, 0);
          for (size_t i = 0; i < no_clusters; ++i) {
            for (const auto& id : clusters[i]) {
              for (size_t j = 0; j < vectors_dim; ++j) {
                new_centers[i * vectors_dim + j] += dataset_vectors[id * vectors_dim + j];
              }
            }
            for (size_t j = 0; j < vectors_dim; ++j) {
              new_centers[i * vectors_dim + j] /= clusters[i].size();
            }
          }
          return new_centers;
        }
      }

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
          std::vector<int> offsets(no_clusters,-1);
          // Calculate all distances using DTW Distance
          T total_cost{};
          for (size_t i = 0; i < no_clusters; i++) {
            for (const auto& center : clusters[i]) {
              T cost{};
              for (const auto& id : clusters[i]) {
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
          std::vector<std::pair<T,T>> prev_centers = std::get<0>(centers);
          std::vector<int> prev_offsets = std::get<2>(centers);
          std::vector<int> prev_lengths = std::get<1>(centers);
          for (size_t i = 0; i < no_clusters; ++i) {
            if(offsets[i] != -1) {
              new_centers_lengths[i] = dataset_curves_lengths[offsets[i]];
              for (size_t j = 0; j < new_centers_lengths[i]; ++j) {
                new_centers.push_back(dataset_curves[dataset_curves_offsets[offsets[i]] + j]);
              }
            } else {
              new_centers_lengths[i] = prev_lengths[i];
              for (size_t j = 0; j < new_centers_lengths[i]; ++j) {
                new_centers.push_back(prev_centers[prev_offsets[i] + j]);
              }
            }
            new_centers_offsets[i] = off;
            off += new_centers_lengths[i];
          }
          return std::make_tuple(new_centers,new_centers_lengths,new_centers_offsets);
      }
      /**
        \brief Lloyd's Approach for Update
      **/
      template <typename T>
      //std::tuple<std::vector<std::pair<T,T>>,std::vector<int>,std::vector<int>>
      void
        LloydsUpdate (const std::vector<std::pair<T,T>>& dataset_curves,
          std::tuple<std::vector<std::pair<T,T>>,std::vector<int>,std::vector<int>>& centers,
          const std::vector<int>& dataset_curves_lengths,
          const std::vector<int>& dataset_curves_offsets,
          const int& no_curves, const int& no_clusters,
          const std::vector<std::vector<size_t>>& clusters,
          std::vector<T>& costs) {

          // Compute new centers' lengths λ
          // center's length equals to the avg 
          // lengths of curves in cluster aka 
          // λ(i) = (1/n) Sum [S(i)]
          std::vector<int> lamdas(no_clusters, 0);
          for (size_t i = 0; i < no_clusters; i++) {
            int n = 0;
            for (const auto& center : clusters[i]) {
              lamdas[i] += dataset_curves_lengths[center];
              n++;
            }
            lamdas[i] /= n;
          }
          // Choose randomly a curve with length 
          // greater or equal than lambda[i] 
          // and set it as center C
          std::vector<int> rand_centers(no_clusters);
          for (size_t i = 0; i < no_clusters; i++) {
            std::vector<int> lamda_curves;
            for (const auto& center : clusters[i]) {
              if (dataset_curves_lengths[center] >= lamdas[i]) {
                lamda_curves.push_back(center);
              }
            }
            int idx = rand() % lamda_curves.size();
            rand_centers[i] = lamda_curves[idx];
          }

          // for(size_t i = 0; i < no_clusters; i++) {
          //   std::cout << lamdas[i] << " " << dataset_curves_lengths[rand_centers[i]] << std::endl;
          // }

          int off{};
          std::vector<std::pair<T,T>> c;
          std::vector<int> c_offsets(no_clusters);
          std::vector<int> c_lengths(no_clusters);
          for (size_t i = 0; i < no_clusters; ++i) {
            size_t start = rand() % (dataset_curves_lengths[rand_centers[i]]-lamdas[i]);
            c_lengths[i] = lamdas[i];
            //std::cout << start << " " << c_lengths[i]+start << " " << dataset_curves_lengths[rand_centers[i]] << std::endl; 
            for (size_t j = start; j < c_lengths[i]+start; ++j) {
              c.push_back(dataset_curves[dataset_curves_offsets[rand_centers[i]] + j]);
            }
            c_offsets[i] = off;
            off += c_lengths[i];
          }

          for (size_t i = 0; i < no_clusters; i++) {
            while(1) {
              std::vector<std::pair<T,T>> Ai (lamdas[i],std::make_pair(0,0));
              std::vector<int> Ni (lamdas[i],0);
              const std::vector<std::pair<T,T>>& temp_centers = c;
              const std::vector<int>& temp_centers_offsets = c_offsets;
              const std::vector<int>& temp_centers_lengths = c_lengths;
              for (const auto& id : clusters[i]) {
                std::vector<std::pair<T,T>> vec = metric::BestTraversal<T>(
                std::next(temp_centers.begin(),temp_centers_offsets[i]),
                std::next(temp_centers.begin(),
                          temp_centers_offsets[i] + temp_centers_lengths[i]),
                std::next(dataset_curves.begin(),dataset_curves_offsets[id]),
                std::next(dataset_curves.begin(),
                          dataset_curves_offsets[id] + dataset_curves_lengths[id]));
                for(size_t j=0; j<vec.size(); j++) {
                  size_t pos_i = vec[j].first;
                  size_t pos_j = vec[j].second;
                  T curr_val_i = Ai[pos_i].first;
                  T curr_val_j = Ai[pos_i].second; 
                  std::pair<T,T> val = dataset_curves[dataset_curves_offsets[id] + pos_j];
                  Ai[pos_i] = std::make_pair(curr_val_i + val.first, curr_val_j + val.second);
                  Ni[pos_i]++; 
                }
              }
            }
          }

          T prev_dist{};
          while(1) {  
            std::vector<std::vector<std::pair<T,T>>> A;
            std::vector<std::vector<int>> N;
            const std::vector<std::pair<T,T>>& temp_centers = c;
            const std::vector<int>& temp_centers_offsets = c_offsets;
            const std::vector<int>& temp_centers_lengths = c_lengths;
            for (size_t i = 0; i < no_clusters; i++) {
              std::vector<std::pair<T,T>> Ai (lamdas[i],std::make_pair(0,0));
              std::vector<int> Ni (lamdas[i],0);
              for (const auto& id : clusters[i]) {
                std::vector<std::pair<T,T>> vec = metric::BestTraversal<T>(
                std::next(temp_centers.begin(),temp_centers_offsets[i]),
                std::next(temp_centers.begin(),
                          temp_centers_offsets[i] + temp_centers_lengths[i]),
                std::next(dataset_curves.begin(),dataset_curves_offsets[id]),
                std::next(dataset_curves.begin(),
                          dataset_curves_offsets[id] + dataset_curves_lengths[id]));
                for(size_t j=0; j<vec.size(); j++) {
                  size_t pos_i = vec[j].first;
                  size_t pos_j = vec[j].second;
                  T curr_val_i = Ai[pos_i].first;
                  T curr_val_j = Ai[pos_i].second; 
                  std::pair<T,T> val = dataset_curves[dataset_curves_offsets[id] + pos_j];
                  Ai[pos_i] = std::make_pair(curr_val_i + val.first, curr_val_j + val.second);
                  Ni[pos_i]++; 
                }
              }
              A.push_back(Ai);
              N.push_back(Ni);
            }
            c.clear();
            c_offsets.clear();
            c_lengths.clear();
            int off{};
            for (size_t i = 0; i < no_clusters; i++) {
              c_lengths[i] = lamdas[i];
              for (size_t j = 0; j < c_lengths[i]; j++) {
                c.push_back(std::make_pair(A[i][j].first/N[i][j],A[i][j].second/N[i][j]));
              }
              c_offsets[i] = off;
              off += c_lengths[i];
            }
            T total_dist{};
            for (size_t i = 0; i < no_clusters; i++) {
              T dist = metric::DTWDistance<T> (
                std::next(c.begin(),c_offsets[i]),
                std::next(c.begin(),
                          c_offsets[i] + c_lengths[i]),
                std::next(temp_centers.begin(),temp_centers_offsets[i]),
                std::next(temp_centers.begin(),
                          temp_centers_offsets[i] + temp_centers_lengths[i]));
              total_dist += dist;
            }
            if (total_dist <= prev_dist*0.1) break;
            prev_dist = total_dist;
          }
          //return std::make_tuple(new_centers,new_centers_lengths,new_centers_offsets);
      }
    }
  }
}

#endif
