#ifndef METRIC
#define METRIC

#include <utility>
#include <cmath>
#include <tuple>

#include "../utils/utils.h"

namespace metric {
  /** \brief Computes Manhattan Distance of 2 points in R^2
    @par const std::pair<T,T>& p - first point
    @par const std::pair<T,T>& q - second point
  */
  template <typename T>
  T _2DManhattanDistance(const std::pair<T,T>& p, const std::pair<T,T>& q) {
    T x_diff = std::abs((std::get<0>(p)-std::get<0>(q)));
    T y_diff = std::abs((std::get<1>(p)-std::get<1>(q)));
    return x_diff + y_diff;
  }
  /** \brief Computes Manhattan Distance of 2 points in R^d
    @par iterator p - iterator of the dataset point
    @par iterator q - iterator of the query point
    @par iterator q_end - end iterator of query point
      (each point is of dimension D so end of query point signals the
       final coordinate of dataset point as well)
    returns: manhattan_distance type of T
  */
  template <typename T, typename iterator>
  T ManhattanDistance(iterator p, iterator q, iterator q_end) {
    T manhattan_distance{};
    for (; q < q_end; ++q, ++p) {
      manhattan_distance += std::abs(*q - *p);
    }
    return manhattan_distance;
  }
  /** \brief Computes Euclidian Distance of 2 points in R^2
    @par const std::pair<T,T>& p - first point
    @par const std::pair<T,T>& q - second point
  */
  template <typename T>
  T _2DEuclidianDistance(const std::pair<T,T>& p, const std::pair<T,T>& q) {
    T x_diff = (std::get<0>(p)-std::get<0>(q)) * (std::get<0>(p)-std::get<0>(q));
    T y_diff = (std::get<1>(p)-std::get<1>(q)) * (std::get<1>(p)-std::get<1>(q));
    return sqrt(x_diff + y_diff);
  }
  /** \brief Computes Euclidian Distance of 2 points in R^d
    @par iterator p - iterator of the dataset point
    @par iterator q - iterator of the query point
    @par iterator q_end - end iterator of query point
      (each point is of dimension D so end of query point signals the
       final coordinate of dataset point as well)
    returns: euclidian_distance type of T
  */
  template <typename T, typename iterator>
  T SquaredEuclidianDistance(iterator p, iterator q, iterator q_end) {
    T manhattan_distance{};
    for (; q < q_end; ++q, ++p) {
      manhattan_distance += (*q - *p) * (*q - *p);
    }
    return manhattan_distance;
  }
  /** \brief Computes Dynamic Time Warping between two curves
    @par iterator p - iterator of the dataset curve
    @par iterator p_end - end iterator of dataset curve
    @par iterator q - iterator of the query curve
    @par iterator q_end - end iterator of query curve
      (each curve represents as a vector of pairs in which a pair
      represents a 2D-point)
  */
  template <typename T, typename iterator>
  T DTWDistance(iterator p, iterator p_end, iterator q, iterator q_end) {
    T dtw_distance{};
    // save up start iterator of query curve
    iterator q_start = q;
    // Get correspodent lengths of the two curves
    size_t N = std::distance(p,p_end);
    size_t M = std::distance(q,q_end);
    size_t dtw_size = (N + 1) * (M + 1);
    // Initialize the dtw 2D array using 1D notation
    auto dtw_array = new T [dtw_size];
    dtw_array[0] = 0;
    for (size_t i = 1; i < N + 1; ++i) {
      dtw_array[i * (M + 1)] = std::numeric_limits<T>::max();
    }
    for (size_t i = 1; i < M + 1; ++i) {
      dtw_array[i] = std::numeric_limits<T>::max();
    }
    // Compute the dtw distance using dynamic programming
    for (size_t i = 1; i < N + 1, p < p_end; ++i, ++p) {
      for (size_t j = 1; j < M + 1, q < q_end; ++j, ++q) {
        T dist = _2DEuclidianDistance(*p, *q);
        dtw_array[i * (M + 1) + j] = dist +
          utils::min(dtw_array[(i - 1) * (M + 1) + j],      // increment
            dtw_array[i * (M + 1) + j - 1],                 // deletion
            dtw_array[(i - 1) * (M + 1) + j - 1]);          // match

      }
      // reset query curve iterator
      q = q_start;
    }
    // Set result to a new variable and free memory allocated
    dtw_distance = dtw_array[dtw_size - 1];
    delete[] dtw_array;
    // Return dynamic time warping distance
    return dtw_distance;
  }
  /** \brief Computes average and max distance ratio appox_dists / exact_dists
    Each tuple consists of the nearest distance found, the id of the point/curve
    with minimum distance and the time taken to be computed
    @par const std::vector<std::tuple<T,U,double>> &exact - Pass by reference
      the vector of results produced by exact Nearest Neighbor for each query.
    @par const std::vector<std::tuple<T,U,double>> &approx - Pass by reference
      the vector of results produced by approximate Nearest Neighbor for
      each query.
  */
  template <typename T, typename U>
  std::tuple<double,double,int>
    EvaluationMetric(const std::vector<std::tuple<T,U,double>> &exact,
      const std::vector<std::tuple<T,U,double>> &approx) {

    // number of queries failed to be ansewered by A-NN
    int cnt_not_found{};
    int found{};
    // distance difference between exact and approximate NN
    double distance_error;
    // the sum of distance_error for each query
    double af{};
    // the averaged af by number of queries ansewered successfully
    double avg_af{};
    // number of queries
    size_t Q = exact.size();

    double max_af = std::numeric_limits<double>::min();
    for (size_t i = 0; i < Q; ++i) {
      // check the case for zero distance found (math warning by divsion with 0)
      if (std::get<0>(approx[i]) != 0 && std::get<0>(exact[i]) != 0) {
        // max distance found, so query failed to be ansewered
        if (std::get<0>(approx[i]) == std::numeric_limits<T>::max()) {
          cnt_not_found++;
          continue;
        }
        distance_error = (double) std::get<0>(approx[i]) / std::get<0>(exact[i]);
        if (distance_error > max_af) {
          max_af = distance_error;
        }
        af += distance_error;
      }
      if (std::get<0>(approx[i]) == std::get<0>(exact[i])) {
        found++;
      }
    }
    // For total number of queries discard the ones who failed to be ansewered
    avg_af = af / (Q - cnt_not_found) ;
    // return result as tuple of max Af, average Af and total #queries failed
    return std::make_tuple(max_af,avg_af,cnt_not_found);
  }
}

#endif
