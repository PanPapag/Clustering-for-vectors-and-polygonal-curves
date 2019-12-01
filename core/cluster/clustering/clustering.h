#ifndef CLUSTERING
#define CLUSTERING

#include "../initialization/initialization.h"
#include "../assignment/assignment.h"
#include "../update/update.h"

namespace cluster {
  namespace vectors {
    using namespace cluster::initialization::vectors;
    using namespace cluster::assignment::vectors;
    using namespace cluster::update::vectors;
    /**
      \brief Class Cluster representing k-means/k-medoids algorithm
      for vectors
    */
    template <typename T>
    class Cluster {
      private:
        /* class Cluster parameters */
        int no_clusters;
        int max_iter;
        std::string init;
        std::string assign;
        std::string update;
        /* Function pointers for the algorithms to be used */
        static constexpr auto p_init = ParkJunInit<T>;
        static constexpr auto p_assign = LloydsAssignment<T>;
        static constexpr auto p_update = LloydsUpdate<T>;
        /* Dataset info */
        std::vector<T> dataset_vectors;
        int no_vectors;
        int vectors_dim;
      public:
        /** \brief class Cluster constructor
          @par no_clusters: int, optional, default: 8
            The number of clusters to form as well as the number of
            centroids to generate.
          @par max_iter : int, default: 300
            Maximum number of iterations of the k-means/k-medoids algorithm for
            a single run.
          @par init : {‘k-means++’, ‘random’}
            Method for initialization, defaults to ‘k-means++’:
          @par assign : {‘lloyds’, ‘lsh’}
            Method for assignment, defaults to ‘lloyds’:
          @par update : {‘mean’, ‘pam’}
            Method for update, defaults to 'mean'
        */
        Cluster(int no_clusters = 8, int max_iter = 300,
          std::string init = "k-means++", std::string assign = "lloyds",
          std::string update = "mean") : no_clusters(no_clusters),
          max_iter(max_iter), init(init), assign(assign), update(update) {
          /* Pass function pointer to the initialization method */
          if (init == "random") {
            auto p_init = cluster::initialization::vectors::RandomInit<T>;
          }
          //TODO pass others too
        }
        /**
          \brief  Class Cluster default destructor
        */
        ~Cluster() = default;
        /** \brief Fit method stores dataset info for clustering
          @par[in] dv : vectors given from dataset
          @par[in] no_v : number of vectors
          @par[in] v_dim : vectors' dimensions (all vectors are dimensionally equal)
        */
        void Fit(const std::vector<T>& dv, const int& no_v, const int& v_dim) {
          dataset_vectors = dv;
          no_vectors = no_v;
          vectors_dim = v_dim;
          /**
            If assignment step is executed using lsh range search, map dataset
            to corresponding lsh structures
          */
          if (assign == "lsh") {
            std::cout << "TODO" << std::endl;
          }
        }
        /**
          \brief Predict the closest cluster each sample in dataset belongs to
          returns a vector of centroids coordinates and a vector of vectors
          which stores in each position the indexes of dataset_vectors assigned
          in this cluster
        */
        std::pair<std::vector<T>,std::vector<std::vector<size_t>>> Predict(void) {
          /* Declare types */
          std::vector<T> centroids;
          std::tuple<std::vector<std::vector<size_t>>,std::vector<T>> clusters;
          /* At first initialize centroids */
          centroids = p_init(dataset_vectors, no_vectors, vectors_dim,
                             no_clusters);
          /* Calculate clusters and update centroids max_iter times */
          for (size_t i = 0; i < max_iter; ++i) {
            /* Assigment step */
            clusters = p_assign(dataset_vectors, centroids, no_vectors,
                                vectors_dim, no_clusters);
            /* Update step */
            centroids = p_update(dataset_vectors, centroids,
                                 no_vectors, vectors_dim, no_clusters,
                                 std::get<0>(clusters), std::get<1>(clusters));
          }
          return std::make_pair(centroids,std::get<0>(clusters));
        }
    };
  }
}

#endif
