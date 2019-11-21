#include <algorithm>
#include <chrono>
#include <cmath>
#include <iostream>
#include <iterator>
#include <random>
#include <string>
#include <tuple>
#include <utility>
#include <vector>

#include "../../core/cluster/initialization/initialization.h"
#include "../../core/cluster/assignment/assignment.h"
#include "../../core/cluster/update/update.h"
#include "../../core/metric/metric.h"
#include "../../core/utils/utils.h"

#include "../includes/args_utils.h"
#include "../includes/io_utils.h"
#include "../includes/report_utils.h"
#include "../includes/utils.h"

using namespace std::chrono;

int main(int argc, char **argv) {
  utils::InputInfo input_info;
  utils::ExitCode status;
  std::string input_buffer, clustering_object;
  int exit_code;

  /* Get arguments */
  exit_code = utils::args::ReadArguments(argc, argv, input_info, status);
  switch (exit_code) {
    case utils::SUCCESS:
      std::cout << "\nArguments provided correctly" << std::endl;
      break;
    case utils::FAIL:
      if (status == utils::NO_ARGS) {
        std::cout << "\nNo arguments provided" << std::endl;
        std::cout << "Proceding to input them.." << std::endl;
        exit_code = utils::args::ScanArguments(input_info, status);
        switch (exit_code) {
          case utils::SUCCESS:
            std::cout << "Arguments provided correctly" << std::endl;
            break;
          case utils::FAIL:
            utils::report::ReportError(status);
            break;
          default:
            break;
        }
      } else {
        utils::report::ReportError(status);
      }
      break;
    default:
      break;
  }

  /**
    Read the first line of the input file to define clustering for either
    vectors or curves
  */
  auto start = high_resolution_clock::now();
  std::cout << "\nGetting clustering object.." << std::endl;
  exit_code = utils::io::GetFirstLine(input_info.input_file,
                                      clustering_object, status);
  if (exit_code != utils::SUCCESS) {
    utils::report::ReportError(status);
  }
  auto stop = high_resolution_clock::now();
  duration <double> total_time = duration_cast<duration<double>>(stop - start);
  std::cout << "Getting clustering object completed successfully." << std::endl;
  std::cout << "Time elapsed: " << total_time.count() << " seconds"
            << std::endl;
  std::cout << "\nClustering: " << clustering_object << std::endl;

  /* Read configuration file and pass parameters to the input_info struct */
  start = high_resolution_clock::now();
  std::cout << "\nReading configuration file.." << std::endl;
  exit_code = utils::io::ReadConfig(input_info.config_file, input_info.K,
                                    input_info.grids, input_info.L,
                                    input_info.k, status);
  if (exit_code != utils::SUCCESS) {
    utils::report::ReportError(status);
  }
  stop = high_resolution_clock::now();
  total_time = duration_cast<duration<double>>(stop - start);
  std::cout << "Reading configuration file completed successfully." << std::endl;
  std::cout << "Time elapsed: " << total_time.count() << " seconds"
            << std::endl;
  
  std::srand(std::time(0));

  if (clustering_object == "vectors") {
    #define T double
    #define U std::string
    /**
     Preprocessing input file to get number of dataset vectors
     and their dimension
    */
    start = high_resolution_clock::now();
    std::cout << "\nGetting number of dataset vectors.." << std::endl;
    exit_code = utils::io::vectors::GetNoDataVectors(input_info.input_file,
                                                        input_info.N, status);
    if (exit_code != utils::SUCCESS) {
      utils::report::ReportError(status);
    }
    stop = high_resolution_clock::now();
    total_time = duration_cast<duration<double>>(stop - start);
    std::cout << "Getting number of dataset vectors completed successfully."
              << std::endl;
    std::cout << "Time elapsed: " << total_time.count() << " seconds"
              << std::endl;

    start = high_resolution_clock::now();
    std::cout << "\nGetting dataset vectors' dimension.." << std::endl;
    exit_code = utils::io::vectors::GetVectorsDim(input_info.input_file,
                                                  input_info.D, status);
    if (exit_code != utils::SUCCESS) {
      utils::report::ReportError(status);
    }
    stop = high_resolution_clock::now();
    total_time = duration_cast<duration<double>>(stop - start);
    std::cout << "Getting dataset vectors' dimension completed successfully."
              << std::endl;
    std::cout << "Time elapsed: " << total_time.count() << " seconds"
              << std::endl;

    /*
      Read dataset and create 1D vector which represents the d-dimensional
      vectors of N vectors. Also create 1D vector that stores vectors' ids.
      1D vector of points representation support cache efficiency and as a result
      faster computations
    */
    start = high_resolution_clock::now();
    std::cout << "\nReading dataset.." << std::endl;
    std::vector<T> dataset_vectors(input_info.N * input_info.D);
    std::vector<U> dataset_vectors_ids(input_info.N);
    exit_code = utils::io::vectors::ReadFile<T,U>(input_info.input_file,
                                                  input_info.N, input_info.D,
                                                  dataset_vectors,
                                                  dataset_vectors_ids, status);
    if (exit_code != utils::SUCCESS) {
      utils::report::ReportError(status);
    }
    stop = high_resolution_clock::now();
    total_time = duration_cast<duration<double>>(stop - start);
    std::cout << "Reading dataset completed successfully." << std::endl;
    std::cout << "Time elapsed: " << total_time.count() << " seconds"
              << std::endl;

    /* Print info */
    input_info.Print(clustering_object);

    /**
      Τest space - NOTE whichever test we operate, delete it before commit
      to master branch in order to avoid conflicts
    */
    using namespace cluster::initialization::vectors;
    auto centers = RandomInit(dataset_vectors,input_info.N,input_info.D,input_info.K);
    for(auto& i:centers) {
      std::cout << i << " "; 
    }
    std::cout << std::endl;
    
    using namespace cluster::assignment::vectors;
    auto clusters = LloydsAssignment(dataset_vectors,centers,input_info.N,input_info.D,input_info.K);

    using namespace cluster::update::vectors;
    auto new_centers = LloydsUpdate(dataset_vectors,centers,input_info.N,input_info.D,
                                    input_info.K,std::get<0>(clusters),std::get<1>(clusters));

    for(auto& i:new_centers) {
      std::cout << i << " "; 
    }
    std::cout << std::endl;

  } else if (clustering_object == "curves") {
    #define T double
    #define U int
    /* Preprocessing input file to get number of dataset curves */
    start = high_resolution_clock::now();
    std::cout << "\nGetting number of dataset curves.." << std::endl;
    exit_code = utils::io::curves::GetNoDataCurves(input_info.input_file,
                                                   input_info.N, status);
    if (exit_code != utils::SUCCESS) {
      utils::report::ReportError(status);
    }
    stop = high_resolution_clock::now();
    total_time = duration_cast<duration<double>>(stop - start);
    std::cout << "Getting number of dataset curves completed successfully."
              << std::endl;
    std::cout << "Time elapsed: " << total_time.count() << " seconds"
              << std::endl;

    /*
      Read dataset and create 1D vector of pairs which stores sequentially each
      curve of lenth m_i. Also create 1D vector that stores curves' ids and
      curves' length. 1D vector of curves representation support cache
      efficiency and as a result, faster computations
    */
    start = high_resolution_clock::now();
    std::cout << "\nReading dataset.." << std::endl;
    std::vector<std::pair<T,T>> dataset_curves;
    std::vector<U> dataset_curves_ids(input_info.N);
    std::vector<int> dataset_curves_lengths(input_info.N);
    std::vector<int> dataset_curves_offsets(input_info.N);
    exit_code = utils::io::curves::ReadFile<T,U>(input_info.input_file,
                                        input_info.N, dataset_curves,
                                        dataset_curves_ids,
                                        dataset_curves_lengths,
                                        dataset_curves_offsets, status);
    if (exit_code != utils::SUCCESS) {
      utils::report::ReportError(status);
    }
    stop = high_resolution_clock::now();
    total_time = duration_cast<duration<double>>(stop - start);
    std::cout << "Reading dataset completed successfully." << std::endl;
    std::cout << "Time elapsed: " << total_time.count() << " seconds"
              << std::endl;

    /* Print info */
    input_info.Print(clustering_object);

    /**
      Τest space - NOTE whichever test we operate, delete it before commit
      to master branch in order to avoid conflicts
    */

    /* TEST: Initializing curves centroids */
    using namespace cluster::initialization::curves;
    auto curves_centers = RandomInit(dataset_curves, dataset_curves_lengths,
                                      dataset_curves_offsets, input_info.N,
                                      input_info.K);
    for(auto& i:curves_centers) {
      std::cout << i << " "; 
    }
    std::cout << std::endl;
  
    using namespace cluster::assignment::curves;
    auto clusters = LloydsAssignment(dataset_curves,curves_centers, dataset_curves_lengths,
                                dataset_curves_offsets, input_info.N,
                                input_info.K);
    
    using namespace cluster::update::curves;
    auto new_centers = LloydsUpdate(dataset_curves,curves_centers, dataset_curves_lengths,
                                dataset_curves_offsets, input_info.N,
                                input_info.K, std::get<0>(clusters),std::get<1>(clusters));
    for(auto& i:new_centers) {
      std::cout << i << " "; 
    }
    std::cout << std::endl;
  }

  return EXIT_SUCCESS;
}
