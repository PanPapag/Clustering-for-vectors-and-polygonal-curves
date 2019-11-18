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

#include "../../core/metric/metric.h"
#include "../../core/utils/utils.h"

#include "../includes/args_utils.h"
#include "../includes/io_utils.h"
#include "../includes/report_utils.h"
#include "../includes/utils.h"

#define T double
#define U std::string

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

  if (clustering_object == "vectors") {
    /**
     Preprocessing input file to get number of dataset vectors
     and their dimension
    */
    start = high_resolution_clock::now();
    std::cout << "\nGetting number of dataset vectors.." << std::endl;
    exit_code = utils::io::vectors::GetNoDatasetVectors(input_info.input_file,
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
      Read dataset and create 1D vector which represents the d-dimensional points
      of N vectors. Also create 1D vector that stores vectors' ids.
      1D vector of points representation support cache efficiency and as a result
      faster computations
    */
    start = high_resolution_clock::now();
    std::cout << "\nReading dataset.." << std::endl;
    std::vector<T> dataset_vectors;
    dataset_vectors.reserve(input_info.N * input_info.D);
    std::vector<U> dataset_ids;
    dataset_ids.reserve(input_info.N);
    exit_code = utils::io::vectors::ReadFile<T,U>(input_info.input_file,
                                                  input_info.N, input_info.D,
                                                  dataset_vectors, dataset_ids,
                                                  status);
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
  } else if (clustering_object == "curves") {
    std::cout << "curves" << std::endl;
  }

  return EXIT_SUCCESS;
}
