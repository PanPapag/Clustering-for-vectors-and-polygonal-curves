#include <gtest/gtest.h>

#include "../core/utils/utils.h"
#include "../core/metric/metric.h"

TEST(Modulo, ModEq) {
  EXPECT_EQ(0, utils::mod(1,1));
  EXPECT_EQ(1, utils::mod(1,2));
  EXPECT_EQ(2, utils::mod(5,3));
}

TEST(Mean, MeanComputation) {
  std::vector<double> v{1,1,1,1};
  EXPECT_EQ(1, utils::ComputeMean(v,1,4));
  EXPECT_EQ(1, utils::ComputeMean(v,2,2));
}

TEST(Max, FindMax) {
  EXPECT_EQ(4, utils::max(1,2,4));
  EXPECT_EQ(1, utils::max(1,1,1));
  EXPECT_EQ(7, utils::max(1,2,4,3,7,5));
}

TEST(Min, FindMin) {
  EXPECT_EQ(1, utils::min(1,2,4));
  EXPECT_EQ(1, utils::min(1,1,1));
  EXPECT_EQ(1, utils::min(1,2,4,3,7,5));
}

TEST(Distance, Manhattan) {
  std::vector<double> vec{3,5,2,10};
  const std::vector<double>& v = vec;
  EXPECT_EQ(0.0, metric::ManhattanDistance<double> (
  std::next(v.begin(), 0),std::next(v.begin(), 0),
  std::next(v.begin(), 0)));
}

TEST(Distance, Euclidean) {
  std::vector<double> vec{3,5,2,10};
  const std::vector<double>& v = vec;
  EXPECT_EQ(0.0, metric::SquaredEuclidianDistance<double> (
  std::next(v.begin(), 0),std::next(v.begin(), 0),
  std::next(v.begin(), 0)));
}


