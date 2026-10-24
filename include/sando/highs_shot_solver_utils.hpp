/* ----------------------------------------------------------------------------
 * Copyright 2026, Kota Kondo, Aerospace Controls Laboratory
 * Massachusetts Institute of Technology
 * All Rights Reserved
 * Authors: Kota Kondo, et al.
 * See LICENSE file for the license information
 * -------------------------------------------------------------------------- */

// Counterpart of gurobi_solver_utils.hpp for the HiGHS/SHOT modeling layer. The helpers live in
// namespace hs so that this header can be included alongside gurobi_solver_utils.hpp.

#pragma once

#include <algorithm>
#include <cassert>
#include <functional>
#include <iterator>
#include <vector>
#include <Eigen/Dense>
#include <sando/highs_shot_model.hpp>

namespace hs {

/** @brief Find the smallest positive element in a vector.
 *  @param v Vector of doubles to search.
 *  @return The minimum positive value, or 0 if none exists.
 */
inline double MinPositiveElement(std::vector<double> v) {
  std::sort(v.begin(), v.end());  // sorted in ascending order
  double min_value = 0;
  for (size_t i = 0; i < v.size(); i++) {
    if (v[i] > 0) {
      min_value = v[i];
      break;
    }
  }
  return min_value;
}

/** @brief Compute the squared L2 norm of a vector as a quadratic expression.
 *  @tparam T Element type (Var or LinExpr).
 *  @param x Input vector.
 *  @return Quadratic expression representing sum of squared elements.
 */
template <typename T>
QuadExpr GetNorm2(const std::vector<T>& x) {
  QuadExpr result = 0;
  for (size_t i = 0; i < x.size(); i++) {
    result = result + LinExpr(x[i]) * LinExpr(x[i]);
  }
  return result;
}

/** @brief Multiply a dense matrix by a vector of expressions.
 *  @tparam T Element type of the vector (Var or LinExpr).
 *  @param A Dense matrix stored as vector of rows.
 *  @param x Input vector.
 *  @return Result vector of linear expressions A*x.
 */
template <typename T>
std::vector<LinExpr> MatrixMultiply(
    const std::vector<std::vector<double>>& A, const std::vector<T>& x) {
  std::vector<LinExpr> result;

  for (size_t i = 0; i < A.size(); i++) {
    LinExpr lin_exp = 0;
    for (size_t m = 0; m < x.size(); m++) {
      lin_exp += A[i][m] * LinExpr(x[m]);
    }
    result.push_back(lin_exp);
  }
  return result;
}

/** @brief Element-wise addition of two std::vectors.
 *  @tparam T Element type.
 *  @param a First vector.
 *  @param b Second vector (must have same size as a).
 *  @return Vector containing element-wise sums.
 */
template <typename T>
std::vector<T> operator+(const std::vector<T>& a, const std::vector<T>& b) {
  assert(a.size() == b.size());

  std::vector<T> result;
  result.reserve(a.size());

  std::transform(a.begin(), a.end(), b.begin(), std::back_inserter(result), std::plus<T>());
  return result;
}

/** @brief Element-wise subtraction of two std::vectors.
 *  @tparam T Element type.
 *  @param a First vector.
 *  @param b Second vector (must have same size as a).
 *  @return Vector containing element-wise differences.
 */
template <typename T>
std::vector<T> operator-(const std::vector<T>& a, const std::vector<T>& b) {
  assert(a.size() == b.size());

  std::vector<T> result;
  result.reserve(a.size());

  std::transform(a.begin(), a.end(), b.begin(), std::back_inserter(result), std::minus<T>());
  return result;
}

/** @brief Subtract a double vector from an expression vector element-wise.
 *  @tparam T Expression type (Var or LinExpr).
 *  @param x Expression vector.
 *  @param b Double vector to subtract.
 *  @return Vector of linear expressions representing x - b.
 */
template <typename T>
std::vector<LinExpr> operator-(const std::vector<T>& x, const std::vector<double>& b) {
  std::vector<LinExpr> result;
  for (size_t i = 0; i < x.size(); i++) {
    LinExpr tmp = LinExpr(x[i]) - b[i];
    result.push_back(tmp);
  }
  return result;
}

/** @brief Convert a dynamic-size Eigen column vector to a std::vector.
 *  @tparam T Scalar type.
 *  @param x Eigen column vector.
 *  @return Equivalent std::vector.
 */
template <typename T>
std::vector<T> eigenVector2std(const Eigen::Matrix<T, -1, 1>& x) {
  std::vector<T> result;
  for (int i = 0; i < x.rows(); i++) {
    result.push_back(x(i));
  }
  return result;
}

/** @brief Convert a dynamic-size Eigen matrix to a nested std::vector.
 *  @tparam T Scalar type.
 *  @param x Eigen matrix.
 *  @return Nested vector where each inner vector is a row.
 */
template <typename T>
std::vector<std::vector<T>> eigenMatrix2std(const Eigen::Matrix<T, -1, -1>& x) {
  std::vector<std::vector<T>> result;

  for (int i = 0; i < x.rows(); i++) {
    std::vector<T> row;
    for (int j = 0; j < x.cols(); j++) {
      row.push_back(x(i, j));
    }
    result.push_back(row);
  }
  return result;
}

/** @brief Convert a fixed-size Eigen matrix to a nested std::vector.
 *  @tparam T Scalar type.
 *  @tparam Rows Number of rows.
 *  @tparam Cols Number of columns.
 *  @param x Eigen matrix.
 *  @return Nested vector where each inner vector is a row.
 */
template <typename T, int Rows, int Cols>
std::vector<std::vector<T>> eigenMatrix2std(const Eigen::Matrix<T, Rows, Cols>& x) {
  std::vector<std::vector<T>> result;

  for (int i = 0; i < x.rows(); i++) {
    std::vector<T> row;
    for (int j = 0; j < x.cols(); j++) {
      row.push_back(x(i, j));
    }
    result.push_back(row);
  }
  return result;
}

/** @brief Extract a single column from a 2D nested std::vector.
 *  @tparam T Element type.
 *  @param x 2D nested vector (matrix).
 *  @param column Column index to extract.
 *  @return Vector containing the specified column.
 */
template <typename T>
std::vector<T> GetColumn(std::vector<std::vector<T>> x, int column) {
  std::vector<T> result;

  for (size_t i = 0; i < x.size(); i++) {
    result.push_back(x[i][column]);
  }
  return result;
}

}  // namespace hs
