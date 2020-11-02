// MatrixStatistics.swift
//
// Copyright (c) 2017 Alexander Taraymovich <taraymovich@me.com>
// All rights reserved.
//
// This software may be modified and distributed under the terms
// of the BSD license. See the LICENSE file for details.

import CLAPACK

// MARK: - Statistic functions on matrix

/// Return vector of largest elements of matrix in a specified dimension.
///
/// - Parameters
///     - A: matrix
///     - d: dimension (.Column or .Row, defaults to .Row)
public func max(_ A: Matrix, _ d: Dim = .Row) -> Vector {
    var maxes: [Double] = []
    switch d {
    case .Row:
        for r in 0..<A.rows {
            maxes.append(A.flat[r*A.cols..<(r+1)*A.cols].max()!)
        }
    case .Column:
        let B = transpose(A)
        for r in 0..<B.rows {
            maxes.append(B.flat[r*B.cols..<(r+1)*B.cols].max()!)
        }
    }
    return maxes
}

/// Return vector of indices of largest elements of matrix in a specified dimension.
///
/// - Parameters
///     - A: matrix
///     - d: dimension (.Column or .Row, defaults to .Row)
public func maxi(_ A: Matrix, _ d: Dim = .Row) -> [Int] {
    var maxes: [Int] = []
    switch d {
    case .Row:
        for r in 0..<A.rows {
            let slice = A.flat[r*A.cols..<(r+1)*A.cols]
            maxes.append(slice.firstIndex(of: slice.max()!)! - r*A.cols)
        }
    case .Column:
        let B = transpose(A)
        for r in 0..<B.rows {
            let slice = B.flat[r*B.cols..<(r+1)*B.cols]
            maxes.append(slice.firstIndex(of: slice.max()!)! - r*B.cols)
        }
    }
    return maxes
}

/// Return vector of smallest elements of matrix in a specified dimension.
///
/// - Parameters
///     - A: matrix
///     - d: dimension (.Column or .Row, defaults to .Row)
public func min(_ A: Matrix, _ d: Dim = .Row) -> Vector {
    var mins: [Double] = []
    switch d {
    case .Row:
        for r in 0..<A.rows {
            mins.append(A.flat[r*A.cols..<(r+1)*A.cols].min()!)
        }
    case .Column:
        let B = transpose(A)
        for r in 0..<B.rows {
            mins.append(B.flat[r*B.cols..<(r+1)*B.cols].min()!)
        }
    }
    return mins
}

/// Return vector of indices of smallest elements of matrix in a specified dimension.
///
/// - Parameters
///     - A: matrix
///     - d: dimension (.Column or .Row, defaults to .Row)
public func mini(_ A: Matrix, _ d: Dim = .Row) -> [Int] {
    var mins: [Int] = []
    switch d {
    case .Row:
        for r in 0..<A.rows {
            let slice = A.flat[r*A.cols..<(r+1)*A.cols]
            mins.append(slice.firstIndex(of: slice.min()!)! - r*A.cols)
        }
    case .Column:
        let B = transpose(A)
        for r in 0..<B.rows {
            let slice = B.flat[r*B.cols..<(r+1)*B.cols]
            mins.append(slice.firstIndex(of: slice.min()!)! - r*B.cols)
        }
    }
    return mins
}

/// Return vector of mean values of matrix in a specified dimension.
///
/// - Parameters
///     - A: matrix
///     - d: dimension (.Column or .Row, defaults to .Row)
public func mean(_ A: Matrix, _ d: Dim = .Row) -> Vector {
    var means: [Double] = []
    switch d {
    case .Row:
        for r in 0..<A.rows {
            let slice = A.flat[r*A.cols..<(r+1)*A.cols]
            means.append(slice.reduce(0, +) / Double(slice.count))
        }
    case .Column:
        let B = transpose(A)
        for r in 0..<B.rows {
            let slice = B.flat[r*B.cols..<(r+1)*B.cols]
            means.append(slice.reduce(0, +) / Double(slice.count))
        }
    }
    return means
}

/// Return vector of standard deviation values of matrix in a specified dimension.
///
/// - Parameters
///     - A: matrix
///     - d: dimension (.Column or .Row, defaults to .Row)
public func std(_ A: Matrix, _ d: Dim = .Row) -> Vector {
    var stds: [Double] = []
    switch d {
    case .Row:
        for r in 0..<A.rows {
            let slice = Array(A.flat[r*A.cols..<(r+1)*A.cols])
            let μ = slice.reduce(0, +) / Double(slice.count)
            stds.append(Double.sqrt(sumsq(slice - μ)/Double(slice.count)))
        }
    case .Column:
        let B = transpose(A)
        for r in 0..<B.rows {
            let slice = Array(B.flat[r*B.cols..<(r+1)*B.cols])
            let μ = slice.reduce(0, +) / Double(slice.count)
            stds.append(Double.sqrt(sumsq(slice - μ)/Double(slice.count)))
        }
    }
    return stds
}

/// Return normalized matrix (substract mean value and divide by standard deviation)
/// in dpecified dimension.
///
/// - Parameters
///     - A: matrix
///     - d: dimension (.Column or .Row, defaults to .Row)
/// - Returns: normalized matrix A
public func normalize(_ A: Matrix, _ d: Dim = .Row) -> Matrix {
    switch d {
    case .Row:
        return Matrix(A.map { normalize(Vector($0)) })
    case .Column:
        return Matrix(transpose(A).map { normalize(Vector($0)) })
    }
}

/// Return vector of sums of values of matrix in a specified dimension.
///
/// - Parameters
///     - A: matrix
///     - d: dimension (.Column or .Row, defaults to .Row)
public func sum(_ A: Matrix, _ d: Dim = .Row) -> Vector {
    var sums: [Double] = []
    switch d {
    case .Row:
        for r in 0..<A.rows {
            let slice = A.flat[r*A.cols..<(r+1)*A.cols]
            sums.append(slice.reduce(0, +))
        }
    case .Column:
        let B = transpose(A)
        for r in 0..<B.rows {
            let slice = B.flat[r*B.cols..<(r+1)*B.cols]
            sums.append(slice.reduce(0, +))
        }
    }
    return sums
}

/// Return vector of sums of squared values of matrix in a specified dimension.
///
/// - Parameters
///     - A: matrix
///     - d: dimension (.Column or .Row, defaults to .Row)
public func sumsq(_ A: Matrix, _ d: Dim = .Row) -> Vector {
    var sums: [Double] = []
    switch d {
    case .Row:
        for r in 0..<A.rows {
            let slice = Array(A.flat[r*A.cols..<(r+1)*A.cols])
            sums.append(square(slice).reduce(0, +))
        }
    case .Column:
        let B = transpose(A)
        for r in 0..<B.rows {
            let slice = Array(B.flat[r*B.cols..<(r+1)*B.cols])
            sums.append(square(slice).reduce(0, +))
        }
    }
    return sums
}
