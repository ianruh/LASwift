// VectorStatistics.swift
//
// Copyright (c) 2017 Alexander Taraymovich <taraymovich@me.com>
// All rights reserved.
//
// This software may be modified and distributed under the terms
// of the BSD license. See the LICENSE file for details.

import CLAPACK

// MARK: - Statistical operations on vector

/// Return the largest element of vector.
///
/// - Parameters
///     - a: vector
/// - Returns: largest element of vector a
public func max(_ a: Vector) -> Double {
    return a.max()!
}

/// Return the index of largest element of vector.
///
/// - Parameters
///     - a: vector
/// - Returns: index of largest element of vector a
public func maxi(_ a: Vector) -> Int {
    return a.firstIndex(of: a.max()!)!
}

/// Return the smallest element of vector.
///
/// - Parameters
///     - a: vector
/// - Returns: smallest element of vector a
public func min(_ a: Vector) -> Double {
    return a.min()!
}

/// Return the index of smallest element of vector.
///
/// - Parameters
///     - a: vector
/// - Returns: index of smallest element of vector a
public func mini(_ a: Vector) -> Int {
    return a.firstIndex(of: a.min()!)!
}

/// Return mean (statistically average) value of vector.
///
/// - Parameters
///     - a: vector
/// - Returns: mean value of vector a
public func mean(_ a: Vector) -> Double {
    return sum(a)/Double(a.count)
}

/// Return standard deviation value of vector.
///
/// - Parameters
///     - a: vector
/// - Returns: standard deviation value of vector a
public func std(_ a: Vector) -> Double {
    let μ = mean(a)
    return Double.sqrt(sumsq(a - μ)/Double(a.count))
}

/// Return normalized vector (substract mean value and divide by standard deviation).
///
/// - Parameters
///     - a: vector
/// - Returns: normalized vector a
public func normalize(_ a: Vector) -> Vector {
    let dev = std(a)
    let avg = mean(a)
    return a.map({($0-avg) / dev})
}

/// Return sum of vector's elements.
///
/// - Parameters
///     - a: vector
/// - Returns: sum of elements of vector a
public func sum(_ a: Vector) -> Double {
    return a.reduce(0, +)
}

/// Return sum of vector's squared elements.
///
/// - Parameters
///     - a: vector
/// - Returns: sum of squared elements of vector a
public func sumsq(_ a: Vector) -> Double {
    let sq = square(a)
    return sq.reduce(0, +)
}
