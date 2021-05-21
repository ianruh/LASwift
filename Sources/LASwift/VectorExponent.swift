// VectorExponent.swift
//
// Copyright (c) 2017 Alexander Taraymovich <taraymovich@me.com>
// All rights reserved.
//
// This software may be modified and distributed under the terms
// of the BSD license. See the LICENSE file for details.

import CLAPACK
import Numerics
#if os(macOS) || os(iOS) || os(tvOS) || os(watchOS)
import Accelerate
#endif

// MARK: - Power and exponential operations on vector

/// Exponentiation function, returning vector raised to power.
///
/// Alternatively, `power(a, p)` can be executed with `a .^ p`.
///
/// Mathematically, `power` would return a complex number when base is negative and
/// power is not an integral value. `power` can’t do that,
/// so instead it signals domain error (returns `±NaN`).
///
/// - Parameters
///     - a: vector
///     - p: power to raise vector to
/// - Returns: elementwise vector power of a raised to p
public func power(_ a: Vector, _ p: Double) -> Vector {
    #if os(macOS) || os(iOS) || os(tvOS) || os(watchOS)
    var c = Vector(repeating: 0.0, count: a.count)
    var l = Int32(a.count)
    var p = p
    vvpows(&c, &p, a, &l)
    return c
    #else
    return a.map({.pow($0, p)})
    #endif
}

/// Exponentiation function, returning vector raised to power.
///
/// Alternatively, `a .^ p` can be executed with `power(a, p)`.
///
/// Mathematically, `power` would return a complex number when base is negative and
/// power is not an integral value. `power` can’t do that, 
/// so instead it signals domain error (returns `±NaN`).
///
/// - Parameters
///     - a: vector
///     - p: power to raise vector to
/// - Returns: elementwise vector power of a raised to p
public func .^ (_ a: Vector, _ p: Double) -> Vector {
    return power(a, p)
}

/// Exponentiation function, returning vector raised to power of 2.
///
/// - Parameters
///     - a: vector
/// - Returns: elementwise vector power of a raised to power of 2
public func square(_ a: Vector) -> Vector {
    #if os(macOS) || os(iOS) || os(tvOS) || os(watchOS)
    return unaryVectorOperation(vDSP_vsqD, a)
    #else
    return a.map({.pow($0, 2)})
    #endif
}

/// Exponentiation function, returning square root of vector.
///
/// Mathematically, `sqrt` would return a complex number when base is negative.
/// `sqrt` can’t do that, so instead it signals domain error (returns `±NaN`).
///
/// - Parameters
///     - a: vector
/// - Returns: elementwise square root of vector a
public func sqrt(_ a: Vector) -> Vector {
    #if os(macOS) || os(iOS) || os(tvOS) || os(watchOS)
    return vectorFunction(vvsqrt, a)
    #else
    return a.map({.sqrt($0)})
    #endif
}

///  Compute `e` (the base of natural logarithms) raised to the power `a`.
///
///  If the magnitude of the result is too large to be representable, exp
///  signals overflow (returns `Inf`).
///
/// - Parameters
///     - a: vector
/// - Returns: elementwise `e` raised to the power of vector a
public func exp(_ a: Vector) -> Vector {
    #if os(macOS) || os(iOS) || os(tvOS) || os(watchOS)
    return vectorFunction(vvexp, a)
    #else
    return a.map({.exp($0)})
    #endif
}

/// Compute the natural logarithm of `a` where `exp(log(a))` equals `a`, exactly in
/// mathematics and approximately in C.
///
/// If x is negative, log signals a domain error (returns `NaN`). If x is zero, it returns
/// negative infinity (`-Inf`); if x is too close to zero, it may signal overflow.
///
/// - Parameters
///     - a: vector
/// - Returns: elementwise natural logarithm of vector a
public func log(_ a: Vector) -> Vector {
    #if os(macOS) || os(iOS) || os(tvOS) || os(watchOS)
    return vectorFunction(vvlog, a)
    #else
    let loge = Double.log(Double.exp(1.0))
    return a.map({.log($0) / loge})
    #endif
}

/// Return the base-2 logarithm of `a`, where `log2(a) = log(a)/log(2)`.
///
/// - Parameters
///     - a: vector
/// - Returns: elementwise base-2 logarithm of vector a
public func log2(_ a: Vector) -> Vector {
    #if os(macOS) || os(iOS) || os(tvOS) || os(watchOS)
    return vectorFunction(vvlog2, a)
    #else
    let log2 = Double.log(2.0)
    return a.map({.log($0) / log2})
    #endif
}

/// Return the base-10 logarithm of `a`, where `log10(a) = log(a)/log(10)`.
///
/// - Parameters
///     - a: vector
/// - Returns: elementwise base-10 logarithm of vector a
public func log10(_ a: Vector) -> Vector {
    #if os(macOS) || os(iOS) || os(tvOS) || os(watchOS)
    return vectorFunction(vvlog10, a)
    #else
    let log10 = Double.log(10)
    return a.map({.log($0) / log10})
    #endif
}
