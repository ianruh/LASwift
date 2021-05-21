// VectorTrigonometry.swift
//
// Copyright (c) 2017 Alexander Taraymovich <taraymovich@me.com>
// All rights reserved.
//
// This software may be modified and distributed under the terms
// of the BSD license. See the LICENSE file for details.

import CLAPACK
#if os(macOS) || os(iOS) || os(tvOS) || os(watchOS)
import Accelerate
#endif

// MARK: - Trigonometric operations on vector

/// Return the sine of `a`, where `a` is given in radians and the return value is
/// in the range -1 to 1.
///
/// - Parameters
///     - a: vector
/// - Returns: sine of a vector values
public func sin(_ a: Vector) -> Vector {
    #if os(macOS) || os(iOS) || os(tvOS) || os(watchOS)
    return vectorFunction(vvsin, a)
    #else
    return a.map({Double.sin($0)})
    #endif
}

/// Return the cosine of `a`, where `a` is given in radians and the return value is
/// in the range -1 to 1.
///
/// - Parameters
///     - a: vector
/// - Returns: cosine of a vector values
public func cos(_ a: Vector) -> Vector {
    #if os(macOS) || os(iOS) || os(tvOS) || os(watchOS)
    return vectorFunction(vvcos, a)
    #else
    return a.map({Double.cos($0)})
    #endif
}

/// Return the tangent of `a`, where `a` is given in radians.
///
/// Mathematically, the tangent function has singularities at odd multiples of
/// pi/2. If the argument x is too close to one of these singularities, tan
/// will return extremely large value.
///
/// - Parameters
///     - a: vector
/// - Returns: tangent of a vector values
public func tan(_ a: Vector) -> Vector {
    #if os(macOS) || os(iOS) || os(tvOS) || os(watchOS)
    return vectorFunction(vvtan, a)
    #else
    return a.map({Double.tan($0)})
    #endif
}
