// MatrixAlgebra.swift
//
// Copyright (c) 2017 Alexander Taraymovich <taraymovich@me.com>
// All rights reserved.
//
// This software may be modified and distributed under the terms
// of the BSD license. See the LICENSE file for details.

import CLAPACK
import CBLAS
#if os(macOS) || os(iOS) || os(tvOS) || os(watchOS)
import Accelerate
#endif

/// Matrix triangular part.
///
/// - Upper: Upper triangular part.
/// - Lower: Lower triangular part.
public enum Triangle {
    case Upper
    case Lower
}

// MARK: - Linear algebra operations on matrices

/// Compute the trace of a matrix.
///
/// - Parameters:
///     - A: matrix to compute trace of
/// - Returns: sum of the elements on the main diagonal
public func trace(_ A: Matrix) -> Double {
    precondition(A.rows == A.cols, "Matrix dimensions must agree")
    return sum((0..<A.rows).map { A[$0, $0] })
}

/// Transpose matrix.
///
/// Alternatively, `transpose(A)` can be executed with `A′`.
///
/// - Parameters
///     - A: matrix
/// - Returns: transposed matrix
public func transpose(_ A: Matrix) -> Matrix {
    let C: Matrix = zeros(A.cols, A.rows)

    #if os(macOS) || os(iOS) || os(tvOS) || os(watchOS)
    vDSP_mtransD(A.flat, 1, &(C.flat), 1, vDSP_Length(A.cols), vDSP_Length(A.rows))
    #else
    for n in 0..<A.flat.count {
        let i: Int = n/A.rows
        let j: Int = n%A.rows
        C.flat[n] = A.flat[A.cols*j + i]
    }
    #endif

    return C
}



/// Transpose matrix.
///
/// Alternatively, `A′` can be executed with `transpose(A)`.
///
/// - Parameters
///     - A: matrix
/// - Returns: transposed matrix
public postfix func ′ (_ a: Matrix) -> Matrix {
    return transpose(a)
}

/// Perform matrix multiplication.
///
/// Alternatively, `mtimes(A, B)` can be executed with `A * B`.
///
/// - Parameters
///     - A: left matrix
///     - B: right matrix
/// - Returns: matrix product of A and B
public func mtimes(_ A: Matrix, _ B: Matrix) -> Matrix {
    precondition(A.cols == B.rows, "Matrix dimensions must agree")

    // CBLAS enum values
    let cblasNoTranspose = CBLAS_TRANSPOSE.init(rawValue: 111)
//    let cblasTranspose = CBLAS_TRANSPOSE.init(rawValue: 112)

    /* Compute matrix product*/

    let A_copy = Matrix(A)
    let B_copy = Matrix(B)

    let order = CBLAS_ORDER.init(101) // Row major ordering
    let M = __CLPK_integer(A_copy.rows)
    let N = __CLPK_integer(B_copy.cols)
    let K = __CLPK_integer(A_copy.cols)
    let alpha: Double = 1.0
    let beta: Double = 0.0
    let LDA = __CLPK_integer(A_copy.cols)
    let LDB = __CLPK_integer(B_copy.cols)

    let C = zeros(A_copy.rows, B_copy.cols)
    let LDC = __CLPK_integer(C.cols)

    let transposeA = cblasNoTranspose
    let transposeB = cblasNoTranspose

    cblas_dgemm(order, transposeA, transposeB, M, N, K, alpha, &A_copy.flat, LDA, &B_copy.flat, LDB, beta, &C.flat, LDC)

    return C
}

/// Perform matrix multiplication.
///
/// Alternatively, `A * B` can be executed with `mtimes(A, B)`.
///
/// - Parameters
///     - A: left matrix
///     - B: right matrix
/// - Returns: matrix product of A and B
public func * (_ A: Matrix, _ B: Matrix) -> Matrix {
    return mtimes(A, B)
}

/// Raise matrix to specified power (integer value).
///
/// If power is 1, source matrix is returned
/// If power is -1, inverted matrix is returned
/// If power is > 1, continuous matrix product result is returned (eg `mpower(A, 2) = mtimes(A, A)`)
/// If power is < -1, continuous matrix product of inverted matrix result is returned
/// All other values are invalid
///
/// Alternatively, `mpower(A, p)` can be executed with `A ^ p`.
///
/// - Parameters
///     - A: matrix
///     - p: power to raise matrix to (integer)
/// - Returns: matrix A raised to power p
public func mpower(_ A: Matrix, _ p: Int) -> Matrix {
    precondition(A.cols == A.rows, "Matrix dimensions must agree")
    switch p {
    case 1:
        return Matrix(A)
    case -1:
        return inv(A)
    case _ where p > 1:
        var C = Matrix(A)
        var p = p
        while (p > 1) {
            C = mtimes(A, C)
            p -= 1
        }
        return C
    case _ where p < -1:
        return inv(mpower(A, -p))
    default:
        return eye(A.rows, A.cols)
    }
}

/// Raise matrix to specified power (integer value).
///
/// If power is 1, source matrix is returned
/// If power is -1, inverted matrix is returned
/// If power is > 1, continuous matrix product result is returned (eg `A ^ 2 = mtimes(A, A)`)
/// If power is < -1, continuous matrix product of inverted matrix result is returned
/// All other values are invalid
///
/// Alternatively, `A ^ p` can be executed with `mpower(A, p)`.
///
/// - Parameters
///     - A: matrix
///     - p: power to raise matrix to (integer)
/// - Returns: matrix A raised to power p
public func ^ (_ a: Matrix, _ p: Int) -> Matrix {
    return mpower(a, p)
}

/// Compute the inverse of a given square matrix.
///
///	A precondition error is thrown if the given matrix is singular or algorithm fails to converge.
///
/// - Parameters:
///     - A: square matrix to invert
/// - Returns: inverse of A matrix
public func inv(_ A: Matrix) -> Matrix {
    precondition(A.rows == A.cols, "Matrix dimensions must agree")
    let B = Matrix(A)
    
    var M = __CLPK_integer(A.rows)
    var N = M
    var LDA = N
    var pivot = [__CLPK_integer](repeating: 0, count: Int(N))
    
    var wkOpt = __CLPK_doublereal(0.0)
    var lWork = __CLPK_integer(-1)
    
    var error: __CLPK_integer = 0
    
    dgetrf_(&M, &N, &(B.flat), &LDA, &pivot, &error)
    
    precondition(error == 0, "Matrix is non invertible")
    
    /* Query and allocate the optimal workspace */
    
    dgetri_(&N, &(B.flat), &LDA, &pivot, &wkOpt, &lWork, &error)
    
    lWork = __CLPK_integer(wkOpt)
    var work = Vector(repeating: 0.0, count: Int(lWork))
    
    /* Compute inversed matrix */
    
    dgetri_(&N, &(B.flat), &LDA, &pivot, &work, &lWork, &error)
    
    precondition(error == 0, "Matrix is non invertible")
    
    return B
}

/// Compute eigen values and vectors of a given square matrix.
///
///	A precondition error is thrown if the algorithm fails to converge.
///
/// - Parameters:
///     - A: square matrix to calculate eigen values and vectors of
/// - Returns: eigenvectors matrix (by rows) and diagonal matrix with eigenvalues on the main diagonal
public func eig(_ A: Matrix) -> (V: Matrix, D: Matrix) {
    precondition(A.rows == A.cols, "Matrix dimensions must agree")
    
    let V = Matrix(A)
    
    var N = __CLPK_integer(A.rows)
    
    var LDA = N
    
    var wkOpt = __CLPK_doublereal(0.0)
    var lWork = __CLPK_integer(-1)
    
    var jobvl: Int8 = 86 // 'V'
    var jobvr: Int8 = 86 // 'V'
    
    var error = __CLPK_integer(0)
    
    // Real parts of eigenvalues
    var wr = Vector(repeating: 0.0, count: Int(N))
    // Imaginary parts of eigenvalues
    var wi = Vector(repeating: 0.0, count: Int(N))
    // Left eigenvectors
    var vl = [__CLPK_doublereal](repeating: 0.0, count: Int(N * N))
    // Right eigenvectors
    var vr = [__CLPK_doublereal](repeating: 0.0, count: Int(N * N))
    
    var ldvl = N
    var ldvr = N
    
    /* Query and allocate the optimal workspace */
    
    dgeev_(&jobvl, &jobvr, &N, &V.flat, &LDA, &wr, &wi, &vl, &ldvl, &vr, &ldvr, &wkOpt, &lWork, &error)
    
    lWork = __CLPK_integer(wkOpt)
    var work = Vector(repeating: 0.0, count: Int(lWork))
    
    /* Compute eigen vectors */
    
    dgeev_(&jobvl, &jobvr, &N, &V.flat, &LDA, &wr, &wi, &vl, &ldvl, &vr, &ldvr, &work, &lWork, &error)
    
    precondition(error == 0, "Failed to compute eigen vectors")
    
    return (toRows(Matrix(A.rows, A.cols, vl), .Column), diag(wr))
}

/// Perform a singular value decomposition of a given matrix.
///
/// - Parameters:
///    - A: matrix to find singular values of
/// - Returns: matrices U, S, and V such that `A = U * S * transpose(V)`
public func svd(_ A: Matrix) -> (U: Matrix, S: Matrix, V: Matrix) {
    /* LAPACK is using column-major order */
    let _A = toCols(A, .Row)
    
    var jobz: Int8 = 65 // 'A'
    
    var M = __CLPK_integer(A.rows);
    var N = __CLPK_integer(A.cols);
    
    var LDA = M;
    var LDU = M;
    var LDVT = N;
    
    var wkOpt = __CLPK_doublereal(0.0)
    var lWork = __CLPK_integer(-1)
    var iWork = [__CLPK_integer](repeating: 0, count: Int(8 * min(M, N)))
    
    var error = __CLPK_integer(0)
    
    var s = Vector(repeating: 0.0, count: Int(min(M, N)))
    let U = Matrix(Int(LDU), Int(M))
    let VT = Matrix(Int(LDVT), Int(N))
    
    /* Query and allocate the optimal workspace */
    dgesdd_(&jobz, &M, &N, &_A.flat, &LDA, &s, &U.flat, &LDU, &VT.flat, &LDVT, &wkOpt, &lWork, &iWork, &error)
    
    lWork = __CLPK_integer(wkOpt)
    var work = Vector(repeating: 0.0, count: Int(lWork))
    
    /* Compute SVD */
    dgesdd_(&jobz, &M, &N, &_A.flat, &LDA, &s, &U.flat, &LDU, &VT.flat, &LDVT, &work, &lWork, &iWork, &error)
    
    precondition(error == 0, "Failed to compute SVD")
    
    return (toRows(U, .Column), diag(Int(M), Int(N), s), VT)
}

/// Perform a generalized singular value decomposition of 2 given matrices.
///
/// - Parameters:
///    - A: first matrix
///    - B: second matrix
/// - Returns: matrices U, V, and Q, plus vectors alpha and beta
public func gsvd(_ A: Matrix, _ B: Matrix) -> (U: Matrix, V: Matrix, Q: Matrix, alpha: Vector, beta: Vector, success: Bool) {
    /* LAPACK is using column-major order */
    let _A = toCols(A, .Row)
    let _B = toCols(B, .Row)
    
    var jobu:Int8 = Int8(Array("U".utf8).first!)
    var jobv:Int8 = Int8(Array("V".utf8).first!)
    var jobq:Int8 = Int8(Array("Q".utf8).first!)
    
    var M = __CLPK_integer(A.rows)
    var N = __CLPK_integer(A.cols)
    var P = __CLPK_integer(B.rows)
    
    var LDA = M
    var LDB = P
    var LDU = M
    var LDV = P
    var LDQ = N
    
    var k = __CLPK_integer()
    var l = __CLPK_integer()
    
    let U = Matrix(Int(LDU), Int(M))
    let V = Matrix(Int(LDV), Int(P))
    let Q = Matrix(Int(LDQ), Int(N))
    var alpha = Vector(repeating: 0.0, count: Int(N))
    var beta = Vector(repeating: 0.0, count: Int(N))

    // dggsvd_ is depricated by BLAS, but is still included in Accelerate. However, dggsvd3_ is not in Accelrate.
    #if os(macOS) || os(iOS) || os(tvOS) || os(watchOS)
    var lWork = __CLPK_integer(max(max(Int(3*N),Int(M)),Int(P))+Int(N))
    var iWork = [__CLPK_integer](repeating: 0, count: Int(N))
    var work = Vector(repeating: 0.0, count: Int(lWork) * 4)
    var error = __CLPK_integer(0)

    dggsvd_(&jobu, &jobv, &jobq, &M, &N, &P, &k, &l, &_A.flat, &LDA, &_B.flat, &LDB, &alpha, &beta, &U.flat, &LDU, &V.flat, &LDV, &Q.flat, &LDQ, &work, &iWork, &error)
    #else
    var lWork = __CLPK_integer(-1)
    var iWork = [__CLPK_integer](repeating: 0, count: Int(N))
    var wkOpt = __CLPK_doublereal(0.0)
    var error = __CLPK_integer(0)

    dggsvd3_(&jobu, &jobv, &jobq, &M, &N, &P, &k, &l, &_A.flat, &LDA, &_B.flat, &LDB, &alpha, &beta, &U.flat, &LDU, &V.flat, &LDV, &Q.flat, &LDQ, &wkOpt, &lWork, &iWork, &error)

    lWork = __CLPK_integer(wkOpt)
    var work = Vector(repeating: 0.0, count: Int(lWork) * 4)

    dggsvd3_(&jobu, &jobv, &jobq, &M, &N, &P, &k, &l, &_A.flat, &LDA, &_B.flat, &LDB, &alpha, &beta, &U.flat, &LDU, &V.flat, &LDV, &Q.flat, &LDQ, &work, &lWork, &iWork, &error)

    #endif
    
    //precondition(error == 0, "Failed to compute SVD")
    if error != 0 {
        return (toRows(U, .Column), toRows(V, .Column), toRows(Q, .Column), Vector(alpha[Int(k)...Int(k+l)-1]), Vector(beta[Int(k)...Int(k+l)-1]), false)
    } else {
        return (toRows(U, .Column), toRows(V, .Column), toRows(Q, .Column), Vector(alpha[Int(k)...Int(k+l)-1]), Vector(beta[Int(k)...Int(k+l)-1]), true)
    }
}

/// Compute the Cholesky factorization of a real symmetric positive definite matrix.
///
///	A precondition error is thrown if the algorithm fails to converge.
///
/// - Parameters:
///     - A: square matrix to compute Cholesky factorization of
///     - t: Triangle value (.Upper, .Lower)
/// - Returns: upper triangular matrix U so that `A = U' * U` or 
///            lower triangular matrix L so that `A = L * L'`
public func chol(_ A: Matrix, _ t: Triangle = .Upper) -> Matrix {
    precondition(A.rows == A.cols, "Matrix dimensions must agree")
    
    var uplo: Int8
    switch t {
    case .Upper:
        uplo = 85 // 'U'
    case .Lower:
        uplo = 76 // 'L'
    }
    
    var N = __CLPK_integer(A.rows)
    
    /* LAPACK is using column-major order */
    var U = toCols(A, .Row)
    
    var LDA = N
    
    var error = __CLPK_integer(0)
    
    /* Compute Cholesky decomposition */
    
    dpotrf_(&uplo, &N, &U.flat, &LDA, &error)
    
    precondition(error == 0, "Failed to compute Cholesky decomposition")
    
    U = toRows(U, .Column)
    
    return tri(U, t)
}

/// Return the upper/lower triangular part of a given matrix.
///
/// - Parameters:
///     - A: matrix
///     - t: Triangle value (.Upper, .Lower)
/// - Returns: upper/lower triangular part
public func tri(_ A: Matrix, _ t: Triangle) -> Matrix {
    let _A = zeros(A.rows, A.cols)
    switch t {
    case .Upper:
        for i in (0..<A.rows) {
            if(i == A.cols) {
                break
            }
            _A[i...i, i..<A.cols] = A[i...i, i..<A.cols]
        }
    case .Lower:
        for i in (0..<A.rows) {
            _A[i...i, 0...min(i, A.cols-1)] = A[i...i, 0...min(i, A.cols-1)]
        }
    }
    return _A
}

//public func dgels_(_ __trans: UnsafeMutablePointer<Int8>!, _ __m: UnsafeMutablePointer<__CLPK_integer>!, _ __n: UnsafeMutablePointer<__CLPK_integer>!, _ __nrhs: UnsafeMutablePointer<__CLPK_integer>!, _ __a: UnsafeMutablePointer<__CLPK_doublereal>!, _ __lda: UnsafeMutablePointer<__CLPK_integer>!, _ __b: UnsafeMutablePointer<__CLPK_doublereal>!, _ __ldb: UnsafeMutablePointer<__CLPK_integer>!, _ __work: UnsafeMutablePointer<__CLPK_doublereal>!, _ __lwork: UnsafeMutablePointer<__CLPK_integer>!, _ __info: UnsafeMutablePointer<__CLPK_integer>!) -> Int32

/// Compute the solution to the system Ax = b for each column in B. With A(m x n): if m>=n, then it finds the least squares solution
/// to the overdetermined probelm by minimizing ||B - A*x||; if n>m, then it finds the minimum norm solution to the underdetermined system Ax = B.
///
///    A precondition error is thrown if the given matrix is singular, algorithm fails to converge or .
///
/// - Parameters:
///     - A: the A matrix in Ax = B
///     - B: the B matrix in Ax = B
/// - Returns: the solution(s) to the system Ax = B
public func linsolve(_ A: Matrix, _ B: Matrix) throws -> Matrix {
    guard A.rows == B.rows else {
        throw LASwiftError.misc("Matrix dimensions must agree")
    }

    let A_copy = Matrix(A)
    let B_copy = Matrix(B)
    
    var TRANS = "T".utf8CString[0]
    var M = __CLPK_integer(A_copy.rows)
    var N = __CLPK_integer(A_copy.cols)
    var NRHS = __CLPK_integer(B_copy.cols)
    
    var LDA = M
    var LDB = M
    
    var wkOpt = __CLPK_doublereal(0.0)
    var lWork = __CLPK_integer(-1)
    
    var error: __CLPK_integer = 0
    
    /* Query and allocate the optimal workspace */
    
    dgels_(&TRANS, &M, &N, &NRHS, &(A_copy.flat), &LDA, &(B_copy.flat), &LDB, &wkOpt, &lWork, &error)
    
    lWork = __CLPK_integer(wkOpt)
    var work = Vector(repeating: 0.0, count: Int(lWork))
    
    /* Compute inversed matrix */
    
    dgels_(&TRANS, &M, &N, &NRHS, &(A_copy.flat), &LDA, &(B_copy.flat), &LDB, &work, &lWork, &error)
    
    guard error == 0 else {
        throw LASwiftError.misc("Matrix A is non invertible")
    }
    
    return B_copy[0..<A.cols, 0..<B.cols]
}
