import XCTest
import Quick

@testable import LASwiftTests

//QCKMain([
//    MatrixSpec.self,
//    VectorSpec.self
//])

XCTMain([
     testCase(PerformanceTests.allTests)
])
