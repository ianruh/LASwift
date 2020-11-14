import XCTest
import Quick

@testable import LASwiftTests

QCKMain([
    MatrixSpec.self,
    VectorSpec.self
], testCases: [
    testCase(PerformanceTests.allTests)
])

//XCTMain([
//     testCase(PerformanceTests.allTests)
//])
