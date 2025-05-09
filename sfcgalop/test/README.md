# SFCGALOP Unit Tests

Unit tests for SFCGALOP using Boost.Test framework.

## Building Tests

```bash
cd build
cmake -DBUILD_TESTING=ON ..
make test_sfcgalop
```

## Running Tests

```bash
# Run all tests
./sfcgalop/test/test_sfcgalop

# List available tests
./sfcgalop/test/test_sfcgalop --list_content

# Run specific test
./sfcgalop/test/test_sfcgalop --run_test=test_area

# Verbose output
./sfcgalop/test/test_sfcgalop --log_level=all

# Generate report
./sfcgalop/test/test_sfcgalop --report_level=detailed
```

## Test Coverage

The test suite includes one test per major functionality:

1. **test_load_wkt** - Load geometry from WKT string
2. **test_validate** - Validate polygon geometry
3. **test_area** - Calculate polygon area
4. **test_distance** - Calculate distance between points
5. **test_operations_list** - Get all operations list
6. **test_exception** - Exception handling
7. **test_output_wkt** - Output formatting
8. **test_intersection** - Intersection operation
9. **test_convexhull** - Convex hull operation
10. **test_is_valid** - Validity predicate

## Test Results

All 10 tests pass successfully:

```text
Running 10 test cases...
*** No errors detected
```

## Test Architecture

- **test_minimal.cpp** - Main test file with 10 focused tests
- Uses Boost.Test for test framework
- Links against sfcgalop libraries
- Tests core functionality without external dependencies

## Adding New Tests

To add a new test:

```cpp
BOOST_AUTO_TEST_CASE(test_new_feature)
{
    // Setup
    SFCGAL::Point point(1, 2);

    // Execute
    auto result = Operations::execute_operation("new_op", "", &point, nullptr);

    // Verify
    BOOST_CHECK(result.has_value());
}
```

## Test Philosophy

Each test follows the principle of "one test per function":
- Single responsibility per test
- Clear test names describing what is tested
- Minimal setup and teardown
- Fast execution (< 1ms per test)
