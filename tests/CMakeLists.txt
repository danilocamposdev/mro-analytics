include(FetchContent)
FetchContent_Declare(
    catch2
    GIT_REPOSITORY https://github.com/catchorg/Catch2.git
    GIT_TAG        v3.6.0
)
FetchContent_MakeAvailable(catch2)

add_executable(mro_analytics_tests
    test_croston.cpp
    test_croston_professor_example.cpp
    test_moments.cpp
    test_classic_ad_hoc.cpp
)

target_link_libraries(mro_analytics_tests
    PRIVATE
        ${PROJECT_NAME}::${PROJECT_NAME}
        Catch2::Catch2WithMain
)

include(CTest)
include(Catch)
catch_discover_tests(mro_analytics_tests)
