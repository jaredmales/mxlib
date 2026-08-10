/** \file error_t_test.cpp
 * \brief Tests error-code names, messages, and conversions.
 */
#include "../../catch2/catch.hpp"

#include "../../../include/error/error.hpp"

#define MXLIBTEST_ERROR_T_TESTS
#include <error/error_t.hpp>

namespace unitTest
{
namespace errorTest
{
namespace error_tTest
{

/// Test error_t code names
/**
 * \ingroup error_error_t_unit_tests
 */
TEST_CASE( "Test error_t code names", "[error::error_t]" )
{
    std::vector<mx::error_t> errcs;
    mx::error_t_vector( errcs );

    REQUIRE(errcs.size() > 0);

    size_t fail = 0;
    for(auto & errc : errcs)
    {
        std::string name = mx::errorName(errc);
        if(name == "" || name == "unknown error_t (bug)")
        {
            ++fail;
        }
    }

    REQUIRE(fail == 0);

    int nv = static_cast<int>(*std::max_element(errcs.begin(), errcs.end()));

    REQUIRE(nv != std::numeric_limits<int>::max()); //if this fails it means the following test is invalid

    REQUIRE(mx::errorName(static_cast<mx::error_t>(std::numeric_limits<int>::max())) == "unknown error_t (bug)");


}

/// Test error_t messages
/**
 * \ingroup error_error_t_unit_tests
 */
TEST_CASE( "Test error_t messages", "[error::error_t]" )
{
    std::vector<mx::error_t> errcs;
    mx::error_t_vector( errcs );

    REQUIRE(errcs.size() > 0);

    size_t fail = 0;
    for(auto & errc : errcs)
    {
        std::string msg = mx::errorMessage(errc);
        if(msg == "" || msg == "unknown error_t (bug)")
        {
            ++fail;
        }
    }

    REQUIRE(fail == 0);

    int nv = static_cast<int>(*std::max_element(errcs.begin(), errcs.end()));

    REQUIRE(nv != std::numeric_limits<int>::max()); //if this fails it means the following test is invalid

    REQUIRE(mx::errorMessage(static_cast<mx::error_t>(std::numeric_limits<int>::max())) == "unknown error_t (bug)");
}

/// Test errno conversions
/**
 * \ingroup error_error_t_unit_tests
 */
TEST_CASE( "Test errno conversions", "[error::error_t]" )
{
    std::vector<int> errnos;
    mx::errno_vector( errnos );

    REQUIRE(errnos.size() > 0);

    size_t fail = 0;
    for(auto & en : errnos)
    {
        mx::error_t errc = mx::errno2error_t(en);
        if(errc == mx::error_t::error)
        {
            ++fail;
        }
    }

    REQUIRE(fail == 0);

    int nv = *std::max_element(errnos.begin(), errnos.end());

    REQUIRE(nv != std::numeric_limits<int>::max()); //if this fails it means the following test is invalid

    REQUIRE(mx::errno2error_t(std::numeric_limits<int>::max()) == mx::error_t::error);
}

/// Test FITS error conversions
/**
 * \ingroup error_error_t_unit_tests
 */
TEST_CASE( "Test FITS error conversions", "[error::error_t]" )
{
    std::vector<int> fitserrs;
    mx::fitserr_vector( fitserrs );

    REQUIRE(fitserrs.size() > 0);

    size_t fail = 0;
    for(auto & fe : fitserrs)
    {
        mx::error_t errc = mx::fits_status2error_t(fe);
        if(errc == mx::error_t::error)
        {
            ++fail;
        }
    }

    REQUIRE(fail == 0);

    int nv = *std::max_element(fitserrs.begin(), fitserrs.end());

    REQUIRE(nv != std::numeric_limits<int>::max()); //if this fails it means the following test is invalid

    REQUIRE(mx::fits_status2error_t(std::numeric_limits<int>::max()) == mx::error_t::error);
}

} // namespace error_tTest
} // namespace errorTest
} // namespace unitTest
