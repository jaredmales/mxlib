/** \file error_test.cpp
 */
#include "../../catch2/catch.hpp"

#include "../../../include/error/error.hpp"

namespace unitTest
{
namespace errorTest
{
namespace errorTest
{
/// Test error_t boolean operators and functions
/**
 * \ingroup error_error_unit_tests
 */
TEST_CASE( "Test error_t boolean operators and functions", "[error::error]" )
{
    #if MXLIBTEST_DOXYGEN_REF //doxygen links
    bool tf = mx::operator!(mx::error_t::noerror);
    tf = mx::operator==(mx::error_t::noerror, true);
    tf = mx::operator!=(mx::error_t::noerror, true);
    #endif

    SECTION( "operator!" )
    {
        mx::error_t errc = mx::error_t::noerror;

        REQUIRE(!errc == true);
        REQUIRE(!!errc == false);

        errc = mx::error_t::error;
        REQUIRE(!errc == false);
        REQUIRE(!!errc == true);
    }

    SECTION( "operator==" )
    {
        mx::error_t errc = mx::error_t::noerror;

        REQUIRE((errc == false));
        REQUIRE(!(errc == true));

        errc = mx::error_t::error;
        REQUIRE((errc == true));
        REQUIRE(!(errc == false));
    }

    SECTION( "operator!=" )
    {
        mx::error_t errc = mx::error_t::noerror;

        REQUIRE(!(errc != false));
        REQUIRE((errc != true));

        errc = mx::error_t::error;
        REQUIRE(!(errc != true));
        REQUIRE((errc != false));
    }

    SECTION( "isError" )
    {
        mx::error_t errc = mx::error_t::noerror;

        REQUIRE(!mx::isError(errc));

        errc = mx::error_t::error;
        REQUIRE(mx::isError(errc));
    }
}

/// Internal error messages
/**
 * \ingroup error_error_unit_tests
 */
TEST_CASE( "Internal error messages", "[error::error]" )
{
    #ifdef MXLIBTEST_DOXYGEN_REF //doxygen link
    mx::internal::mxlib_error_message(mx::error_t::noerror, "");
    #endif

    SECTION("verbose::o")
    {
        std::string msg = mx::internal::mxlib_error_message<mx::verbose::o>(mx::error_t::noerror, "test");
        REQUIRE(msg == "");

        msg = mx::internal::mxlib_error_message<mx::verbose::o>(mx::error_t::error, "test");
        REQUIRE(msg == "");

        msg = mx::internal::mxlib_error_message<mx::verbose::o>(mx::error_t::noerror);
        REQUIRE(msg == "");

        msg = mx::internal::mxlib_error_message<mx::verbose::o>(mx::error_t::error);
        REQUIRE(msg == "");
    }

    SECTION("verbose::v")
    {
        std::string msg = mx::internal::mxlib_error_message<mx::verbose::v>(mx::error_t::noerror, "test");
        REQUIRE(msg != "");

        msg = mx::internal::mxlib_error_message<mx::verbose::v>(mx::error_t::error, "test");
        REQUIRE(msg != "");

        msg = mx::internal::mxlib_error_message<mx::verbose::v>(mx::error_t::noerror);
        REQUIRE(msg != "");

        msg = mx::internal::mxlib_error_message<mx::verbose::v>(mx::error_t::error);
        REQUIRE(msg != "");
    }

    SECTION("verbose::vv")
    {
        std::string msg = mx::internal::mxlib_error_message<mx::verbose::vv>(mx::error_t::noerror, "test");
        REQUIRE(msg != "");

        msg = mx::internal::mxlib_error_message<mx::verbose::vv>(mx::error_t::error, "test");
        REQUIRE(msg != "");

        msg = mx::internal::mxlib_error_message<mx::verbose::vv>(mx::error_t::noerror);
        REQUIRE(msg != "");

        msg = mx::internal::mxlib_error_message<mx::verbose::vv>(mx::error_t::error);
        REQUIRE(msg != "");
    }

    SECTION("verbose::vvv")
    {
        std::string msg = mx::internal::mxlib_error_message<mx::verbose::vvv>(mx::error_t::noerror, "test");
        REQUIRE(msg != "");

        msg = mx::internal::mxlib_error_message<mx::verbose::vvv>(mx::error_t::error, "test");
        REQUIRE(msg != "");

        msg = mx::internal::mxlib_error_message<mx::verbose::vvv>(mx::error_t::noerror);
        REQUIRE(msg != "");

        msg = mx::internal::mxlib_error_message<mx::verbose::vvv>(mx::error_t::error);
        REQUIRE(msg != "");
    }
}

/// Internal error reports
/**
 * \ingroup error_error_unit_tests
 */
TEST_CASE( "Internal error reports", "[error::error]" )
{
    #ifdef MXLIBTEST_DOXYGEN_REF //doxygen link
    mx::error_t errc = mx::internal::mxlib_error_report(mx::error_t::noerror, "");
    #endif

    SECTION("verbose::o")
    {
        mx::error_t errc = mx::internal::mxlib_error_report<mx::verbose::o>(mx::error_t::noerror, "test");
        REQUIRE(errc == mx::error_t::noerror);

        errc = mx::internal::mxlib_error_report<mx::verbose::o>(mx::error_t::error, "test");
        REQUIRE(errc == mx::error_t::error);

        errc = mx::internal::mxlib_error_report<mx::verbose::o>(mx::error_t::noerror);
        REQUIRE(errc == mx::error_t::noerror);

        errc = mx::internal::mxlib_error_report<mx::verbose::o>(mx::error_t::error);
        REQUIRE(errc == mx::error_t::error);
    }

    SECTION("verbose::v")
    {
        mx::error_t errc = mx::internal::mxlib_error_report<mx::verbose::v>(mx::error_t::noerror, "test");
        REQUIRE(errc == mx::error_t::noerror);

        errc = mx::internal::mxlib_error_report<mx::verbose::v>(mx::error_t::error, "test");
        REQUIRE(errc == mx::error_t::error);

        errc = mx::internal::mxlib_error_report<mx::verbose::v>(mx::error_t::noerror);
        REQUIRE(errc == mx::error_t::noerror);

        errc = mx::internal::mxlib_error_report<mx::verbose::v>(mx::error_t::error);
        REQUIRE(errc == mx::error_t::error);
    }

    SECTION("verbose::vv")
    {
        mx::error_t errc = mx::internal::mxlib_error_report<mx::verbose::vv>(mx::error_t::noerror, "test");
        REQUIRE(errc == mx::error_t::noerror);

        errc = mx::internal::mxlib_error_report<mx::verbose::vv>(mx::error_t::error, "test");
        REQUIRE(errc == mx::error_t::error);

        errc = mx::internal::mxlib_error_report<mx::verbose::vv>(mx::error_t::noerror);
        REQUIRE(errc == mx::error_t::noerror);

        errc = mx::internal::mxlib_error_report<mx::verbose::vv>(mx::error_t::error);
        REQUIRE(errc == mx::error_t::error);
    }

    SECTION("verbose::vvv")
    {
        mx::error_t errc = mx::internal::mxlib_error_report<mx::verbose::vvv>(mx::error_t::noerror, "test");
        REQUIRE(errc == mx::error_t::noerror);

        errc = mx::internal::mxlib_error_report<mx::verbose::vvv>(mx::error_t::error, "test");
        REQUIRE(errc == mx::error_t::error);

        errc = mx::internal::mxlib_error_report<mx::verbose::vvv>(mx::error_t::noerror);
        REQUIRE(errc == mx::error_t::noerror);

        errc = mx::internal::mxlib_error_report<mx::verbose::vvv>(mx::error_t::error);
        REQUIRE(errc == mx::error_t::error);
    }
}

/// Error messages
/**
 * \ingroup error_error_unit_tests
 */
TEST_CASE( "Error messages", "[error::error]" )
{
    #ifdef MXLIBTEST_DOXYGEN_REF //doxygen link
    mx::error_message(mx::error_t::noerror, "");
    #endif

    SECTION("verbose::o")
    {
        std::string msg = mx::error_message<mx::verbose::o>(mx::error_t::noerror, "test");
        REQUIRE(msg == "");

        msg = mx::error_message<mx::verbose::o>(mx::error_t::error, "test");
        REQUIRE(msg == "");

        msg = mx::error_message<mx::verbose::o>(mx::error_t::noerror);
        REQUIRE(msg == "");

        msg = mx::error_message<mx::verbose::o>(mx::error_t::error);
        REQUIRE(msg == "");
    }

    SECTION("verbose::v")
    {
        std::string msg = mx::error_message<mx::verbose::v>(mx::error_t::noerror, "test");
        REQUIRE(msg != "");

        msg = mx::error_message<mx::verbose::v>(mx::error_t::error, "test");
        REQUIRE(msg != "");

        msg = mx::error_message<mx::verbose::v>(mx::error_t::noerror);
        REQUIRE(msg != "");

        msg = mx::error_message<mx::verbose::v>(mx::error_t::error);
        REQUIRE(msg != "");
    }

    SECTION("verbose::vv")
    {
        std::string msg = mx::error_message<mx::verbose::vv>(mx::error_t::noerror, "test");
        REQUIRE(msg != "");

        msg = mx::error_message<mx::verbose::vv>(mx::error_t::error, "test");
        REQUIRE(msg != "");

        msg = mx::error_message<mx::verbose::vv>(mx::error_t::noerror);
        REQUIRE(msg != "");

        msg = mx::error_message<mx::verbose::vv>(mx::error_t::error);
        REQUIRE(msg != "");
    }

    SECTION("verbose::vvv")
    {
        std::string msg = mx::error_message<mx::verbose::vvv>(mx::error_t::noerror, "test");
        REQUIRE(msg != "");

        msg = mx::error_message<mx::verbose::vvv>(mx::error_t::error, "test");
        REQUIRE(msg != "");

        msg = mx::error_message<mx::verbose::vvv>(mx::error_t::noerror);
        REQUIRE(msg != "");

        msg = mx::error_message<mx::verbose::vvv>(mx::error_t::error);
        REQUIRE(msg != "");
    }
}

/// Error reports
/**
 * \ingroup error_error_unit_tests
 */
TEST_CASE( "Error reports", "[error::error]" )
{
    #ifdef MXLIBTEST_DOXYGEN_REF //doxygen link
    mx::error_t errc = mx::error_report(mx::error_t::noerror, "");
    #endif

    SECTION("verbose::o")
    {
        mx::error_t errc = mx::error_report<mx::verbose::o>(mx::error_t::noerror, "test");
        REQUIRE(errc == mx::error_t::noerror);

        errc = mx::error_report<mx::verbose::o>(mx::error_t::error, "test");
        REQUIRE(errc == mx::error_t::error);

        errc = mx::error_report<mx::verbose::o>(mx::error_t::noerror);
        REQUIRE(errc == mx::error_t::noerror);

        errc = mx::error_report<mx::verbose::o>(mx::error_t::error);
        REQUIRE(errc == mx::error_t::error);
    }

    SECTION("verbose::v")
    {
        mx::error_t errc = mx::error_report<mx::verbose::v>(mx::error_t::noerror, "test");
        REQUIRE(errc == mx::error_t::noerror);

        errc = mx::error_report<mx::verbose::v>(mx::error_t::error, "test");
        REQUIRE(errc == mx::error_t::error);

        errc = mx::error_report<mx::verbose::v>(mx::error_t::noerror);
        REQUIRE(errc == mx::error_t::noerror);

        errc = mx::error_report<mx::verbose::v>(mx::error_t::error);
        REQUIRE(errc == mx::error_t::error);
    }

    SECTION("verbose::vv")
    {
        mx::error_t errc = mx::error_report<mx::verbose::vv>(mx::error_t::noerror, "test");
        REQUIRE(errc == mx::error_t::noerror);

        errc = mx::error_report<mx::verbose::vv>(mx::error_t::error, "test");
        REQUIRE(errc == mx::error_t::error);

        errc = mx::error_report<mx::verbose::vv>(mx::error_t::noerror);
        REQUIRE(errc == mx::error_t::noerror);

        errc = mx::error_report<mx::verbose::vv>(mx::error_t::error);
        REQUIRE(errc == mx::error_t::error);
    }

    SECTION("verbose::vvv")
    {
        mx::error_t errc = mx::error_report<mx::verbose::vvv>(mx::error_t::noerror, "test");
        REQUIRE(errc == mx::error_t::noerror);

        errc = mx::error_report<mx::verbose::vvv>(mx::error_t::error, "test");
        REQUIRE(errc == mx::error_t::error);

        errc = mx::error_report<mx::verbose::vvv>(mx::error_t::noerror);
        REQUIRE(errc == mx::error_t::noerror);

        errc = mx::error_report<mx::verbose::vvv>(mx::error_t::error);
        REQUIRE(errc == mx::error_t::error);
    }
}

/// \cond
mx::error_t mx_error_check_test_fxn(mx::error_t errc)
{
    return errc;
}

mx::error_t mx_error_check_test(mx::error_t errc)
{
    typedef mx::verbose::vv verboseT;

    mx_error_check(mx_error_check_test_fxn(errc));

    return mx::error_t::exception; //a sentinel
}

mx::error_t mx_error_check_code_test(mx::error_t errc)
{
    typedef mx::verbose::vv verboseT;

    mx_error_check_code(errc);

    return mx::error_t::exception; //a sentinel
}

int mx_error_check_rv_test(mx::error_t errc)
{
    typedef mx::verbose::vv verboseT;

    mx_error_check_rv(mx_error_check_test_fxn(errc), -1);

    return 0;
}

int mx_error_check_code_rv_test(mx::error_t errc)
{
    typedef mx::verbose::vv verboseT;

    mx_error_check_code_rv(mx_error_check_test_fxn(errc), -1);

    return 0;
}

mx::error_t mx_error_return_test(mx::error_t errc)
{
    typedef mx::verbose::vv verboseT;

    mx_error_return(mx_error_check_test_fxn(errc));

}

mx::error_t mx_error_return_code_test(mx::error_t errc)
{
    typedef mx::verbose::vv verboseT;

    mx_error_return_code(errc);

}

/// \endcond

/// Error macros
/**
 * \ingroup error_error_unit_tests
 */
TEST_CASE( "Error macros", "[error::error]" )
{
    #if MXLIBTEST_DOXYGEN_REF //doxygen links
    mx_error_check(errc);
    mx_error_check_code(errc);
    mx_error_check_rv(errc);
    mx_error_check_code_rv(errc);
    mx_error_return(errc);
    mx_error_return_code(errc);
    #endif

    SECTION("mx_error_check")
    {
        //use a sentinel to test that noerror goes through the macro
        REQUIRE(mx_error_check_test(mx::error_t::noerror) == mx::error_t::exception);
        REQUIRE(mx_error_check_test(mx::error_t::error) == mx::error_t::error);
        REQUIRE(mx_error_check_test(mx::error_t::eacces) == mx::error_t::eacces);
    }

    SECTION("mx_error_check_code")
    {
        //use a sentinel to test that noerror goes through the macro
        REQUIRE(mx_error_check_code_test(mx::error_t::noerror) == mx::error_t::exception);
        REQUIRE(mx_error_check_code_test(mx::error_t::error) == mx::error_t::error);
        REQUIRE(mx_error_check_code_test(mx::error_t::eacces) == mx::error_t::eacces);
    }

    SECTION("mx_error_check_rv")
    {
        REQUIRE(mx_error_check_rv_test(mx::error_t::noerror) == 0);
        REQUIRE(mx_error_check_rv_test(mx::error_t::error) == -1);
        REQUIRE(mx_error_check_rv_test(mx::error_t::eacces) == -1);
    }

    SECTION("mx_error_check_code_rv")
    {
        REQUIRE(mx_error_check_code_rv_test(mx::error_t::noerror) == 0);
        REQUIRE(mx_error_check_code_rv_test(mx::error_t::error) == -1);
        REQUIRE(mx_error_check_code_rv_test(mx::error_t::eacces) == -1);
    }

    SECTION("mx_error_return")
    {
        REQUIRE(mx_error_return_test(mx::error_t::noerror) == mx::error_t::noerror);
        REQUIRE(mx_error_return_test(mx::error_t::error) == mx::error_t::error);
        REQUIRE(mx_error_return_test(mx::error_t::eacces) == mx::error_t::eacces);
    }

    SECTION("mx_error_return_code")
    {
        REQUIRE(mx_error_return_code_test(mx::error_t::noerror) == mx::error_t::noerror);
        REQUIRE(mx_error_return_code_test(mx::error_t::error) == mx::error_t::error);
        REQUIRE(mx_error_return_code_test(mx::error_t::eacces) == mx::error_t::eacces);
    }
}
} // namespace errorTest
} // namespace errorTest
} // namespace unitTest
