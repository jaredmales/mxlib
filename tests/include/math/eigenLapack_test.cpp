/** \file eigenLapack_test.cpp
 * \brief Tests Eigen-compatible BLAS/LAPACK helpers.
 */
#include "../../catch2/catch.hpp"

#include "../../../include/math/eigenLapack.hpp"

/** \cond */
/// Compile the symmetric-eigensolver workspace implementation for a representative scalar type.
template struct mx::math::syevrMem<double>;
/** \endcond */

/** \defgroup eigenLapack_unit_tests eigenLapack Unit Tests
 * \ingroup math_unit_tests
 */

namespace unitTest::math_eigenLapack_test
{

/** \brief Verifies covariance construction and KL-mode calculation for a small reference matrix.
 *
 * \ingroup eigenLapack_unit_tests
 */
TEST_CASE( "eigenSYRK and calcKLModes form normalized KL modes", "[math::eigenLapack]" )
{
    Eigen::Array<double, Eigen::Dynamic, Eigen::Dynamic> references( 3, 2 );
    references << 1.0, 0.0, 0.0, 2.0, 0.0, 0.0;

    Eigen::Array<double, Eigen::Dynamic, Eigen::Dynamic> covariance;
    mx::math::eigenSYRK( covariance, references );

    REQUIRE( covariance.rows() == 2 );
    REQUIRE( covariance.cols() == 2 );
    REQUIRE( covariance( 0, 0 ) == Approx( 1.0 ) );
    REQUIRE( covariance( 1, 0 ) == Approx( 0.0 ) );
    REQUIRE( covariance( 1, 1 ) == Approx( 4.0 ) );

    Eigen::Array<double, Eigen::Dynamic, Eigen::Dynamic> modes;
    double eigenTime{ -1 };
    double modeTime{ -1 };
    REQUIRE( mx::math::calcKLModes<double>( modes, covariance, references, 0, nullptr, &eigenTime, &modeTime ) == 0 );
    REQUIRE( modes.rows() == 2 );
    REQUIRE( modes.cols() == 3 );
    REQUIRE( modes.matrix().row( 0 ).norm() == Approx( 1.0 ).margin( 1e-12 ) );
    REQUIRE( modes.matrix().row( 1 ).norm() == Approx( 1.0 ).margin( 1e-12 ) );
    REQUIRE( modes.matrix().row( 0 ).dot( modes.matrix().row( 1 ) ) == Approx( 0.0 ).margin( 1e-12 ) );
    REQUIRE( eigenTime >= 0 );
    REQUIRE( modeTime >= 0 );

    mx::math::syevrMem<double> workspace;
    REQUIRE( mx::math::calcKLModes( modes, covariance, references, 1, &workspace ) == 0 );
    REQUIRE( modes.rows() == 1 );
    REQUIRE( modes.cols() == 3 );
    REQUIRE( modes.matrix().row( 0 ).norm() == Approx( 1.0 ).margin( 1e-12 ) );

    const auto *work = workspace.work;
    const auto *iWork = workspace.iWork;
    REQUIRE( mx::math::calcKLModes( modes, covariance, references, 1, &workspace ) == 0 );
    REQUIRE( workspace.work == work );
    REQUIRE( workspace.iWork == iWork );
}

/** \brief Verifies that calcKLModes rejects incompatible and non-square covariance matrices without leaking workspaces.
 *
 * \ingroup eigenLapack_unit_tests
 */
TEST_CASE( "calcKLModes rejects invalid covariance geometry", "[math::eigenLapack]" )
{
    Eigen::Array<double, Eigen::Dynamic, Eigen::Dynamic> modes;

    Eigen::Array<double, Eigen::Dynamic, Eigen::Dynamic> covariance( 2, 2 );
    covariance.matrix().setIdentity();
    Eigen::Array<double, Eigen::Dynamic, Eigen::Dynamic> wrongReferenceCount( 3, 3 );
    wrongReferenceCount.setZero();
    REQUIRE( mx::math::calcKLModes( modes, covariance, wrongReferenceCount ) == -1 );

    Eigen::Array<double, Eigen::Dynamic, Eigen::Dynamic> nonSquareCovariance( 2, 3 );
    nonSquareCovariance.setZero();
    Eigen::Array<double, Eigen::Dynamic, Eigen::Dynamic> references( 3, 2 );
    references.setZero();
    REQUIRE( mx::math::calcKLModes( modes, nonSquareCovariance, references ) == -1 );
}

} // namespace unitTest::math_eigenLapack_test
