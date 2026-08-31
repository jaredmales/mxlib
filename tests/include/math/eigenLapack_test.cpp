/** \file eigenLapack_test.cpp
 * \brief Tests Eigen-compatible BLAS/LAPACK helpers.
 */
#include "../../catch2/catch.hpp"

#define MX_NO_ERROR_REPORTS

#include "../../../include/math/eigenLapack.hpp"

#include <algorithm>
#include <cstdlib>
#include <limits>
#include <type_traits>
#include <vector>

/** \cond */
/// Compile the symmetric-eigensolver workspace implementation for a representative scalar type.
template struct mx::math::syevrMem<double>;

/** Verify that the owning mx::math::syevrMem workspace cannot be copied. */
TEST_CASE( "syevrMem workspaces are non-copyable", "[math::syevrMem]" )
{
    REQUIRE_FALSE( std::is_copy_constructible_v<mx::math::syevrMem<double>> );
    REQUIRE_FALSE( std::is_copy_assignable_v<mx::math::syevrMem<double>> );
}

using matrixT = Eigen::Array<double, Eigen::Dynamic, Eigen::Dynamic>;

int allocationCall{ 0 };
int failAllocationOnCall{ 0 };
int queryCall{ 0 };
int solveCall{ 0 };
MXLAPACK_INT queryInfo{ 0 };
MXLAPACK_INT solveInfo{ 0 };
bool overrideQuerySizes{ false };
double injectedWorkSize{ 0.0 };
MXLAPACK_INT injectedIntegerWorkSize{ 0 };
std::vector<double> injectedEigenvalues;
int nonfiniteEigenvectorRow{ -1 };
int nonfiniteEigenvectorColumn{ -1 };

/// Reset all deterministic eigenLapack failure controls.
void resetControls()
{
    allocationCall = 0;
    failAllocationOnCall = 0;
    queryCall = 0;
    solveCall = 0;
    queryInfo = 0;
    solveInfo = 0;
    overrideQuerySizes = false;
    injectedWorkSize = 0.0;
    injectedIntegerWorkSize = 0;
    injectedEigenvalues.clear();
    nonfiniteEigenvectorRow = -1;
    nonfiniteEigenvectorColumn = -1;
    mx::math::detail::eigenLapackTestHooks<double>::allocator = nullptr;
    mx::math::detail::eigenLapackTestHooks<double>::solver = nullptr;
}

/// Restore production allocation and solver functions at scope exit.
struct hookReset
{
    /// Start a test scope from the production hook state.
    hookReset()
    {
        resetControls();
    }

    /// Restore the production hook state even when a Catch2 assertion aborts the scope.
    ~hookReset()
    {
        resetControls();
    }
};

/// Return a diagonal square matrix with the supplied diagonal entries.
matrixT diagonalMatrix( const std::vector<double> &diagonal /**< [in] diagonal entries */ )
{
    matrixT matrix( diagonal.size(), diagonal.size() );
    matrix.setZero();
    for( std::size_t index = 0; index < diagonal.size(); ++index )
    {
        matrix( index, index ) = diagonal[index];
    }
    return matrix;
}

/// Allocate workspace storage unless the configured allocation call must fail.
void *controlledAllocation( std::size_t bytes /**< [in] requested byte count */ )
{
    ++allocationCall;
    if( allocationCall == failAllocationOnCall )
    {
        return nullptr;
    }
    return std::malloc( bytes );
}

/// Provide controlled SYEVR query, failure, and result behavior.
MXLAPACK_INT controlledSyevr( char jobz,              /**< [in] eigenvector request */
                              char range,             /**< [in] eigenvalue selection mode */
                              char uplo,              /**< [in] populated input triangle */
                              MXLAPACK_INT n,         /**< [in] matrix order */
                              double *matrix,         /**< [in,out] input matrix */
                              MXLAPACK_INT lda,       /**< [in] input leading dimension */
                              double valueLower,      /**< [in] lower value bound */
                              double valueUpper,      /**< [in] upper value bound */
                              MXLAPACK_INT indexLow,  /**< [in] one-based lower index */
                              MXLAPACK_INT indexHigh, /**< [in] one-based upper index */
                              double tolerance,       /**< [in] convergence tolerance */
                              MXLAPACK_INT *found,    /**< [out] selected eigenvalue count */
                              double *eigenvalues,    /**< [out] eigenvalues */
                              double *eigenvectors,   /**< [out] eigenvectors */
                              MXLAPACK_INT ldz,       /**< [in] eigenvector leading dimension */
                              MXLAPACK_INT *support,  /**< [out] eigenvector support */
                              double *work,           /**< [in,out] floating workspace */
                              MXLAPACK_INT lwork,     /**< [in] floating workspace size */
                              MXLAPACK_INT *iwork,    /**< [in,out] integer workspace */
                              MXLAPACK_INT liwork /**< [in] integer workspace size */ )
{
    static_cast<void>( jobz );
    static_cast<void>( uplo );
    static_cast<void>( matrix );
    static_cast<void>( lda );
    static_cast<void>( valueLower );
    static_cast<void>( valueUpper );
    static_cast<void>( tolerance );
    static_cast<void>( support );

    if( lwork == -1 || liwork == -1 )
    {
        ++queryCall;
        if( queryInfo != 0 )
        {
            return queryInfo;
        }

        work[0] = overrideQuerySizes ? injectedWorkSize : std::max<MXLAPACK_INT>( 1, 26 * n );
        iwork[0] = overrideQuerySizes ? injectedIntegerWorkSize : std::max<MXLAPACK_INT>( 1, 10 * n );
        return 0;
    }

    ++solveCall;
    if( solveInfo != 0 )
    {
        return solveInfo;
    }

    const MXLAPACK_INT first = range == 'I' ? indexLow - 1 : 0;
    const MXLAPACK_INT count = range == 'I' ? indexHigh - indexLow + 1 : n;
    *found = count;

    for( MXLAPACK_INT column = 0; column < count; ++column )
    {
        const MXLAPACK_INT sourceColumn = first + column;
        eigenvalues[column] = injectedEigenvalues.empty() ? sourceColumn + 1 : injectedEigenvalues[sourceColumn];

        for( MXLAPACK_INT row = 0; row < n; ++row )
        {
            eigenvectors[column * ldz + row] = row == sourceColumn ? 1.0 : 0.0;
        }
    }

    if( nonfiniteEigenvectorRow >= 0 && nonfiniteEigenvectorColumn >= 0 )
    {
        eigenvectors[nonfiniteEigenvectorColumn * ldz + nonfiniteEigenvectorRow] =
            std::numeric_limits<double>::infinity();
    }

    return 0;
}
/** \endcond */

/** \defgroup eigenLapack_unit_tests eigenLapack Unit Tests
 * \ingroup math_unit_tests
 */

namespace unitTest::math_eigenLapack_test
{

/** \brief Verifies eigenSYEVR full/ranged solves, triangle selection, workspace reuse, and explicit release.
 *
 * \ingroup eigenLapack_unit_tests
 */
TEST_CASE( "eigenSYEVR solves selected ranges and reuses workspaces", "[math::eigenLapack][eigenSYEVR]" )
{
    hookReset reset;

    matrixT eigenvectors;
    matrixT eigenvalues;

    matrixT covariance = diagonalMatrix( { 1.0, 4.0, 9.0 } );
    REQUIRE( mx::math::eigenSYEVR( eigenvectors, eigenvalues, covariance ) == 0 );
    REQUIRE( eigenvectors.rows() == 3 );
    REQUIRE( eigenvectors.cols() == 3 );
    REQUIRE( eigenvalues.rows() == 3 );
    REQUIRE( eigenvalues( 0 ) == Approx( 1.0 ) );
    REQUIRE( eigenvalues( 1 ) == Approx( 4.0 ) );
    REQUIRE( eigenvalues( 2 ) == Approx( 9.0 ) );

    mx::math::syevrMem<double> workspace;
    covariance = diagonalMatrix( { 1.0, 4.0, 9.0 } );
    REQUIRE( mx::math::eigenSYEVR( eigenvectors, eigenvalues, covariance, 1, 3, 'L', &workspace ) == 0 );
    REQUIRE( eigenvectors.rows() == 3 );
    REQUIRE( eigenvectors.cols() == 2 );
    REQUIRE( eigenvalues( 0 ) == Approx( 4.0 ) );
    REQUIRE( eigenvalues( 1 ) == Approx( 9.0 ) );

    const auto *support = workspace.iSuppZ;
    const auto *minimumWork = workspace.minWork;
    const auto *work = workspace.work;
    const auto *minimumIntegerWork = workspace.minIWork;
    const auto *integerWork = workspace.iWork;

    covariance = diagonalMatrix( { 1.0, 4.0, 9.0 } );
    REQUIRE( mx::math::eigenSYEVR( eigenvectors, eigenvalues, covariance, 1, 3, 'L', &workspace ) == 0 );
    REQUIRE( workspace.iSuppZ == support );
    REQUIRE( workspace.minWork == minimumWork );
    REQUIRE( workspace.work == work );
    REQUIRE( workspace.minIWork == minimumIntegerWork );
    REQUIRE( workspace.iWork == integerWork );

    covariance = diagonalMatrix( { 1.0, 2.0, 3.0, 4.0, 5.0, 6.0 } );
    REQUIRE( mx::math::eigenSYEVR( eigenvectors, eigenvalues, covariance, 0, -1, 'L', &workspace ) == 0 );
    REQUIRE( workspace.sizeISuppZ >= 12 );
    REQUIRE( workspace.sizeMinWork >= 156 );
    REQUIRE( workspace.sizeMinIWork >= 60 );

    matrixT upperCovariance( 2, 2 );
    upperCovariance( 0, 0 ) = 2.0;
    upperCovariance( 0, 1 ) = 1.0;
    upperCovariance( 1, 0 ) = 99.0;
    upperCovariance( 1, 1 ) = 2.0;
    REQUIRE( mx::math::eigenSYEVR( eigenvectors, eigenvalues, upperCovariance, 0, -1, 'U', &workspace ) == 0 );
    REQUIRE( eigenvalues( 0 ) == Approx( 1.0 ) );
    REQUIRE( eigenvalues( 1 ) == Approx( 3.0 ) );

    workspace.free();
    REQUIRE( workspace.iSuppZ == nullptr );
    REQUIRE( workspace.minWork == nullptr );
    REQUIRE( workspace.work == nullptr );
    REQUIRE( workspace.minIWork == nullptr );
    REQUIRE( workspace.iWork == nullptr );
    REQUIRE( workspace.sizeISuppZ == 0 );
    REQUIRE( workspace.sizeMinWork == 0 );
    REQUIRE( workspace.sizeWork == 0 );
    REQUIRE( workspace.sizeMinIWork == 0 );
    REQUIRE( workspace.sizeIWork == 0 );
    REQUIRE( workspace.n == 0 );

    covariance = diagonalMatrix( { 2.0, 5.0 } );
    REQUIRE( mx::math::eigenSYEVR( eigenvectors, eigenvalues, covariance, 0, -1, 'L', &workspace ) == 0 );
    REQUIRE( eigenvalues( 0 ) == Approx( 2.0 ) );
    REQUIRE( eigenvalues( 1 ) == Approx( 5.0 ) );

    workspace.free();
    workspace.free();
}

/** \brief Verifies eigenSYEVR rejects empty, non-square, and invalid half-open eigenvalue ranges before LAPACK.
 *
 * \ingroup eigenLapack_unit_tests
 */
TEST_CASE( "eigenSYEVR rejects invalid matrix and range geometry", "[math::eigenLapack][eigenSYEVR]" )
{
    hookReset reset;

    matrixT eigenvectors;
    matrixT eigenvalues;

    matrixT empty( 0, 0 );
    REQUIRE( mx::math::eigenSYEVR( eigenvectors, eigenvalues, empty ) == -1 );

    matrixT rectangular( 2, 3 );
    rectangular.setZero();
    REQUIRE( mx::math::eigenSYEVR( eigenvectors, eigenvalues, rectangular ) == -1 );

    matrixT covariance = diagonalMatrix( { 1.0, 2.0 } );
    REQUIRE( mx::math::eigenSYEVR( eigenvectors, eigenvalues, covariance, 0, 0 ) == -1 );
    REQUIRE( mx::math::eigenSYEVR( eigenvectors, eigenvalues, covariance, -1, 1 ) == -1 );
    REQUIRE( mx::math::eigenSYEVR( eigenvectors, eigenvalues, covariance, 0, 3 ) == -1 );
    REQUIRE( mx::math::eigenSYEVR( eigenvectors, eigenvalues, covariance, 0, -2 ) == -1 );
}

/** \brief Verifies eigenSYEVR propagates every workspace-allocation failure and both LAPACK solver failures.
 *
 * \ingroup eigenLapack_unit_tests
 */
TEST_CASE( "eigenSYEVR reports allocation and solver failures", "[math::eigenLapack][eigenSYEVR]" )
{
    hookReset reset;
    matrixT eigenvectors;
    matrixT eigenvalues;

    mx::math::detail::eigenLapackTestHooks<double>::allocator = &controlledAllocation;
    failAllocationOnCall = 1;
    matrixT covariance = diagonalMatrix( { 1.0, 2.0 } );
    REQUIRE( mx::math::eigenSYEVR( eigenvectors, eigenvalues, covariance ) == -1000 );

    resetControls();
    mx::math::detail::eigenLapackTestHooks<double>::allocator = &controlledAllocation;
    failAllocationOnCall = 5;
    covariance = diagonalMatrix( { 1.0, 2.0 } );
    REQUIRE( mx::math::eigenSYEVR( eigenvectors, eigenvalues, covariance ) == -1000 );

    for( int failedCall = 1; failedCall <= 5; ++failedCall )
    {
        resetControls();
        mx::math::detail::eigenLapackTestHooks<double>::allocator = &controlledAllocation;
        failAllocationOnCall = failedCall;

        mx::math::syevrMem<double> workspace;
        covariance = diagonalMatrix( { 1.0, 2.0 } );
        REQUIRE( mx::math::eigenSYEVR( eigenvectors, eigenvalues, covariance, 0, -1, 'L', &workspace ) == -1000 );
        REQUIRE( allocationCall == failedCall );
    }

    resetControls();
    mx::math::detail::eigenLapackTestHooks<double>::solver = &controlledSyevr;
    queryInfo = 17;
    covariance = diagonalMatrix( { 1.0, 2.0 } );
    REQUIRE( mx::math::eigenSYEVR( eigenvectors, eigenvalues, covariance ) == 17 );
    REQUIRE( queryCall == 1 );
    REQUIRE( solveCall == 0 );

    resetControls();
    mx::math::detail::eigenLapackTestHooks<double>::solver = &controlledSyevr;
    solveInfo = 23;
    mx::math::syevrMem<double> workspace;
    covariance = diagonalMatrix( { 1.0, 2.0 } );
    REQUIRE( mx::math::eigenSYEVR( eigenvectors, eigenvalues, covariance, 0, -1, 'L', &workspace ) == 23 );
    REQUIRE( queryCall == 1 );
    REQUIRE( solveCall == 1 );
}

/** \brief Verifies eigenSYEVR rejects nonpositive, nonfinite, and unrepresentable queried workspace sizes.
 *
 * \ingroup eigenLapack_unit_tests
 */
TEST_CASE( "eigenSYEVR validates queried workspace sizes", "[math::eigenLapack][eigenSYEVR]" )
{
    hookReset reset;
    matrixT eigenvectors;
    matrixT eigenvalues;
    matrixT covariance;

    const std::vector<double> invalidFloatingSizes{ 0.0,
                                                    -1.0,
                                                    std::numeric_limits<double>::quiet_NaN(),
                                                    std::numeric_limits<double>::infinity(),
                                                    static_cast<double>( std::numeric_limits<MXLAPACK_INT>::max() ) +
                                                        1.0 };

    for( const double invalidSize : invalidFloatingSizes )
    {
        resetControls();
        mx::math::detail::eigenLapackTestHooks<double>::solver = &controlledSyevr;
        overrideQuerySizes = true;
        injectedWorkSize = invalidSize;
        injectedIntegerWorkSize = 10;
        covariance = diagonalMatrix( { 1.0, 2.0 } );
        REQUIRE( mx::math::eigenSYEVR( eigenvectors, eigenvalues, covariance ) == -1 );
        REQUIRE( solveCall == 0 );
    }

    for( const MXLAPACK_INT invalidSize : { MXLAPACK_INT{ 0 }, MXLAPACK_INT{ -1 } } )
    {
        resetControls();
        mx::math::detail::eigenLapackTestHooks<double>::solver = &controlledSyevr;
        overrideQuerySizes = true;
        injectedWorkSize = 10;
        injectedIntegerWorkSize = invalidSize;
        mx::math::syevrMem<double> workspace;
        covariance = diagonalMatrix( { 1.0, 2.0 } );
        REQUIRE( mx::math::eigenSYEVR( eigenvectors, eigenvalues, covariance, 0, -1, 'L', &workspace ) == -1 );
        REQUIRE( solveCall == 0 );
    }
}

/** \brief Verifies calcEigenVecs mode clamping, raw eigenvectors, scaling, and caller-owned workspace reuse.
 *
 * \ingroup eigenLapack_unit_tests
 */
TEST_CASE( "calcEigenVecs clamps mode requests and controls normalization", "[math::eigenLapack][calcEigenVecs]" )
{
    hookReset reset;
    matrixT eigenvectors;
    matrixT eigenvalues;
    matrixT covariance = diagonalMatrix( { 1.0, 4.0, 9.0 } );

    double eigenTime{ -1.0 };
    REQUIRE( mx::math::calcEigenVecs<
                 double>( eigenvectors, eigenvalues, covariance, 1, false, false, nullptr, &eigenTime ) == 0 );
    REQUIRE( eigenvectors.rows() == 3 );
    REQUIRE( eigenvectors.cols() == 1 );
    REQUIRE( eigenvalues( 0 ) == Approx( 9.0 ) );
    REQUIRE( eigenvectors.matrix().col( 0 ).norm() == Approx( 1.0 ).margin( 1e-12 ) );
    REQUIRE( eigenTime >= 0.0 );

    mx::math::syevrMem<double> workspace;
    REQUIRE( mx::math::calcEigenVecs( eigenvectors, eigenvalues, covariance, 0, false, false, &workspace ) == 0 );
    REQUIRE( eigenvectors.cols() == 3 );
    REQUIRE( eigenvalues( 0 ) == Approx( 1.0 ) );
    REQUIRE( eigenvalues( 1 ) == Approx( 4.0 ) );
    REQUIRE( eigenvalues( 2 ) == Approx( 9.0 ) );

    REQUIRE( mx::math::calcEigenVecs( eigenvectors, eigenvalues, covariance, -2, false, true, &workspace ) == 0 );
    REQUIRE( eigenvectors.cols() == 3 );

    REQUIRE( mx::math::calcEigenVecs( eigenvectors, eigenvalues, covariance, 8, false, false, &workspace ) == 0 );
    REQUIRE( eigenvectors.cols() == 3 );

    REQUIRE( mx::math::calcEigenVecs( eigenvectors, eigenvalues, covariance, 2, true, false, &workspace ) == 0 );
    REQUIRE( eigenvectors.cols() == 2 );
    REQUIRE( eigenvalues( 0 ) == Approx( 4.0 ) );
    REQUIRE( eigenvalues( 1 ) == Approx( 9.0 ) );
    REQUIRE( eigenvectors.matrix().col( 0 ).norm() == Approx( 0.5 ).margin( 1e-12 ) );
    REQUIRE( eigenvectors.matrix().col( 1 ).norm() == Approx( 1.0 / 3.0 ).margin( 1e-12 ) );
}

/** \brief Verifies calcEigenVecs check-mode handling of zero, negative, tiny, and nonfinite solver results.
 *
 * \ingroup eigenLapack_unit_tests
 */
TEST_CASE( "calcEigenVecs validates exceptional normalized results", "[math::eigenLapack][calcEigenVecs]" )
{
    hookReset reset;
    mx::math::detail::eigenLapackTestHooks<double>::solver = &controlledSyevr;
    injectedEigenvalues = { 0.0, -1.0, std::numeric_limits<double>::quiet_NaN(), 1e-20, 4.0 };
    nonfiniteEigenvectorRow = 0;
    nonfiniteEigenvectorColumn = 4;

    matrixT covariance = diagonalMatrix( { 1.0, 2.0, 3.0, 4.0, 5.0 } );
    matrixT eigenvectors;
    matrixT eigenvalues;
    REQUIRE( mx::math::calcEigenVecs( eigenvectors, eigenvalues, covariance, 0, true, true ) == 0 );

    REQUIRE( eigenvectors.rows() == 5 );
    REQUIRE( eigenvectors.cols() == 5 );
    REQUIRE( eigenvectors.matrix().col( 0 ).norm() == Approx( 0.0 ) );
    REQUIRE( eigenvectors.matrix().col( 1 ).norm() == Approx( 0.0 ) );
    REQUIRE( eigenvectors.matrix().col( 2 ).norm() == Approx( 0.0 ) );
    REQUIRE( eigenvectors.matrix().col( 3 ).norm() == Approx( 1e10 ).epsilon( 1e-12 ) );
    REQUIRE( eigenvectors.matrix().col( 4 ).norm() == Approx( 0.0 ) );
    REQUIRE( eigenvectors.isFinite().all() );
    REQUIRE( eigenvalues( 0 ) == Approx( 0.0 ) );
    REQUIRE( eigenvalues( 1 ) == Approx( -1.0 ) );
    REQUIRE( mx::math::isNan( eigenvalues( 2 ) ) );
    REQUIRE( eigenvalues( 3 ) == Approx( 1e-20 ) );
    REQUIRE( eigenvalues( 4 ) == Approx( 4.0 ) );
}

/** \brief Verifies calcEigenVecs rejects invalid covariance geometry and maps query/solve failures to -1.
 *
 * \ingroup eigenLapack_unit_tests
 */
TEST_CASE( "calcEigenVecs rejects invalid inputs and solver failures", "[math::eigenLapack][calcEigenVecs]" )
{
    hookReset reset;
    matrixT eigenvectors;
    matrixT eigenvalues;

    matrixT rectangular( 2, 3 );
    rectangular.setZero();
    REQUIRE( mx::math::calcEigenVecs( eigenvectors, eigenvalues, rectangular ) == -1 );

    matrixT empty( 0, 0 );
    mx::math::syevrMem<double> workspace;
    REQUIRE( mx::math::calcEigenVecs( eigenvectors, eigenvalues, empty, 0, false, false, &workspace ) == -1 );

    mx::math::detail::eigenLapackTestHooks<double>::solver = &controlledSyevr;
    queryInfo = 29;
    matrixT covariance = diagonalMatrix( { 1.0, 2.0 } );
    REQUIRE( mx::math::calcEigenVecs( eigenvectors, eigenvalues, covariance ) == -1 );

    resetControls();
    mx::math::detail::eigenLapackTestHooks<double>::solver = &controlledSyevr;
    solveInfo = 31;
    REQUIRE( mx::math::calcEigenVecs( eigenvectors, eigenvalues, covariance, 0, false, false, &workspace ) == -1 );
}

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
