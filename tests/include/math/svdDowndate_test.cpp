/** \file
 * \brief Tests reusable thin-SVD row and column deletion updates.
 *
 * \ingroup gen_math_test_files
 */

#include "../../catch2/catch.hpp"

#include "../../../include/math/svdDowndate.hpp"

#include <Eigen/SVD>

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <limits>
#include <new>
#include <stdexcept>
#include <string>
#include <tuple>
#include <type_traits>
#include <utility>
#include <vector>

/** \cond */
namespace
{

using realT = double;
using matrixT = mx::math::svdDeletionMatrix<realT>;
using vectorT = mx::math::svdDeletionVector<realT>;
using backendT = mx::math::svdDeletionBackend;
using statusT = mx::math::svdDeletionStatus;

// Deterministic orthonormal columns from the discrete sine transform.
matrixT sineFactor( Eigen::Index rows, /* [in] number of rows */
                    Eigen::Index cols /* [in] number of columns */ )
{
    matrixT factor( rows, cols );
    const realT scale = std::sqrt( realT( 2 ) / static_cast<realT>( rows + 1 ) );
    const realT pi = std::acos( realT( -1 ) );
    for( Eigen::Index row = 0; row < rows; ++row )
    {
        for( Eigen::Index column = 0; column < cols; ++column )
        {
            factor( row, column ) = scale * std::sin( pi * static_cast<realT>( ( row + 1 ) * ( column + 1 ) ) /
                                                      static_cast<realT>( rows + 1 ) );
        }
    }
    return factor;
}

// Return a dynamic identity array.
matrixT identityMatrix( Eigen::Index size /* [in] square matrix size */ )
{
    matrixT identity( size, size );
    identity.matrix().setIdentity();
    return identity;
}

// Construct a represented matrix from thin singular factors.
matrixT representedMatrix( const matrixT &left,   /* [in] thin left factor */
                           const vectorT &values, /* [in] singular values */
                           const matrixT &right /* [in] thin right factor */ )
{
    return ( left.matrix() * values.matrix().asDiagonal() * right.matrix().transpose() ).array();
}

// Physically remove selected rows from a matrix.
matrixT retainedRows( const matrixT &matrix, /* [in] input matrix */
                      std::span<const Eigen::Index> deleted /* [in] sorted deleted rows */ )
{
    matrixT retained( matrix.rows() - static_cast<Eigen::Index>( deleted.size() ), matrix.cols() );
    Eigen::Index output{ 0 };
    std::size_t nextDeleted{ 0 };
    for( Eigen::Index row = 0; row < matrix.rows(); ++row )
    {
        if( nextDeleted < deleted.size() && row == deleted[nextDeleted] )
        {
            ++nextDeleted;
            continue;
        }
        retained.matrix().row( output++ ) = matrix.matrix().row( row );
    }
    return retained;
}

// Physically remove selected columns from a matrix.
matrixT retainedColumns( const matrixT &matrix, /* [in] input matrix */
                         std::span<const Eigen::Index> deleted /* [in] sorted deleted columns */ )
{
    matrixT retained( matrix.rows(), matrix.cols() - static_cast<Eigen::Index>( deleted.size() ) );
    Eigen::Index output{ 0 };
    std::size_t nextDeleted{ 0 };
    for( Eigen::Index column = 0; column < matrix.cols(); ++column )
    {
        if( nextDeleted < deleted.size() && column == deleted[nextDeleted] )
        {
            ++nextDeleted;
            continue;
        }
        retained.matrix().col( output++ ) = matrix.matrix().col( column );
    }
    return retained;
}

// Return direct full-SVD singular values padded to the represented rank.
vectorT directSingularValues( const matrixT &matrix, /* [in] physically retained matrix */
                              Eigen::Index rank /* [in] represented output rank */ )
{
    Eigen::JacobiSVD<Eigen::Matrix<realT, Eigen::Dynamic, Eigen::Dynamic>> svd( matrix.matrix(),
                                                                                Eigen::ComputeThinU |
                                                                                    Eigen::ComputeThinV );
    vectorT values( rank );
    values.setZero();
    const Eigen::Index copyCount = std::min<Eigen::Index>( rank, svd.singularValues().size() );
    values.matrix().head( copyCount ) = svd.singularValues().head( copyCount );
    return values;
}

// Return the preserved covariance represented by a deletion result.
matrixT representedCovariance( const matrixT &preservedFactor, /* [in] unchanged thin factor */
                               const mx::math::svdDeletionResult<realT> &result /* [in] deletion result */ )
{
    const Eigen::Matrix<realT, Eigen::Dynamic, Eigen::Dynamic> directions =
        preservedFactor.matrix() * result.rotation().matrix();
    return ( directions * result.squaredSingularValues().head( result.outputRank() ).matrix().asDiagonal() *
             directions.transpose() )
        .array();
}

// Named diagnostics returned to a TEST_CASE for local Catch2 assertions.
struct deletionComparison
{
    realT squaredSingularError{ 0 };
    realT squaredSingularTolerance{ 0 };
    realT singularError{ 0 };
    realT singularTolerance{ 0 };
    realT covarianceError{ 0 };
    realT covarianceTolerance{ 0 };
};

// Compare all-mode row-deletion outputs with a direct retained-matrix SVD.
deletionComparison compareRowResult( const matrixT &matrix,                            /* [in] represented matrix */
                                     const matrixT &right,                             /* [in] unchanged right factor */
                                     std::span<const Eigen::Index> deleted,            /* [in] deleted rows */
                                     const mx::math::svdDeletionResult<realT> &result, /* [in] deletion result */
                                     realT tolerance = 5e-11 /* [in] relative tolerance */ )
{
    const matrixT retained = retainedRows( matrix, deleted );
    const vectorT direct = directSingularValues( retained, result.baseRank() );
    const vectorT directSquared = direct.square();
    const matrixT expected = ( retained.matrix().transpose() * retained.matrix() ).array();
    const matrixT actual = representedCovariance( right, result );

    deletionComparison comparison;
    comparison.squaredSingularError = ( result.squaredSingularValues() - directSquared ).matrix().norm();
    comparison.squaredSingularTolerance = tolerance * std::max<realT>( 1, directSquared.matrix().norm() );
    comparison.singularError = ( result.singularValues() - direct ).matrix().norm();
    comparison.singularTolerance = realT( 16 ) * std::sqrt( std::numeric_limits<realT>::epsilon() ) *
                                   std::max<realT>( 1, direct.size() > 0 ? direct( 0 ) : realT( 0 ) );
    comparison.covarianceError = ( actual - expected ).matrix().norm();
    comparison.covarianceTolerance = tolerance * std::max<realT>( 1, expected.matrix().norm() );
    return comparison;
}

// Compare all-mode column-deletion outputs with a direct retained-matrix SVD.
deletionComparison compareColumnResult( const matrixT &matrix,                 /* [in] represented matrix */
                                        const matrixT &left,                   /* [in] unchanged left factor */
                                        std::span<const Eigen::Index> deleted, /* [in] deleted columns */
                                        const mx::math::svdDeletionResult<realT> &result, /* [in] deletion result */
                                        realT tolerance = 5e-11 /* [in] relative tolerance */ )
{
    const matrixT retained = retainedColumns( matrix, deleted );
    const vectorT direct = directSingularValues( retained, result.baseRank() );
    const vectorT directSquared = direct.square();
    const matrixT expected = ( retained.matrix() * retained.matrix().transpose() ).array();
    const matrixT actual = representedCovariance( left, result );

    deletionComparison comparison;
    comparison.squaredSingularError = ( result.squaredSingularValues() - directSquared ).matrix().norm();
    comparison.squaredSingularTolerance = tolerance * std::max<realT>( 1, directSquared.matrix().norm() );
    comparison.singularError = ( result.singularValues() - direct ).matrix().norm();
    comparison.singularTolerance = realT( 16 ) * std::sqrt( std::numeric_limits<realT>::epsilon() ) *
                                   std::max<realT>( 1, direct.size() > 0 ? direct( 0 ) : realT( 0 ) );
    comparison.covarianceError = ( actual - expected ).matrix().norm();
    comparison.covarianceTolerance = tolerance * std::max<realT>( 1, expected.matrix().norm() );
    return comparison;
}

mx::math::detail::svdDeletionTestOperation failingOperation =
    mx::math::detail::svdDeletionTestOperation::prepareWorkspace;

// Deterministic behavior selected for an injected solver call.
enum class solverHookMode
{
    production,
    queryFailure,
    invalidQuery,
    invalidIntegerQuery,
    solveFailure,
    countMismatch,
    nonFiniteValue,
    nonFiniteVector,
    invalidOrdering,
    outsideInterlacing,
    invalidVectorNorm,
    invalidResidual,
    negativeSpectrum,
    tinySpectrum,
    roundoffClamp,
    indefinite
};

solverHookMode syevrMode = solverHookMode::production;
solverHookMode gesvdMode = solverHookMode::production;
solverHookMode laed9Mode = solverHookMode::production;

// Throw a deterministic allocation failure for the selected operation.
void throwAllocation( mx::math::detail::svdDeletionTestOperation operation /* [in] operation being attempted */ )
{
    if( operation == failingOperation )
    {
        throw std::bad_alloc();
    }
}

// Throw a deterministic storage-length failure for the selected operation.
void throwLengthError( mx::math::detail::svdDeletionTestOperation operation /* [in] operation being attempted */ )
{
    if( operation == failingOperation )
    {
        throw std::length_error( "injected SVD deletion storage length" );
    }
}

// Emulate selected SYEVR query, failure, and malformed-output paths.
MXLAPACK_INT syevrHook( char jobz,                 /* [in] eigenvector request */
                        char range,                /* [in] eigenvalue selection mode */
                        char uplo,                 /* [in] populated triangle */
                        MXLAPACK_INT n,            /* [in] matrix order */
                        realT *matrix,             /* [in,out] input matrix */
                        MXLAPACK_INT lda,          /* [in] input leading dimension */
                        realT valueLower,          /* [in] lower value bound */
                        realT valueUpper,          /* [in] upper value bound */
                        MXLAPACK_INT indexLow,     /* [in] first requested index */
                        MXLAPACK_INT indexHigh,    /* [in] last requested index */
                        realT tolerance,           /* [in] convergence tolerance */
                        MXLAPACK_INT *found,       /* [out] returned eigenvalue count */
                        realT *eigenvalues,        /* [out] eigenvalues */
                        realT *eigenvectors,       /* [out] eigenvectors */
                        MXLAPACK_INT ldz,          /* [in] eigenvector leading dimension */
                        MXLAPACK_INT *support,     /* [out] eigenvector support */
                        realT *work,               /* [in,out] floating workspace */
                        MXLAPACK_INT lwork,        /* [in] floating workspace size */
                        MXLAPACK_INT *integerWork, /* [in,out] integer workspace */
                        MXLAPACK_INT integerWorkSize /* [in] integer workspace size */ )
{
    static_cast<void>( jobz );
    static_cast<void>( range );
    static_cast<void>( uplo );
    static_cast<void>( matrix );
    static_cast<void>( lda );
    static_cast<void>( valueLower );
    static_cast<void>( valueUpper );
    static_cast<void>( indexLow );
    static_cast<void>( indexHigh );
    static_cast<void>( tolerance );
    static_cast<void>( support );

    if( lwork == -1 || integerWorkSize == -1 )
    {
        if( syevrMode == solverHookMode::queryFailure )
        {
            return 61;
        }
        work[0] = syevrMode == solverHookMode::invalidQuery ? realT( 0 ) : realT( std::max( 1, 32 * n ) );
        integerWork[0] = syevrMode == solverHookMode::invalidIntegerQuery ? 0 : std::max<MXLAPACK_INT>( 1, 10 * n );
        return 0;
    }
    if( syevrMode == solverHookMode::solveFailure )
    {
        return 71;
    }

    *found = syevrMode == solverHookMode::countMismatch ? std::max<MXLAPACK_INT>( 0, n - 1 ) : n;
    for( MXLAPACK_INT column = 0; column < n; ++column )
    {
        eigenvalues[column] = static_cast<realT>( column + 1 );
        for( MXLAPACK_INT row = 0; row < n; ++row )
        {
            eigenvectors[row + column * ldz] = row == column ? realT( 1 ) : realT( 0 );
        }
    }
    if( syevrMode == solverHookMode::nonFiniteValue )
    {
        eigenvalues[0] = std::numeric_limits<realT>::infinity();
    }
    else if( syevrMode == solverHookMode::nonFiniteVector )
    {
        eigenvectors[0] = std::numeric_limits<realT>::infinity();
    }
    else if( syevrMode == solverHookMode::invalidOrdering && n > 1 )
    {
        eigenvalues[0] = realT( 2 );
        eigenvalues[1] = realT( 1 );
    }
    else if( syevrMode == solverHookMode::roundoffClamp )
    {
        eigenvalues[0] = -std::numeric_limits<realT>::epsilon();
    }
    else if( syevrMode == solverHookMode::indefinite )
    {
        eigenvalues[0] = realT( -0.25 );
    }
    return 0;
}

// Emulate selected GESVD query, failure, and malformed-output paths.
MXLAPACK_INT gesvdHook( char jobu,             /* [in] left-vector request */
                        char jobvt,            /* [in] right-vector request */
                        MXLAPACK_INT rows,     /* [in] matrix row count */
                        MXLAPACK_INT cols,     /* [in] matrix column count */
                        realT *matrix,         /* [in,out] input matrix */
                        MXLAPACK_INT lda,      /* [in] input leading dimension */
                        realT *singular,       /* [out] singular values */
                        realT *left,           /* [out] left vectors */
                        MXLAPACK_INT ldu,      /* [in] left-vector leading dimension */
                        realT *rightTranspose, /* [out] transposed right vectors */
                        MXLAPACK_INT ldvt,     /* [in] right-vector leading dimension */
                        realT *work,           /* [in,out] floating workspace */
                        MXLAPACK_INT lwork /* [in] floating workspace size */ )
{
    static_cast<void>( jobu );
    static_cast<void>( jobvt );
    static_cast<void>( rows );
    static_cast<void>( matrix );
    static_cast<void>( lda );
    static_cast<void>( left );
    static_cast<void>( ldu );

    if( lwork == -1 )
    {
        if( gesvdMode == solverHookMode::queryFailure )
        {
            return 62;
        }
        work[0] = gesvdMode == solverHookMode::invalidQuery ? realT( 0 ) : realT( std::max( 1, 5 * cols ) );
        return 0;
    }
    if( gesvdMode == solverHookMode::solveFailure )
    {
        return 72;
    }

    for( MXLAPACK_INT index = 0; index < cols; ++index )
    {
        singular[index] = static_cast<realT>( cols - index );
        for( MXLAPACK_INT row = 0; row < cols; ++row )
        {
            rightTranspose[row + index * ldvt] = row == index ? realT( 1 ) : realT( 0 );
        }
    }
    if( gesvdMode == solverHookMode::nonFiniteValue )
    {
        singular[0] = std::numeric_limits<realT>::infinity();
    }
    else if( gesvdMode == solverHookMode::nonFiniteVector )
    {
        rightTranspose[0] = std::numeric_limits<realT>::infinity();
    }
    else if( gesvdMode == solverHookMode::invalidOrdering && cols > 1 )
    {
        singular[0] = realT( 1 );
        singular[1] = realT( 2 );
    }
    else if( gesvdMode == solverHookMode::negativeSpectrum )
    {
        singular[cols - 1] = realT( -1 );
    }
    else if( gesvdMode == solverHookMode::tinySpectrum )
    {
        singular[0] = realT( 1e-200 );
    }
    return 0;
}

// Emulate selected LAED9 failure and malformed-output paths.
MXLAPACK_INT laed9Hook( realT *eigenvalues,            /* [out] ascending updated eigenvalues */
                        realT *delta,                  /* [out] secular workspace */
                        realT *eigenvectors,           /* [out] updated eigenvectors */
                        MXLAPACK_INT rank,             /* [in] active secular-system dimension */
                        MXLAPACK_INT leadingDimension, /* [in] output leading dimension */
                        realT rho,                     /* [in] positive rank-one weight */
                        realT *poles,                  /* [in,out] ascending diagonal poles */
                        realT *update /* [in,out] normalized update */ )
{
    static_cast<void>( delta );
    static_cast<void>( rho );
    static_cast<void>( update );

    if( laed9Mode == solverHookMode::solveFailure )
    {
        return 73;
    }

    for( MXLAPACK_INT column = 0; column < rank; ++column )
    {
        eigenvalues[column] = poles[column];
        for( MXLAPACK_INT row = 0; row < rank; ++row )
        {
            eigenvectors[row + column * leadingDimension] = row == column ? realT( 1 ) : realT( 0 );
        }
    }
    if( laed9Mode == solverHookMode::nonFiniteValue )
    {
        eigenvalues[0] = std::numeric_limits<realT>::infinity();
    }
    else if( laed9Mode == solverHookMode::nonFiniteVector )
    {
        eigenvectors[0] = std::numeric_limits<realT>::infinity();
    }
    else if( laed9Mode == solverHookMode::invalidOrdering && rank > 1 )
    {
        eigenvalues[0] = realT( 1 );
        eigenvalues[1] = realT( 0 );
    }
    else if( laed9Mode == solverHookMode::outsideInterlacing )
    {
        eigenvalues[0] = poles[0] - realT( 1 );
    }
    else if( laed9Mode == solverHookMode::invalidVectorNorm )
    {
        eigenvectors[0] = realT( 2 );
    }
    else if( laed9Mode == solverHookMode::roundoffClamp )
    {
        eigenvalues[rank - 1] = std::numeric_limits<realT>::epsilon();
    }
    return 0;
}

// Restore production failure hooks at scope exit.
class hookGuard
{
  public:
    // Install no hooks at the beginning of a test scope.
    hookGuard()
    {
        mx::math::detail::svdDeletionHooks<double>() = {};
        mx::math::detail::svdDeletionHooks<float>() = {};
        syevrMode = solverHookMode::production;
        gesvdMode = solverHookMode::production;
        laed9Mode = solverHookMode::production;
    }

    // Restore production hooks after every test exit path.
    ~hookGuard()
    {
        mx::math::detail::svdDeletionHooks<double>() = {};
        mx::math::detail::svdDeletionHooks<float>() = {};
        syevrMode = solverHookMode::production;
        gesvdMode = solverHookMode::production;
        laed9Mode = solverHookMode::production;
    }
};

} // namespace
/** \endcond */

namespace unitTest::math_svdDowndate_test
{

/// SVD row deletion matches direct full SVDs
/** Verifies that svdRemoveRows reproduces direct full SVDs for complete tall, wide, and square factors.
 *
 * \ingroup svdDowndate_unit_tests
 */
TEST_CASE( "SVD row deletion matches direct full SVDs", "[math::svdDowndate][rows]" )
{
    for( const auto [rows, cols, deleted] :
         std::vector<std::tuple<Eigen::Index, Eigen::Index, std::vector<Eigen::Index>>>{ { 7, 4, { 0, 5 } },
                                                                                         { 4, 7, { 1 } },
                                                                                         { 5, 5, { 2 } },
                                                                                         { 9, 3, { 0, 2, 4, 6 } } } )
    {
        const Eigen::Index rank = std::min( rows, cols );
        const matrixT left = sineFactor( rows, rank );
        const matrixT right = sineFactor( cols, rank );
        vectorT singular( rank );
        for( Eigen::Index index = 0; index < rank; ++index )
        {
            singular( index ) = realT( 9 - 2 * index ) + realT( 0.25 * index );
        }
        const matrixT matrix = representedMatrix( left, singular, right );

        for( const backendT backend : { backendT::leadingCovariance, backendT::stableCore } )
        {
            CAPTURE( rows, cols, deleted );
            INFO( "backend: " << mx::math::svdDeletionBackendName( backend ) );
            mx::math::svdDeletionResult<realT> result;
            mx::math::svdDeletionWorkspace<realT> workspace;
            const statusT status = mx::math::svdRemoveRows( result, singular, left, deleted, rank, workspace, backend );
            REQUIRE( mx::math::svdDeletionSucceeded( status ) );
            REQUIRE( result.backend() == backend );
            const deletionComparison comparison = compareRowResult( matrix, right, deleted, result );
            REQUIRE( comparison.squaredSingularError <= comparison.squaredSingularTolerance );
            REQUIRE( comparison.singularError <= comparison.singularTolerance );
            REQUIRE( comparison.covarianceError <= comparison.covarianceTolerance );
        }
    }
}

/// SVD column deletion matches direct full SVDs
/** Verifies that svdRemoveColumns reproduces direct full SVDs and the row-deletion transpose dual.
 *
 * \ingroup svdDowndate_unit_tests
 */
TEST_CASE( "SVD column deletion matches direct full SVDs", "[math::svdDowndate][columns]" )
{
    for( const auto [rows, cols, deleted] :
         std::vector<std::tuple<Eigen::Index, Eigen::Index, std::vector<Eigen::Index>>>{ { 7, 4, { 1 } },
                                                                                         { 4, 7, { 0, 6 } },
                                                                                         { 5, 5, { 1, 3 } } } )
    {
        const Eigen::Index rank = std::min( rows, cols );
        const matrixT left = sineFactor( rows, rank );
        const matrixT right = sineFactor( cols, rank );
        vectorT singular( rank );
        for( Eigen::Index index = 0; index < rank; ++index )
        {
            singular( index ) = realT( 11 - 2 * index ) + realT( 0.125 * index );
        }
        const matrixT matrix = representedMatrix( left, singular, right );

        for( const backendT backend : { backendT::leadingCovariance, backendT::stableCore } )
        {
            CAPTURE( rows, cols, deleted );
            INFO( "backend: " << mx::math::svdDeletionBackendName( backend ) );
            mx::math::svdDeletionResult<realT> result;
            mx::math::svdDeletionWorkspace<realT> workspace;
            const statusT status =
                mx::math::svdRemoveColumns( result, singular, right, deleted, rank, workspace, backend );
            REQUIRE( mx::math::svdDeletionSucceeded( status ) );
            const deletionComparison comparison = compareColumnResult( matrix, left, deleted, result );
            REQUIRE( comparison.squaredSingularError <= comparison.squaredSingularTolerance );
            REQUIRE( comparison.singularError <= comparison.singularTolerance );
            REQUIRE( comparison.covarianceError <= comparison.covarianceTolerance );
        }
    }
}

/// SVD deletion core entry points agree
/** Verifies svdDeletionLeadingCore, svdDeletionStableCore, and svdDeletionCore against the same direct SVD.
 *
 * \ingroup svdDowndate_unit_tests
 */
TEST_CASE( "SVD deletion core entry points agree", "[math::svdDowndate][core]" )
{
    const matrixT left = sineFactor( 7, 4 );
    const matrixT right = sineFactor( 5, 4 );
    vectorT singular( 4 );
    singular << 10, 6, 3, 0.5;
    const matrixT matrix = representedMatrix( left, singular, right );
    const std::vector<Eigen::Index> deleted{ 1, 5 };
    matrixT deletedRows( deleted.size(), left.cols() );
    for( Eigen::Index row = 0; row < deletedRows.rows(); ++row )
    {
        deletedRows.matrix().row( row ) = left.matrix().row( deleted[row] );
    }

    mx::math::svdDeletionResult<realT> leadingResult;
    mx::math::svdDeletionWorkspace<realT> leadingWorkspace;
    REQUIRE( mx::math::svdDeletionSucceeded(
        mx::math::svdDeletionLeadingCore( leadingResult, singular, deletedRows, 4, leadingWorkspace ) ) );
    const deletionComparison leadingComparison = compareRowResult( matrix, right, deleted, leadingResult );
    REQUIRE( leadingComparison.squaredSingularError <= leadingComparison.squaredSingularTolerance );
    REQUIRE( leadingComparison.singularError <= leadingComparison.singularTolerance );
    REQUIRE( leadingComparison.covarianceError <= leadingComparison.covarianceTolerance );

    mx::math::svdDeletionResult<realT> stableResult;
    mx::math::svdDeletionWorkspace<realT> stableWorkspace;
    REQUIRE( mx::math::svdDeletionSucceeded(
        mx::math::svdDeletionStableCore( stableResult, singular, deletedRows, 4, stableWorkspace ) ) );
    const deletionComparison stableComparison = compareRowResult( matrix, right, deleted, stableResult );
    REQUIRE( stableComparison.squaredSingularError <= stableComparison.squaredSingularTolerance );
    REQUIRE( stableComparison.singularError <= stableComparison.singularTolerance );
    REQUIRE( stableComparison.covarianceError <= stableComparison.covarianceTolerance );

    mx::math::svdDeletionResult<realT> dispatchResult;
    REQUIRE( mx::math::svdDeletionSucceeded( mx::math::svdDeletionCore( dispatchResult,
                                                                        singular,
                                                                        deletedRows,
                                                                        3,
                                                                        stableWorkspace,
                                                                        backendT::stableCore ) ) );
    REQUIRE( dispatchResult.outputRank() == 3 );
    REQUIRE( ( dispatchResult.singularValues() - stableResult.singularValues() ).matrix().norm() < 1e-12 );
    REQUIRE( ( dispatchResult.squaredSingularValues() - stableResult.squaredSingularValues() ).matrix().norm() <
             1e-12 );
    REQUIRE( dispatchResult.singularValues().size() == 4 );
    REQUIRE( dispatchResult.rotation().cols() == 3 );

    const matrixT retained = retainedRows( matrix, deleted );
    Eigen::JacobiSVD<Eigen::Matrix<realT, Eigen::Dynamic, Eigen::Dynamic>> directSvd( retained.matrix(),
                                                                                      Eigen::ComputeThinV );
    const Eigen::Matrix<realT, Eigen::Dynamic, Eigen::Dynamic> directDirections = directSvd.matrixV().leftCols( 3 );
    const Eigen::Matrix<realT, Eigen::Dynamic, Eigen::Dynamic> updatedDirections =
        right.matrix() * dispatchResult.rotation().matrix();
    const Eigen::Matrix<realT, Eigen::Dynamic, Eigen::Dynamic> directProjector =
        directDirections * directDirections.transpose();
    const Eigen::Matrix<realT, Eigen::Dynamic, Eigen::Dynamic> updatedProjector =
        updatedDirections * updatedDirections.transpose();
    REQUIRE( ( updatedProjector - directProjector ).norm() < 1e-10 );
}

/// Rank-one secular SVD deletion matches direct row and column SVDs
/** Verifies that svdRemoveRows and svdRemoveColumns with the rankOneSecular backend reproduce direct full SVDs
 * while publishing either the complete eigensystem or a requested leading prefix.
 *
 * \ingroup svdDowndate_unit_tests
 */
TEST_CASE( "Rank-one secular SVD deletion matches direct row and column SVDs", "[math::svdDowndate][rankOneSecular]" )
{
    const Eigen::Index rank{ 4 };
    const matrixT left = sineFactor( 9, rank );
    const matrixT right = sineFactor( 7, rank );
    vectorT singular( rank );
    singular << 11, 6, 2.5, 0.75;
    const matrixT matrix = representedMatrix( left, singular, right );

    SECTION( "row deletion" )
    {
        const std::vector<Eigen::Index> deleted{ 4 };
        mx::math::svdDeletionResult<realT> result;
        mx::math::svdDeletionWorkspace<realT> workspace;
        REQUIRE(
            mx::math::svdRemoveRows( result, singular, left, deleted, rank, workspace, backendT::rankOneSecular ) ==
            statusT::success );
        REQUIRE( result.backend() == backendT::rankOneSecular );
        const deletionComparison comparison = compareRowResult( matrix, right, deleted, result );
        REQUIRE( comparison.squaredSingularError <= comparison.squaredSingularTolerance );
        REQUIRE( comparison.singularError <= comparison.singularTolerance );
        REQUIRE( comparison.covarianceError <= comparison.covarianceTolerance );
    }

    SECTION( "column deletion" )
    {
        const std::vector<Eigen::Index> deleted{ 2 };
        mx::math::svdDeletionResult<realT> result;
        mx::math::svdDeletionWorkspace<realT> workspace;
        REQUIRE(
            mx::math::svdRemoveColumns( result, singular, right, deleted, rank, workspace, backendT::rankOneSecular ) ==
            statusT::success );
        REQUIRE( result.backend() == backendT::rankOneSecular );
        const deletionComparison comparison = compareColumnResult( matrix, left, deleted, result );
        REQUIRE( comparison.squaredSingularError <= comparison.squaredSingularTolerance );
        REQUIRE( comparison.singularError <= comparison.singularTolerance );
        REQUIRE( comparison.covarianceError <= comparison.covarianceTolerance );
    }

    SECTION( "leading output prefix" )
    {
        constexpr Eigen::Index outputRank{ 2 };
        const std::vector<Eigen::Index> deleted{ 4 };
        const matrixT retained = retainedRows( matrix, deleted );
        Eigen::JacobiSVD<Eigen::Matrix<realT, Eigen::Dynamic, Eigen::Dynamic>> direct( retained.matrix(),
                                                                                       Eigen::ComputeThinV );

        mx::math::svdDeletionResult<realT> result;
        mx::math::svdDeletionWorkspace<realT> workspace;
        REQUIRE( mx::math::svdRemoveRows( result,
                                          singular,
                                          left,
                                          deleted,
                                          outputRank,
                                          workspace,
                                          backendT::rankOneSecular ) == statusT::success );
        REQUIRE( result.outputRank() == outputRank );
        REQUIRE( result.rotation().rows() == rank );
        REQUIRE( result.rotation().cols() == outputRank );
        REQUIRE( result.singularValues().size() == rank );
        REQUIRE( ( result.singularValues().head( outputRank ).matrix() - direct.singularValues().head( outputRank ) )
                     .norm() < 5e-11 );

        const Eigen::Matrix<realT, Eigen::Dynamic, Eigen::Dynamic> expectedDirections =
            direct.matrixV().leftCols( outputRank );
        const Eigen::Matrix<realT, Eigen::Dynamic, Eigen::Dynamic> actualDirections =
            right.matrix() * result.rotation().matrix();
        const Eigen::Matrix<realT, Eigen::Dynamic, Eigen::Dynamic> expectedProjector =
            expectedDirections * expectedDirections.transpose();
        const Eigen::Matrix<realT, Eigen::Dynamic, Eigen::Dynamic> actualProjector =
            actualDirections * actualDirections.transpose();
        REQUIRE( ( actualProjector - expectedProjector ).norm() < 5e-10 );
    }
}

/// Rank-one secular SVD deletion handles deflation and leverage edge cases
/** Verifies exact and clustered repeated singular values, zero singular values, zero-leverage deletion, and
 * high-leverage deletion through svdRemoveRows with the rankOneSecular backend.
 *
 * \ingroup svdDowndate_unit_tests
 */
TEST_CASE( "Rank-one secular SVD deletion handles deflation and leverage edge cases",
           "[math::svdDowndate][rankOneSecular][conditioning]" )
{
    SECTION( "scalar secular system" )
    {
        matrixT left( 2, 1 );
        left << std::sqrt( 0.5 ), std::sqrt( 0.5 );
        vectorT singular( 1 );
        singular << std::sqrt( 2.0 );
        const std::vector<Eigen::Index> deleted{ 0 };

        mx::math::svdDeletionResult<realT> result;
        mx::math::svdDeletionWorkspace<realT> workspace;
        REQUIRE( mx::math::svdRemoveRows( result, singular, left, deleted, 1, workspace, backendT::rankOneSecular ) ==
                 statusT::success );
        REQUIRE( result.singularValues()( 0 ) == Approx( 1.0 ).epsilon( 1e-12 ) );
        REQUIRE( result.squaredSingularValues()( 0 ) == Approx( 1.0 ).epsilon( 1e-12 ) );
        REQUIRE( std::abs( result.rotation()( 0, 0 ) ) == Approx( 1.0 ).epsilon( 1e-12 ) );
    }

    SECTION( "exact repeated spectrum" )
    {
        const matrixT left = sineFactor( 8, 4 );
        const matrixT right = sineFactor( 6, 4 );
        vectorT singular( 4 );
        singular << 9, 9, 3, 3;
        const matrixT matrix = representedMatrix( left, singular, right );
        const std::vector<Eigen::Index> deleted{ 3 };

        mx::math::svdDeletionResult<realT> result;
        mx::math::svdDeletionWorkspace<realT> workspace;
        REQUIRE( mx::math::svdDeletionSucceeded(
            mx::math::svdRemoveRows( result, singular, left, deleted, 4, workspace, backendT::rankOneSecular ) ) );
        const deletionComparison comparison = compareRowResult( matrix, right, deleted, result, 2e-10 );
        REQUIRE( comparison.squaredSingularError <= comparison.squaredSingularTolerance );
        REQUIRE( comparison.singularError <= comparison.singularTolerance );
        REQUIRE( comparison.covarianceError <= comparison.covarianceTolerance );
    }

    SECTION( "clustered spectrum" )
    {
        const matrixT left = sineFactor( 9, 4 );
        const matrixT right = sineFactor( 7, 4 );
        const realT spacing = realT( 8 ) * std::numeric_limits<realT>::epsilon();
        vectorT singular( 4 );
        singular << 10, 10 * ( 1 - spacing ), 2, 2 * ( 1 - spacing );
        const matrixT matrix = representedMatrix( left, singular, right );
        const std::vector<Eigen::Index> deleted{ 5 };

        mx::math::svdDeletionResult<realT> result;
        mx::math::svdDeletionWorkspace<realT> workspace;
        REQUIRE( mx::math::svdDeletionSucceeded(
            mx::math::svdRemoveRows( result, singular, left, deleted, 4, workspace, backendT::rankOneSecular ) ) );
        const deletionComparison comparison = compareRowResult( matrix, right, deleted, result, 5e-10 );
        REQUIRE( comparison.squaredSingularError <= comparison.squaredSingularTolerance );
        REQUIRE( comparison.singularError <= comparison.singularTolerance );
        REQUIRE( comparison.covarianceError <= comparison.covarianceTolerance );
    }

    SECTION( "zero singular spectrum" )
    {
        const matrixT left = sineFactor( 8, 4 );
        const matrixT right = sineFactor( 6, 4 );
        vectorT singular( 4 );
        singular << 8, 3, 0, 0;
        const matrixT matrix = representedMatrix( left, singular, right );
        const std::vector<Eigen::Index> deleted{ 2 };

        mx::math::svdDeletionResult<realT> result;
        mx::math::svdDeletionWorkspace<realT> workspace;
        REQUIRE( mx::math::svdDeletionSucceeded(
            mx::math::svdRemoveRows( result, singular, left, deleted, 4, workspace, backendT::rankOneSecular ) ) );
        const deletionComparison comparison = compareRowResult( matrix, right, deleted, result, 2e-10 );
        REQUIRE( comparison.squaredSingularError <= comparison.squaredSingularTolerance );
        REQUIRE( comparison.singularError <= comparison.singularTolerance );
        REQUIRE( comparison.covarianceError <= comparison.covarianceTolerance );
    }

    SECTION( "zero leverage" )
    {
        matrixT left( 5, 4 );
        left.setZero();
        left.matrix().topRows( 4 ).setIdentity();
        const matrixT right = sineFactor( 6, 4 );
        vectorT singular( 4 );
        singular << 8, 4, 2, 1;
        const matrixT matrix = representedMatrix( left, singular, right );
        const std::vector<Eigen::Index> deleted{ 4 };

        mx::math::svdDeletionResult<realT> result;
        mx::math::svdDeletionWorkspace<realT> workspace;
        REQUIRE( mx::math::svdDeletionSucceeded(
            mx::math::svdRemoveRows( result, singular, left, deleted, 4, workspace, backendT::rankOneSecular ) ) );
        const deletionComparison comparison = compareRowResult( matrix, right, deleted, result );
        REQUIRE( comparison.squaredSingularError <= comparison.squaredSingularTolerance );
        REQUIRE( comparison.singularError <= comparison.singularTolerance );
        REQUIRE( comparison.covarianceError <= comparison.covarianceTolerance );
    }

    SECTION( "high leverage" )
    {
        matrixT left( 6, 3 );
        left.setZero();
        const realT residualLeverage = 1e-10;
        left( 0, 0 ) = std::sqrt( 1 - residualLeverage );
        left( 3, 0 ) = std::sqrt( residualLeverage );
        left( 1, 1 ) = std::sqrt( 0.75 );
        left( 4, 1 ) = 0.5;
        left( 2, 2 ) = std::sqrt( 0.6 );
        left( 5, 2 ) = std::sqrt( 0.4 );
        const matrixT right = sineFactor( 5, 3 );
        vectorT singular( 3 );
        singular << 10, 4, 1;
        const matrixT matrix = representedMatrix( left, singular, right );
        const std::vector<Eigen::Index> deleted{ 0 };

        REQUIRE( mx::math::validateSvdDeletionFactor( left ) == statusT::success );
        mx::math::svdDeletionResult<realT> result;
        mx::math::svdDeletionWorkspace<realT> workspace;
        REQUIRE( mx::math::svdDeletionSucceeded(
            mx::math::svdRemoveRows( result, singular, left, deleted, 3, workspace, backendT::rankOneSecular ) ) );
        const deletionComparison comparison = compareRowResult( matrix, right, deleted, result, 5e-10 );
        REQUIRE( comparison.squaredSingularError <= comparison.squaredSingularTolerance );
        REQUIRE( comparison.singularError <= comparison.singularTolerance );
        REQUIRE( comparison.covarianceError <= comparison.covarianceTolerance );
    }
}

/// Rank-one secular SVD deletion supports float and enforces its deletion contract
/** Verifies float accuracy, the empty identity path, one-row-only rejection in svdDeletionCore and svdRemoveRows,
 * workspace capacity validation, and stable backend/status names for rankOneSecular.
 *
 * \ingroup svdDowndate_unit_tests
 */
TEST_CASE( "Rank-one secular SVD deletion supports float and enforces its deletion contract",
           "[math::svdDowndate][rankOneSecular][float][errors]" )
{
    hookGuard guard;
    REQUIRE( std::string( mx::math::svdDeletionBackendName( backendT::rankOneSecular ) ) == "rankOneSecular" );
    REQUIRE( std::string( mx::math::svdDeletionStatusName( statusT::unsupportedDeletionCount ) ) ==
             "unsupportedDeletionCount" );

    SECTION( "float row deletion" )
    {
        using floatMatrixT = mx::math::svdDeletionMatrix<float>;
        using floatVectorT = mx::math::svdDeletionVector<float>;
        const floatMatrixT left = sineFactor( 7, 3 ).cast<float>();
        const floatMatrixT right = sineFactor( 5, 3 ).cast<float>();
        floatVectorT singular( 3 );
        singular << 7, 3, 0.5F;
        const Eigen::MatrixXf matrix = left.matrix() * singular.matrix().asDiagonal() * right.matrix().transpose();
        Eigen::MatrixXf retained( 6, 5 );
        retained.topRows( 2 ) = matrix.topRows( 2 );
        retained.bottomRows( 4 ) = matrix.bottomRows( 4 );
        Eigen::JacobiSVD<Eigen::MatrixXf> direct( retained, Eigen::ComputeThinV );
        const std::vector<Eigen::Index> deleted{ 2 };

        mx::math::svdDeletionResult<float> result;
        mx::math::svdDeletionWorkspace<float> workspace;
        REQUIRE( mx::math::svdDeletionSucceeded(
            mx::math::svdRemoveRows( result, singular, left, deleted, 3, workspace, backendT::rankOneSecular ) ) );
        REQUIRE( result.backend() == backendT::rankOneSecular );
        REQUIRE( ( result.singularValues().matrix() - direct.singularValues().head( 3 ) ).norm() < 5e-4F );
    }

    SECTION( "empty deletion is identity" )
    {
        vectorT singular( 3 );
        singular << 7, 2, 0.5;
        const matrixT factor = identityMatrix( 3 );
        const std::vector<Eigen::Index> deleted;
        mx::math::svdDeletionResult<realT> result;
        mx::math::svdDeletionWorkspace<realT> workspace;

        REQUIRE( mx::math::svdRemoveRows( result, singular, factor, deleted, 2, workspace, backendT::rankOneSecular ) ==
                 statusT::success );
        REQUIRE( result.backend() == backendT::rankOneSecular );
        REQUIRE( result.outputRank() == 2 );
        REQUIRE( ( result.singularValues() - singular ).matrix().norm() == Approx( 0.0 ) );
        REQUIRE( ( result.rotation() - identityMatrix( 3 ).leftCols( 2 ) ).matrix().norm() == Approx( 0.0 ) );
        REQUIRE( result.minimumPSDValue() == Approx( 0.25 / 49.0 ) );
    }

    SECTION( "direct empty deletion is identity" )
    {
        vectorT singular( 3 );
        singular << 7, 2, 0.5;
        matrixT deletedRows( 0, 3 );
        mx::math::svdDeletionResult<realT> result;
        mx::math::svdDeletionWorkspace<realT> workspace;

        REQUIRE( mx::math::svdDeletionCore( result, singular, deletedRows, 2, workspace, backendT::rankOneSecular ) ==
                 statusT::success );
        REQUIRE( result.backend() == backendT::rankOneSecular );
        REQUIRE( ( result.singularValues() - singular ).matrix().norm() == Approx( 0.0 ) );
        REQUIRE( ( result.rotation() - identityMatrix( 3 ).leftCols( 2 ) ).matrix().norm() == Approx( 0.0 ) );
    }

    SECTION( "zero spectrum returns an arbitrary identity basis" )
    {
        vectorT singular = vectorT::Zero( 2 );
        matrixT deletedRows( 1, 2 );
        deletedRows << 0.25, 0.5;
        mx::math::svdDeletionResult<realT> result;
        mx::math::svdDeletionWorkspace<realT> workspace;

        REQUIRE( mx::math::svdDeletionCore( result, singular, deletedRows, 2, workspace, backendT::rankOneSecular ) ==
                 statusT::success );
        REQUIRE( result.singularValues().matrix().norm() == Approx( 0.0 ) );
        REQUIRE( ( result.rotation() - identityMatrix( 2 ) ).matrix().norm() == Approx( 0.0 ) );
    }

    SECTION( "tiny update is treated as identity" )
    {
        vectorT singular( 2 );
        singular << 1, 0.5;
        matrixT deletedRows( 1, 2 );
        deletedRows << 1e-200, 0;
        mx::math::svdDeletionResult<realT> result;
        mx::math::svdDeletionWorkspace<realT> workspace;

        REQUIRE( mx::math::svdDeletionCore( result, singular, deletedRows, 2, workspace, backendT::rankOneSecular ) ==
                 statusT::success );
        REQUIRE( ( result.singularValues() - singular ).matrix().norm() == Approx( 0.0 ) );
        REQUIRE( ( result.rotation() - identityMatrix( 2 ) ).matrix().norm() == Approx( 0.0 ) );
    }

    SECTION( "invalid direct core shape is rejected" )
    {
        vectorT singular( 2 );
        singular << 3, 1;
        matrixT deletedRows( 1, 1 );
        deletedRows << 0.25;
        mx::math::svdDeletionResult<realT> result;
        mx::math::svdDeletionWorkspace<realT> workspace;

        REQUIRE( mx::math::svdDeletionCore( result, singular, deletedRows, 2, workspace, backendT::rankOneSecular ) ==
                 statusT::invalidInput );
        REQUIRE( result.status() == statusT::invalidInput );
    }

    SECTION( "result preparation failure is propagated" )
    {
        vectorT singular( 2 );
        singular << 3, 1;
        matrixT deletedRows( 1, 2 );
        deletedRows << 0.25, 0.1;
        mx::math::svdDeletionResult<realT> result;
        mx::math::svdDeletionWorkspace<realT> workspace;

        failingOperation = mx::math::detail::svdDeletionTestOperation::prepareResult;
        mx::math::detail::svdDeletionHooks<double>().operation = throwAllocation;
        REQUIRE( mx::math::svdDeletionCore( result, singular, deletedRows, 2, workspace, backendT::rankOneSecular ) ==
                 statusT::allocationFailure );
        REQUIRE( result.status() == statusT::allocationFailure );
    }

    SECTION( "multiple deletions are rejected" )
    {
        const matrixT factor = sineFactor( 6, 3 );
        vectorT singular( 3 );
        singular << 7, 3, 1;
        const std::vector<Eigen::Index> deleted{ 1, 4 };
        matrixT deletedRows( 2, 3 );
        deletedRows.matrix().row( 0 ) = factor.matrix().row( 1 );
        deletedRows.matrix().row( 1 ) = factor.matrix().row( 4 );
        mx::math::svdDeletionResult<realT> result;
        mx::math::svdDeletionWorkspace<realT> workspace;

        REQUIRE( workspace.prepare( 3, 2, backendT::rankOneSecular ) == statusT::unsupportedDeletionCount );
        REQUIRE( !workspace.prepared() );
        REQUIRE( mx::math::svdDeletionCore( result, singular, deletedRows, 3, workspace, backendT::rankOneSecular ) ==
                 statusT::unsupportedDeletionCount );
        REQUIRE( result.status() == statusT::unsupportedDeletionCount );
        REQUIRE( mx::math::svdRemoveRows( result, singular, factor, deleted, 3, workspace, backendT::rankOneSecular ) ==
                 statusT::unsupportedDeletionCount );
        REQUIRE( result.status() == statusT::unsupportedDeletionCount );
    }

    SECTION( "finite oversized update is rejected before squaring" )
    {
        vectorT singular( 2 );
        singular << 1, 0.5;
        matrixT deletedRows( 1, 2 );
        deletedRows << std::numeric_limits<realT>::max() / 4, std::numeric_limits<realT>::max() / 4;
        mx::math::svdDeletionResult<realT> result;
        mx::math::svdDeletionWorkspace<realT> workspace;

        REQUIRE( mx::math::svdDeletionCore( result, singular, deletedRows, 2, workspace, backendT::rankOneSecular ) ==
                 statusT::invalidInput );
        REQUIRE( result.status() == statusT::invalidInput );
    }

    SECTION( "materially indefinite core is rejected" )
    {
        vectorT singular( 2 );
        singular << 1, 0.5;
        matrixT deletedRows( 1, 2 );
        deletedRows << 2, 0;
        mx::math::svdDeletionResult<realT> result;
        mx::math::svdDeletionWorkspace<realT> workspace;

        REQUIRE( mx::math::svdDeletionCore( result, singular, deletedRows, 2, workspace, backendT::rankOneSecular ) ==
                 statusT::nonPositiveSemidefinite );
        REQUIRE( result.minimumPSDValue() == Approx( -3.0 ) );
    }

    SECTION( "finite result that cannot be rescaled is rejected" )
    {
        const realT scale = realT( 2 ) * std::sqrt( std::numeric_limits<realT>::max() );
        vectorT singular( 2 );
        singular << scale, scale / realT( 2 );
        matrixT deletedRows( 1, 2 );
        deletedRows << 0.25, 0;
        mx::math::svdDeletionResult<realT> result;
        mx::math::svdDeletionWorkspace<realT> workspace;

        REQUIRE( mx::math::svdDeletionCore( result, singular, deletedRows, 2, workspace, backendT::rankOneSecular ) ==
                 statusT::rescalingOverflow );
        REQUIRE( result.status() == statusT::rescalingOverflow );
    }
}

/// Rank-one secular SVD deletion reports solver failures
/** Verifies LAED9 failure, non-finite output, invalid ordering and interlacing, vector norm and residual validation,
 * and roundoff clamping through svdDeletionCore with the rankOneSecular backend.
 *
 * \ingroup svdDowndate_unit_tests
 */
TEST_CASE( "Rank-one secular SVD deletion reports solver failures",
           "[math::svdDowndate][rankOneSecular][errors][solver]" )
{
    hookGuard guard;
    vectorT singular( 3 );
    singular << 8, 4, 2;
    matrixT deleted( 1, 3 );
    deleted << 0.25, 0.2, 0.1;

    SECTION( "LAED9 solve failure" )
    {
        mx::math::svdDeletionResult<realT> result;
        mx::math::svdDeletionWorkspace<realT> workspace;
        laed9Mode = solverHookMode::solveFailure;
        mx::math::detail::svdDeletionHooks<double>().laed9 = laed9Hook;
        REQUIRE( mx::math::svdDeletionCore( result, singular, deleted, 3, workspace, backendT::rankOneSecular ) ==
                 statusT::solverFailure );
        REQUIRE( result.lapackInfo() == 73 );
    }

    SECTION( "LAED9 non-finite values and vectors" )
    {
        for( const solverHookMode mode : { solverHookMode::nonFiniteValue, solverHookMode::nonFiniteVector } )
        {
            CAPTURE( static_cast<int>( mode ) );
            mx::math::svdDeletionResult<realT> result;
            mx::math::svdDeletionWorkspace<realT> workspace;
            laed9Mode = mode;
            mx::math::detail::svdDeletionHooks<double>().laed9 = laed9Hook;
            REQUIRE( mx::math::svdDeletionCore( result, singular, deleted, 3, workspace, backendT::rankOneSecular ) ==
                     statusT::nonFiniteOutput );
        }
    }

    SECTION( "LAED9 invalid ordering" )
    {
        mx::math::svdDeletionResult<realT> result;
        mx::math::svdDeletionWorkspace<realT> workspace;
        laed9Mode = solverHookMode::invalidOrdering;
        mx::math::detail::svdDeletionHooks<double>().laed9 = laed9Hook;
        REQUIRE( mx::math::svdDeletionCore( result, singular, deleted, 3, workspace, backendT::rankOneSecular ) ==
                 statusT::invalidSolverOutput );
    }

    SECTION( "LAED9 result outside secular interlacing bounds" )
    {
        mx::math::svdDeletionResult<realT> result;
        mx::math::svdDeletionWorkspace<realT> workspace;
        laed9Mode = solverHookMode::outsideInterlacing;
        mx::math::detail::svdDeletionHooks<double>().laed9 = laed9Hook;
        REQUIRE( mx::math::svdDeletionCore( result, singular, deleted, 3, workspace, backendT::rankOneSecular ) ==
                 statusT::invalidSolverOutput );
    }

    SECTION( "LAED9 non-unit eigenvector" )
    {
        mx::math::svdDeletionResult<realT> result;
        mx::math::svdDeletionWorkspace<realT> workspace;
        laed9Mode = solverHookMode::invalidVectorNorm;
        mx::math::detail::svdDeletionHooks<double>().laed9 = laed9Hook;
        REQUIRE( mx::math::svdDeletionCore( result, singular, deleted, 3, workspace, backendT::rankOneSecular ) ==
                 statusT::invalidSolverOutput );
    }

    SECTION( "LAED9 eigenvector with a large matrix-free residual" )
    {
        mx::math::svdDeletionResult<realT> result;
        mx::math::svdDeletionWorkspace<realT> workspace;
        laed9Mode = solverHookMode::invalidResidual;
        mx::math::detail::svdDeletionHooks<double>().laed9 = laed9Hook;
        REQUIRE( mx::math::svdDeletionCore( result, singular, deleted, 3, workspace, backendT::rankOneSecular ) ==
                 statusT::invalidSolverOutput );
    }

    SECTION( "small positive roundoff eigenvalue is clamped" )
    {
        vectorT clampSingular( 2 );
        clampSingular << 1, 0.5;
        matrixT clampDeleted( 1, 2 );
        clampDeleted << 1, 0;
        mx::math::svdDeletionResult<realT> result;
        mx::math::svdDeletionWorkspace<realT> workspace;
        laed9Mode = solverHookMode::roundoffClamp;
        mx::math::detail::svdDeletionHooks<double>().laed9 = laed9Hook;

        REQUIRE(
            mx::math::svdDeletionCore( result, clampSingular, clampDeleted, 1, workspace, backendT::rankOneSecular ) ==
            statusT::successWithClamping );
        REQUIRE( result.clampedEigenvalues() == 1 );
        REQUIRE( result.squaredSingularValues()( 1 ) == Approx( 0.0 ) );
    }
}

/// SVD deletion preserves the projected-factor contract
/** Verifies that truncated-factor deletion matches a direct SVD of the represented low-rank matrix, not the
 * original full-rank matrix, through svdRemoveRows and svdRemoveColumns.
 *
 * \ingroup svdDowndate_unit_tests
 */
TEST_CASE( "SVD deletion preserves the projected-factor contract", "[math::svdDowndate][projected]" )
{
    const matrixT fullLeft = sineFactor( 7, 5 );
    const matrixT fullRight = sineFactor( 5, 5 );
    vectorT fullSingular( 5 );
    fullSingular << 10, 8, 6, 5, 4;
    const matrixT full = representedMatrix( fullLeft, fullSingular, fullRight );

    const Eigen::Index rank = 3;
    const matrixT left = fullLeft.matrix().leftCols( rank ).array();
    const matrixT right = fullRight.matrix().leftCols( rank ).array();
    const vectorT singular = fullSingular.head( rank );
    const matrixT projected = representedMatrix( left, singular, right );

    const std::vector<Eigen::Index> deletedRows{ 1, 5 };
    mx::math::svdDeletionResult<realT> rowResult;
    mx::math::svdDeletionWorkspace<realT> rowWorkspace;
    REQUIRE( mx::math::svdDeletionSucceeded(
        mx::math::svdRemoveRows( rowResult, singular, left, deletedRows, rank, rowWorkspace, backendT::stableCore ) ) );
    const deletionComparison rowComparison = compareRowResult( projected, right, deletedRows, rowResult );
    REQUIRE( rowComparison.squaredSingularError <= rowComparison.squaredSingularTolerance );
    REQUIRE( rowComparison.singularError <= rowComparison.singularTolerance );
    REQUIRE( rowComparison.covarianceError <= rowComparison.covarianceTolerance );
    const vectorT projectedRowValues = directSingularValues( retainedRows( projected, deletedRows ), rank );
    const vectorT fullRowValues = directSingularValues( retainedRows( full, deletedRows ), rank );
    REQUIRE( ( projectedRowValues - fullRowValues ).matrix().norm() > 0.1 );

    const std::vector<Eigen::Index> deletedColumns{ 0, 4 };
    mx::math::svdDeletionResult<realT> columnResult;
    mx::math::svdDeletionWorkspace<realT> columnWorkspace;
    REQUIRE( mx::math::svdDeletionSucceeded( mx::math::svdRemoveColumns( columnResult,
                                                                         singular,
                                                                         right,
                                                                         deletedColumns,
                                                                         rank,
                                                                         columnWorkspace,
                                                                         backendT::leadingCovariance ) ) );
    const deletionComparison columnComparison = compareColumnResult( projected, left, deletedColumns, columnResult );
    REQUIRE( columnComparison.squaredSingularError <= columnComparison.squaredSingularTolerance );
    REQUIRE( columnComparison.singularError <= columnComparison.singularTolerance );
    REQUIRE( columnComparison.covarianceError <= columnComparison.covarianceTolerance );
    const vectorT projectedColumnValues = directSingularValues( retainedColumns( projected, deletedColumns ), rank );
    const vectorT fullColumnValues = directSingularValues( retainedColumns( full, deletedColumns ), rank );
    REQUIRE( ( projectedColumnValues - fullColumnValues ).matrix().norm() > 0.1 );
}

/// SVD deletion handles repeated and high-leverage systems
/** Verifies repeated spectra, high-leverage deletion, and the one-row exact deletion invariant through
 * svdRemoveRows without comparing ambiguous individual singular vectors.
 *
 * \ingroup svdDowndate_unit_tests
 */
TEST_CASE( "SVD deletion handles repeated and high-leverage systems", "[math::svdDowndate][conditioning]" )
{
    SECTION( "repeated spectrum" )
    {
        const matrixT left = sineFactor( 8, 4 );
        const matrixT right = sineFactor( 6, 4 );
        vectorT singular( 4 );
        singular << 9, 9, 3, 3;
        const matrixT matrix = representedMatrix( left, singular, right );
        const std::vector<Eigen::Index> deleted{ 0, 4, 7 };

        for( const backendT backend : { backendT::leadingCovariance, backendT::stableCore } )
        {
            INFO( "backend: " << mx::math::svdDeletionBackendName( backend ) );
            mx::math::svdDeletionResult<realT> result;
            mx::math::svdDeletionWorkspace<realT> workspace;
            REQUIRE( mx::math::svdDeletionSucceeded(
                mx::math::svdRemoveRows( result, singular, left, deleted, 4, workspace, backend ) ) );
            const deletionComparison comparison = compareRowResult( matrix, right, deleted, result );
            REQUIRE( comparison.squaredSingularError <= comparison.squaredSingularTolerance );
            REQUIRE( comparison.singularError <= comparison.singularTolerance );
            REQUIRE( comparison.covarianceError <= comparison.covarianceTolerance );
        }
    }

    SECTION( "high leverage" )
    {
        matrixT left( 6, 3 );
        left.setZero();
        const realT residualLeverage = 1e-8;
        left( 0, 0 ) = std::sqrt( 1 - residualLeverage );
        left( 3, 0 ) = std::sqrt( residualLeverage );
        left( 1, 1 ) = std::sqrt( 0.75 );
        left( 4, 1 ) = 0.5;
        left( 2, 2 ) = std::sqrt( 0.6 );
        left( 5, 2 ) = std::sqrt( 0.4 );
        const matrixT right = sineFactor( 5, 3 );
        vectorT singular( 3 );
        singular << 10, 4, 1;
        const matrixT matrix = representedMatrix( left, singular, right );
        const std::vector<Eigen::Index> deleted{ 0, 1 };

        REQUIRE( mx::math::validateSvdDeletionFactor( left ) == statusT::success );
        for( const backendT backend : { backendT::leadingCovariance, backendT::stableCore } )
        {
            INFO( "backend: " << mx::math::svdDeletionBackendName( backend ) );
            mx::math::svdDeletionResult<realT> result;
            mx::math::svdDeletionWorkspace<realT> workspace;
            REQUIRE( mx::math::svdDeletionSucceeded(
                mx::math::svdRemoveRows( result, singular, left, deleted, 3, workspace, backend ) ) );
            const deletionComparison comparison = compareRowResult( matrix, right, deleted, result, 2e-10 );
            REQUIRE( comparison.squaredSingularError <= comparison.squaredSingularTolerance );
            REQUIRE( comparison.singularError <= comparison.singularTolerance );
            REQUIRE( comparison.covarianceError <= comparison.covarianceTolerance );
        }
    }

    SECTION( "one-row deletion" )
    {
        matrixT left( 2, 1 );
        left << std::sqrt( 0.5 ), std::sqrt( 0.5 );
        vectorT singular( 1 );
        singular << std::sqrt( 2.0 );
        const std::vector<Eigen::Index> deleted{ 0 };

        for( const backendT backend : { backendT::leadingCovariance, backendT::stableCore } )
        {
            mx::math::svdDeletionResult<realT> result;
            mx::math::svdDeletionWorkspace<realT> workspace;
            REQUIRE( mx::math::svdDeletionSucceeded(
                mx::math::svdRemoveRows( result, singular, left, deleted, 1, workspace, backend ) ) );
            REQUIRE( result.singularValues()( 0 ) == Approx( 1.0 ).epsilon( 1e-12 ) );
            REQUIRE( result.squaredSingularValues()( 0 ) == Approx( 1.0 ).epsilon( 1e-12 ) );
        }

        const realT literalLongMalesValue = singular( 0 ) * ( realT( 1 ) - left( 0, 0 ) * left( 0, 0 ) );
        REQUIRE( literalLongMalesValue == Approx( realT( 1 ) / std::sqrt( realT( 2 ) ) ).epsilon( 1e-12 ) );
        REQUIRE( std::abs( literalLongMalesValue - realT( 1 ) ) > 0.25 );
    }

    SECTION( "square deleted-side factor" )
    {
        const matrixT left = identityMatrix( 2 );
        vectorT singular( 2 );
        singular << 3, 1;
        const std::vector<Eigen::Index> deleted{ 0 };
        matrixT literalLongMalesCore = identityMatrix( 2 );
        literalLongMalesCore.matrix().noalias() -= left.matrix().row( 0 ).transpose() * left.matrix().row( 0 );
        literalLongMalesCore.matrix() = singular.matrix().asDiagonal() * literalLongMalesCore.matrix();
        const vectorT literalLongMalesValues = directSingularValues( literalLongMalesCore, 2 );

        mx::math::svdDeletionResult<realT> result;
        mx::math::svdDeletionWorkspace<realT> workspace;
        REQUIRE( mx::math::svdDeletionSucceeded(
            mx::math::svdRemoveRows( result, singular, left, deleted, 2, workspace, backendT::stableCore ) ) );
        REQUIRE( ( result.singularValues() - literalLongMalesValues ).matrix().norm() < 1e-12 );
    }
}

/// SVD deletion is invariant across finite scales
/** Verifies exponent-safe scale invariance of svdRemoveRows for very large and very small finite spectra.
 *
 * \ingroup svdDowndate_unit_tests
 */
TEST_CASE( "SVD deletion is invariant across finite scales", "[math::svdDowndate][scaling]" )
{
    const matrixT left = sineFactor( 7, 3 );
    vectorT baseSingular( 3 );
    baseSingular << 8, 3, 1;
    const std::vector<Eigen::Index> deleted{ 1, 5 };

    for( const backendT backend : { backendT::leadingCovariance, backendT::stableCore } )
    {
        mx::math::svdDeletionResult<realT> reference;
        mx::math::svdDeletionWorkspace<realT> referenceWorkspace;
        REQUIRE( mx::math::svdDeletionSucceeded(
            mx::math::svdRemoveRows( reference, baseSingular, left, deleted, 3, referenceWorkspace, backend ) ) );

        for( const realT scale : { realT( 1e140 ), realT( 1e-140 ) } )
        {
            const vectorT scaledSingular = baseSingular * scale;
            mx::math::svdDeletionResult<realT> scaled;
            mx::math::svdDeletionWorkspace<realT> workspace;
            REQUIRE( mx::math::svdDeletionSucceeded(
                mx::math::svdRemoveRows( scaled, scaledSingular, left, deleted, 3, workspace, backend ) ) );
            for( Eigen::Index index = 0; index < 3; ++index )
            {
                REQUIRE( scaled.singularValues()( index ) ==
                         Approx( reference.singularValues()( index ) * scale ).epsilon( 5e-11 ) );
                REQUIRE( scaled.squaredSingularValues()( index ) ==
                         Approx( reference.squaredSingularValues()( index ) * scale * scale ).epsilon( 1e-10 ) );
            }
        }
    }
}

/// SVD deletion handles structural and numerical edge cases
/** Verifies rank loss, repeated spectra, zero singular values, empty deletion, and factor validation through
 * the public SVD deletion APIs.
 *
 * \ingroup svdDowndate_unit_tests
 */
TEST_CASE( "SVD deletion handles structural and numerical edge cases", "[math::svdDowndate][rank]" )
{
    matrixT left( 3, 2 );
    left << 1, 0, 0, std::sqrt( 0.5 ), 0, std::sqrt( 0.5 );
    const matrixT right = identityMatrix( 2 );
    vectorT singular( 2 );
    singular << 7, 2;
    const std::vector<Eigen::Index> deleted{ 0 };

    REQUIRE( mx::math::validateSvdDeletionFactor( left ) == statusT::success );
    matrixT invalidFactor = left;
    invalidFactor( 1, 1 ) *= 2;
    REQUIRE( mx::math::validateSvdDeletionFactor( invalidFactor ) == statusT::factorNotOrthonormal );

    mx::math::svdDeletionResult<realT> result;
    mx::math::svdDeletionWorkspace<realT> workspace;
    for( const backendT backend : { backendT::leadingCovariance, backendT::stableCore } )
    {
        REQUIRE( mx::math::svdDeletionSucceeded(
            mx::math::svdRemoveRows( result, singular, left, deleted, 2, workspace, backend ) ) );
        REQUIRE( result.squaredSingularValues()( 0 ) == Approx( 4.0 ) );
        REQUIRE( result.squaredSingularValues()( 1 ) == Approx( 0.0 ).margin( 1e-12 ) );
    }

    const std::vector<Eigen::Index> noDeletion;
    REQUIRE( mx::math::svdRemoveColumns( result, singular, right, noDeletion, 2, workspace, backendT::stableCore ) ==
             statusT::success );
    REQUIRE( ( result.rotation() - identityMatrix( 2 ) ).matrix().norm() == Approx( 0.0 ) );
    REQUIRE( result.minimumPSDValue() == Approx( 1.0 ) );

    matrixT noDeletedRows( 0, 2 );
    REQUIRE( mx::math::svdDeletionLeadingCore( result, singular, noDeletedRows, 2, workspace ) == statusT::success );
    REQUIRE( result.minimumPSDValue() == Approx( 4.0 / 49.0 ) );
    REQUIRE( mx::math::svdDeletionStableCore( result, singular, noDeletedRows, 2, workspace ) == statusT::success );
    REQUIRE( result.minimumPSDValue() == Approx( 1.0 ) );

    vectorT zeroSingular( 2 );
    zeroSingular << 5, 0;
    matrixT deletedRow( 1, 2 );
    deletedRow << std::sqrt( 0.5 ), std::sqrt( 0.5 );
    REQUIRE( mx::math::svdDeletionSucceeded(
        mx::math::svdDeletionCore( result, zeroSingular, deletedRow, 2, workspace, backendT::stableCore ) ) );
    REQUIRE( result.singularValues()( 1 ) == Approx( 0.0 ).margin( 1e-13 ) );

    vectorT allZero( 2 );
    allZero.setZero();
    for( const backendT backend : { backendT::leadingCovariance, backendT::stableCore } )
    {
        REQUIRE( mx::math::svdDeletionSucceeded(
            mx::math::svdDeletionCore( result, allZero, deletedRow, 2, workspace, backend ) ) );
        REQUIRE( result.singularValues().matrix().norm() == Approx( 0.0 ) );
        REQUIRE( ( result.rotation() - identityMatrix( 2 ) ).matrix().norm() == Approx( 0.0 ) );
    }

    matrixT zeroLeverageFactor( 3, 2 );
    zeroLeverageFactor << 1, 0, 0, 1, 0, 0;
    vectorT zeroLeverageSingular( 2 );
    zeroLeverageSingular << 2, 1;
    const matrixT zeroLeverageMatrix =
        representedMatrix( zeroLeverageFactor, zeroLeverageSingular, identityMatrix( 2 ) );
    const std::vector<Eigen::Index> zeroLeverageDeletion{ 2 };
    for( const backendT backend : { backendT::leadingCovariance, backendT::stableCore } )
    {
        INFO( "backend: " << mx::math::svdDeletionBackendName( backend ) );
        mx::math::svdDeletionResult<realT> zeroLeverageResult;
        mx::math::svdDeletionWorkspace<realT> zeroLeverageWorkspace;
        REQUIRE( mx::math::svdDeletionSucceeded( mx::math::svdRemoveRows( zeroLeverageResult,
                                                                          zeroLeverageSingular,
                                                                          zeroLeverageFactor,
                                                                          zeroLeverageDeletion,
                                                                          2,
                                                                          zeroLeverageWorkspace,
                                                                          backend ) ) );
        const deletionComparison comparison =
            compareRowResult( zeroLeverageMatrix, identityMatrix( 2 ), zeroLeverageDeletion, zeroLeverageResult );
        REQUIRE( comparison.squaredSingularError <= comparison.squaredSingularTolerance );
        REQUIRE( comparison.singularError <= comparison.singularTolerance );
        REQUIRE( comparison.covarianceError <= comparison.covarianceTolerance );
    }
}

/// SVD deletion supports float and workspace reuse
/** Verifies float specialization accuracy and reusable result/workspace capacity through svdDeletionCore.
 *
 * \ingroup svdDowndate_unit_tests
 */
TEST_CASE( "SVD deletion supports float and workspace reuse", "[math::svdDowndate][float][workspace]" )
{
    STATIC_REQUIRE( !std::is_copy_constructible_v<mx::math::svdDeletionResult<float>> );
    STATIC_REQUIRE( !std::is_copy_assignable_v<mx::math::svdDeletionResult<float>> );
    STATIC_REQUIRE( std::is_nothrow_move_constructible_v<mx::math::svdDeletionResult<float>> );
    STATIC_REQUIRE( std::is_nothrow_move_assignable_v<mx::math::svdDeletionResult<float>> );

    using floatMatrixT = mx::math::svdDeletionMatrix<float>;
    using floatVectorT = mx::math::svdDeletionVector<float>;
    floatVectorT singular( 3 );
    singular << 6, 3, 1;
    floatMatrixT deleted( 1, 3 );
    deleted << 0.5f, 0.25f, 0.125f;

    mx::math::svdDeletionResult<float> result;
    mx::math::svdDeletionWorkspace<float> workspace;
    REQUIRE( workspace.prepare( 3, 4, backendT::leadingCovariance ) == statusT::success );
    REQUIRE( mx::math::svdDeletionCore( result, singular, deleted, 3, workspace, backendT::leadingCovariance ) ==
             statusT::success );
    REQUIRE( workspace.baseRank() == 3 );
    REQUIRE( workspace.maximumDeleted() == 4 );
    REQUIRE( workspace.backend() == backendT::leadingCovariance );
    REQUIRE( workspace.prepared() );
    REQUIRE( result.baseRank() == 3 );
    REQUIRE( result.outputRank() == 3 );
    REQUIRE( result.lapackInfo() == 0 );
    REQUIRE( result.minimumPSDValue() >= -1e-4f );

    const floatVectorT leadingSquared = result.squaredSingularValues();
    REQUIRE( workspace.prepare( 3, 4, backendT::stableCore ) == statusT::success );
    REQUIRE( mx::math::svdDeletionSucceeded(
        mx::math::svdDeletionCore( result, singular, deleted, 3, workspace, backendT::stableCore ) ) );
    REQUIRE( ( result.squaredSingularValues() - leadingSquared ).matrix().norm() < 2e-4f );
    floatMatrixT floatIdentity( 3, 3 );
    floatIdentity.matrix().setIdentity();
    REQUIRE( mx::math::validateSvdDeletionFactor( floatIdentity ) == statusT::success );

    mx::math::svdDeletionWorkspace<float> movedWorkspace( std::move( workspace ) );
    REQUIRE( !workspace.prepared() );
    REQUIRE( movedWorkspace.prepared() );
    REQUIRE( movedWorkspace.maximumDeleted() == 4 );
    mx::math::svdDeletionWorkspace<float> assignedWorkspace;
    assignedWorkspace = std::move( movedWorkspace );
    REQUIRE( !movedWorkspace.prepared() );
    REQUIRE( movedWorkspace.baseRank() == 0 );
    REQUIRE( movedWorkspace.maximumDeleted() == 0 );
    REQUIRE( movedWorkspace.backend() == backendT::stableCore );
    REQUIRE( movedWorkspace.lapackInfo() == 0 );
    movedWorkspace.clear();
    REQUIRE( !movedWorkspace.prepared() );
    REQUIRE( assignedWorkspace.prepared() );
    REQUIRE( mx::math::svdDeletionSucceeded(
        mx::math::svdDeletionCore( result, singular, deleted, 3, assignedWorkspace, backendT::stableCore ) ) );
    REQUIRE( movedWorkspace.prepare( 3, 1, backendT::leadingCovariance ) == statusT::success );
    REQUIRE( movedWorkspace.prepared() );
    movedWorkspace.clear();

    const floatMatrixT physicalLeft = sineFactor( 6, 3 ).cast<float>();
    const floatMatrixT physicalRight = sineFactor( 5, 3 ).cast<float>();
    floatVectorT physicalSingular( 3 );
    physicalSingular << 7, 3, 0.5f;
    const Eigen::MatrixXf physicalMatrix =
        physicalLeft.matrix() * physicalSingular.matrix().asDiagonal() * physicalRight.matrix().transpose();
    Eigen::MatrixXf retainedMatrix( 4, 5 );
    retainedMatrix.row( 0 ) = physicalMatrix.row( 1 );
    retainedMatrix.row( 1 ) = physicalMatrix.row( 2 );
    retainedMatrix.row( 2 ) = physicalMatrix.row( 3 );
    retainedMatrix.row( 3 ) = physicalMatrix.row( 5 );
    Eigen::JacobiSVD<Eigen::MatrixXf> directFloat( retainedMatrix, Eigen::ComputeThinV );
    const std::vector<Eigen::Index> physicalDeleted{ 0, 4 };
    for( const backendT backend : { backendT::leadingCovariance, backendT::stableCore } )
    {
        mx::math::svdDeletionResult<float> physicalResult;
        mx::math::svdDeletionWorkspace<float> physicalWorkspace;
        REQUIRE( mx::math::svdDeletionSucceeded( mx::math::svdRemoveRows( physicalResult,
                                                                          physicalSingular,
                                                                          physicalLeft,
                                                                          physicalDeleted,
                                                                          3,
                                                                          physicalWorkspace,
                                                                          backend ) ) );
        REQUIRE( ( physicalResult.singularValues().matrix() - directFloat.singularValues().head( 3 ) ).norm() < 2e-4f );
    }

    mx::math::svdDeletionResult<float> movedResult( std::move( result ) );
    REQUIRE( movedResult.baseRank() == 3 );
    mx::math::svdDeletionResult<float> assignedResult;
    assignedResult = std::move( movedResult );
    REQUIRE( assignedResult.baseRank() == 3 );
    REQUIRE( movedResult.baseRank() == 0 );
    REQUIRE( movedResult.outputRank() == 0 );
    REQUIRE( movedResult.maximumOutputRank() == 0 );
    REQUIRE( movedResult.status() == statusT::notComputed );
    REQUIRE( movedResult.backend() == backendT::stableCore );
    REQUIRE( movedResult.clampedEigenvalues() == 0 );
    REQUIRE( movedResult.minimumPSDValue() == 0 );
    REQUIRE( movedResult.lapackInfo() == 0 );
    REQUIRE( movedResult.singularValues().size() == 0 );
    REQUIRE( movedResult.squaredSingularValues().size() == 0 );
    REQUIRE( movedResult.rotation().size() == 0 );
    REQUIRE( mx::math::svdDeletionCore( movedResult, singular, deleted, 0, assignedWorkspace, backendT::stableCore ) ==
             statusT::invalidInput );
    REQUIRE( movedResult.status() == statusT::invalidInput );
    REQUIRE( movedResult.prepare( 3, 2 ) == statusT::success );
    REQUIRE( movedResult.baseRank() == 3 );
    REQUIRE( movedResult.outputRank() == 2 );

    assignedWorkspace.clear();
    REQUIRE( !assignedWorkspace.prepared() );
    REQUIRE( assignedWorkspace.baseRank() == 0 );
    REQUIRE( assignedWorkspace.maximumDeleted() == 0 );
    REQUIRE( assignedWorkspace.lapackInfo() == 0 );
}

/// SVD deletion accepts views and reuses prepared capacity
/** Verifies non-owning Eigen views and capacity reuse without reallocation through svdDeletionCore and
 * validateSvdDeletionFactor.
 *
 * \ingroup svdDowndate_unit_tests
 */
TEST_CASE( "SVD deletion accepts views and reuses prepared capacity", "[math::svdDowndate][workspace][views]" )
{
    hookGuard guard;
    vectorT singularStorage( 4 );
    singularStorage << 8, 4, 2, 0.5;
    matrixT factorStorage( 6, 4 );
    factorStorage.setZero();
    factorStorage.matrix().leftCols( 3 ) = sineFactor( 6, 3 ).matrix();
    matrixT deletedStorage( 3, 4 );
    deletedStorage.setZero();
    deletedStorage.matrix().topLeftCorner( 1, 3 ) = factorStorage.matrix().block( 2, 0, 1, 3 );

    REQUIRE( mx::math::validateSvdDeletionFactor( factorStorage.leftCols( 3 ) ) == statusT::success );
    mx::math::svdDeletionResult<realT> result;
    mx::math::svdDeletionWorkspace<realT> workspace;
    REQUIRE( result.prepare( 3, 3 ) == statusT::success );
    REQUIRE( workspace.prepare( 3, 3, backendT::stableCore ) == statusT::success );

    failingOperation = mx::math::detail::svdDeletionTestOperation::prepareWorkspace;
    mx::math::detail::svdDeletionHooks<double>().operation = throwAllocation;
    REQUIRE( mx::math::svdDeletionSucceeded( mx::math::svdDeletionCore( result,
                                                                        singularStorage.head( 3 ),
                                                                        deletedStorage.topLeftCorner( 1, 3 ),
                                                                        2,
                                                                        workspace,
                                                                        backendT::stableCore ) ) );
    REQUIRE( workspace.maximumDeleted() == 3 );
    REQUIRE( result.maximumOutputRank() == 3 );
    REQUIRE( result.outputRank() == 2 );
    REQUIRE( result.rotation().cols() == 2 );
    REQUIRE( result.singularValues().size() == 3 );

    matrixT tooManyDeleted( 4, 3 );
    tooManyDeleted.setZero();
    REQUIRE( mx::math::svdDeletionCore( result,
                                        singularStorage.head( 3 ),
                                        tooManyDeleted,
                                        3,
                                        workspace,
                                        backendT::stableCore ) == statusT::allocationFailure );
}

/// SVD deletion accepts under-aligned consumer storage
/** Verifies that validateSvdDeletionFactor and svdRemoveRows use unaligned raw views at the shared-library boundary
 * when the consumer's Eigen packet-alignment setting is smaller than mxlib's setting.
 *
 * \ingroup svdDowndate_unit_tests
 */
TEST_CASE( "SVD deletion accepts under-aligned consumer storage", "[math::svdDowndate][abi][alignment]" )
{
#ifdef MXLIB_SVD_DELETION_TEST_CONSUMER_ALIGNMENT
    STATIC_REQUIRE( EIGEN_MAX_ALIGN_BYTES == 16 );
    STATIC_REQUIRE( EIGEN_MAX_STATIC_ALIGN_BYTES == 16 );
#endif

    constexpr Eigen::Index rows{ 7 };
    constexpr Eigen::Index rank{ 3 };
    constexpr std::uintptr_t testedAlignment{ 32 };

    std::vector<realT> factorStorage( static_cast<std::size_t>( rows * rank ) + 4 );
    realT *factorData = factorStorage.data();
    while( reinterpret_cast<std::uintptr_t>( factorData ) % testedAlignment != 16 )
    {
        ++factorData;
    }
    REQUIRE( reinterpret_cast<std::uintptr_t>( factorData ) % testedAlignment == 16 );

    using unalignedMatrixMap = Eigen::Map<matrixT, Eigen::Unaligned>;
    unalignedMatrixMap factor( factorData, rows, rank );
    factor = sineFactor( rows, rank );
    REQUIRE( mx::math::validateSvdDeletionFactor( factor ) == statusT::success );

    std::vector<realT> singularStorage( static_cast<std::size_t>( rank ) + 4 );
    realT *singularData = singularStorage.data();
    while( reinterpret_cast<std::uintptr_t>( singularData ) % testedAlignment != 16 )
    {
        ++singularData;
    }
    REQUIRE( reinterpret_cast<std::uintptr_t>( singularData ) % testedAlignment == 16 );

    using unalignedVectorMap = Eigen::Map<vectorT, Eigen::Unaligned>;
    unalignedVectorMap singularValues( singularData, rank );
    singularValues << 5, 3, 1;

    const std::vector<Eigen::Index> deletedRows{ 2 };
    mx::math::svdDeletionResult<realT> result;
    mx::math::svdDeletionWorkspace<realT> workspace;
    REQUIRE(
        mx::math::svdRemoveRows( result, singularValues, factor, deletedRows, rank, workspace, backendT::stableCore ) ==
        statusT::success );
    REQUIRE( result.rotation().rows() == rank );
    REQUIRE( result.rotation().cols() == rank );
    REQUIRE( result.rotation().matrix().squaredNorm() == Approx( static_cast<realT>( rank ) ).margin( 1e-12 ) );
    REQUIRE( result.singularValues().matrix().squaredNorm() > 0 );
}

/// SVD deletion accepts a default empty ABI index descriptor
/** Verifies that the raw svdRemoveRowsAbiV2 entry point treats a default-constructed empty index descriptor as an
 * identity deletion.
 *
 * \ingroup svdDowndate_unit_tests
 */
TEST_CASE( "SVD deletion accepts a default empty ABI index descriptor", "[math::svdDowndate][abi][identity]" )
{
    vectorT singular( 2 );
    singular << 4, 1;
    matrixT factor = identityMatrix( 2 );
    mx::math::svdDeletionResult<realT> result;
    mx::math::svdDeletionWorkspace<realT> workspace;

    const mx::math::svdDeletionConstVectorViewV2<realT> singularView{ singular.data(), singular.size() };
    const mx::math::svdDeletionConstMatrixViewV2<realT> factorView{ factor.data(),
                                                                    factor.rows(),
                                                                    factor.cols(),
                                                                    factor.outerStride() };
    const mx::math::svdDeletionConstIndexViewV2 noDeletedIndices;
    REQUIRE( mx::math::detail::svdRemoveRowsAbiV2<realT>( result,
                                                          singularView,
                                                          factorView,
                                                          noDeletedIndices,
                                                          2,
                                                          workspace,
                                                          backendT::stableCore ) == statusT::success );
    REQUIRE( ( result.singularValues() - singular ).matrix().norm() == Approx( 0.0 ) );
    REQUIRE( ( result.rotation() - identityMatrix( 2 ) ).matrix().norm() == Approx( 0.0 ) );
}

/// SVD deletion rejects malformed ABI descriptors
/** Verifies the ABI-v2 raw entry points reject vector, matrix, and index descriptors whose claimed storage cannot be
 * addressed safely.
 *
 * \ingroup svdDowndate_unit_tests
 */
TEST_CASE( "SVD deletion rejects malformed ABI descriptors", "[math::svdDowndate][abi][errors]" )
{
    vectorT singular( 2 );
    singular << 4, 1;
    matrixT factor = identityMatrix( 2 );
    const Eigen::Index deletedIndex{ 0 };
    mx::math::svdDeletionResult<realT> result;
    mx::math::svdDeletionWorkspace<realT> workspace;

    const mx::math::svdDeletionConstVectorViewV2<realT> singularView{ singular.data(), singular.size() };
    const mx::math::svdDeletionConstMatrixViewV2<realT> factorView{ factor.data(),
                                                                    factor.rows(),
                                                                    factor.cols(),
                                                                    factor.outerStride() };
    const auto *misalignedScalar =
        reinterpret_cast<const realT *>( reinterpret_cast<const unsigned char *>( singular.data() ) + 1 );

    REQUIRE( result.prepare( -1, 1 ) == statusT::invalidInput );
    REQUIRE( mx::math::detail::validateSvdDeletionFactorAbiV2( { factor.data(), -1, 2, factor.outerStride() }, 0 ) ==
             statusT::invalidInput );
    REQUIRE( mx::math::detail::validateSvdDeletionFactorAbiV2(
                 mx::math::svdDeletionConstMatrixViewV2<realT>{ nullptr, 2, 2, 2 },
                 0.0 ) == statusT::invalidInput );
    REQUIRE( mx::math::detail::validateSvdDeletionFactorAbiV2( { misalignedScalar, 2, 2, 2 }, 0.0 ) ==
             statusT::invalidInput );
    REQUIRE( mx::math::detail::validateSvdDeletionFactorAbiV2( { factor.data(), 2, 2, 1 }, 0 ) ==
             statusT::invalidInput );
    REQUIRE(
        mx::math::detail::validateSvdDeletionFactorAbiV2(
            { factor.data(), std::numeric_limits<std::int64_t>::max(), 1, std::numeric_limits<std::int64_t>::max() },
            0 ) == statusT::invalidInput );
    REQUIRE( mx::math::detail::validateSvdDeletionFactorAbiV2(
                 { factor.data(), 2, 2, std::numeric_limits<std::int64_t>::max() },
                 0 ) == statusT::invalidInput );
    REQUIRE( mx::math::detail::validateSvdDeletionFactorAbiV2(
                 mx::math::svdDeletionConstMatrixViewV2<float>{ nullptr, 1, 1, 1 },
                 0.0F ) == statusT::invalidInput );
    REQUIRE( mx::math::detail::svdDeletionLeadingCoreAbiV2<realT>(
                 result,
                 { singular.data(), std::numeric_limits<std::int64_t>::max() },
                 factorView,
                 2,
                 workspace ) == statusT::invalidInput );
    REQUIRE( mx::math::detail::svdDeletionLeadingCoreAbiV2<realT>( result,
                                                                   { misalignedScalar, singular.size() },
                                                                   factorView,
                                                                   2,
                                                                   workspace ) == statusT::invalidInput );
    REQUIRE( mx::math::detail::svdDeletionStableCoreAbiV2<realT>(
                 result,
                 { singular.data(), std::numeric_limits<std::int64_t>::max() },
                 factorView,
                 2,
                 workspace ) == statusT::invalidInput );
    REQUIRE( mx::math::detail::svdRemoveRowsAbiV2<realT>(
                 result,
                 singularView,
                 factorView,
                 { &deletedIndex, std::numeric_limits<std::int64_t>::max(), sizeof( deletedIndex ) },
                 2,
                 workspace,
                 backendT::stableCore ) == statusT::invalidInput );
    REQUIRE( mx::math::detail::svdRemoveRowsAbiV2<realT>( result,
                                                          singularView,
                                                          factorView,
                                                          { &deletedIndex, 1, 3 },
                                                          2,
                                                          workspace,
                                                          backendT::stableCore ) == statusT::invalidInput );
    REQUIRE( mx::math::detail::svdRemoveRowsAbiV2<realT>( result,
                                                          singularView,
                                                          factorView,
                                                          { nullptr, 1, sizeof( deletedIndex ) },
                                                          2,
                                                          workspace,
                                                          backendT::stableCore ) == statusT::invalidInput );
    REQUIRE( mx::math::detail::svdRemoveColumnsAbiV2<realT>( result,
                                                             singularView,
                                                             factorView,
                                                             { nullptr, 1, sizeof( deletedIndex ) },
                                                             2,
                                                             workspace,
                                                             backendT::stableCore ) == statusT::invalidInput );

    const std::int8_t deleted8{ 0 };
    const std::int16_t deleted16{ 0 };
    const std::int32_t deleted32{ 0 };
    for( const mx::math::svdDeletionConstIndexViewV2 indices :
         { mx::math::svdDeletionConstIndexViewV2{ &deleted8, 1, sizeof( deleted8 ) },
           mx::math::svdDeletionConstIndexViewV2{ &deleted16, 1, sizeof( deleted16 ) },
           mx::math::svdDeletionConstIndexViewV2{ &deleted32, 1, sizeof( deleted32 ) } } )
    {
        REQUIRE( mx::math::svdDeletionSucceeded(
            mx::math::detail::svdRemoveRowsAbiV2<
                realT>( result, singularView, factorView, indices, 2, workspace, backendT::leadingCovariance ) ) );
    }
}

/// SVD deletion reports invalid inputs and allocation failures
/** Verifies public validation, status text, allocation failures, and complete-side rejection in the SVD
 * deletion system.
 *
 * \ingroup svdDowndate_unit_tests
 */
TEST_CASE( "SVD deletion reports invalid inputs and allocation failures", "[math::svdDowndate][errors]" )
{
    hookGuard guard;
    REQUIRE( std::string( mx::math::svdDeletionBackendName( backendT::leadingCovariance ) ) == "leadingCovariance" );
    REQUIRE( std::string( mx::math::svdDeletionBackendName( backendT::stableCore ) ) == "stableCore" );
    REQUIRE( std::string( mx::math::svdDeletionBackendName( static_cast<backendT>( 99 ) ) ) == "unknown" );

    for( const statusT status : { statusT::notComputed,
                                  statusT::success,
                                  statusT::successWithClamping,
                                  statusT::invalidInput,
                                  statusT::allocationFailure,
                                  statusT::workspaceQueryFailure,
                                  statusT::solverFailure,
                                  statusT::nonFiniteOutput,
                                  statusT::invalidSolverOutput,
                                  statusT::rescalingOverflow,
                                  statusT::nonPositiveSemidefinite,
                                  statusT::factorNotOrthonormal } )
    {
        REQUIRE( std::string( mx::math::svdDeletionStatusName( status ) ) != "unknown" );
    }
    REQUIRE( std::string( mx::math::svdDeletionStatusName( static_cast<statusT>( 99 ) ) ) == "unknown" );
    REQUIRE( mx::math::svdDeletionSucceeded( statusT::success ) );
    REQUIRE( mx::math::svdDeletionSucceeded( statusT::successWithClamping ) );
    REQUIRE( !mx::math::svdDeletionSucceeded( statusT::solverFailure ) );

    vectorT singular( 2 );
    singular << 4, 1;
    matrixT factor = identityMatrix( 2 );
    mx::math::svdDeletionResult<realT> result;
    mx::math::svdDeletionWorkspace<realT> workspace;

    mx::math::svdDeletionResult<realT> movedFromResult;
    mx::math::svdDeletionResult<realT> resultOwner( std::move( movedFromResult ) );
    failingOperation = mx::math::detail::svdDeletionTestOperation::prepareResult;
    mx::math::detail::svdDeletionHooks<double>().operation = throwAllocation;
    REQUIRE( mx::math::svdDeletionCore( movedFromResult, singular, factor, 0, workspace, backendT::stableCore ) ==
             statusT::allocationFailure );
    mx::math::detail::svdDeletionHooks<double>().operation = throwLengthError;
    REQUIRE( movedFromResult.prepare( 2, 2 ) == statusT::allocationFailure );
    mx::math::detail::svdDeletionHooks<double>().operation = nullptr;
    REQUIRE( movedFromResult.prepare( 2, 2 ) == statusT::success );

    mx::math::svdDeletionWorkspace<realT> movedFromWorkspace;
    mx::math::svdDeletionWorkspace<realT> workspaceOwner( std::move( movedFromWorkspace ) );
    failingOperation = mx::math::detail::svdDeletionTestOperation::prepareWorkspace;
    mx::math::detail::svdDeletionHooks<double>().operation = throwAllocation;
    REQUIRE( movedFromWorkspace.prepare( 2, 1, backendT::stableCore ) == statusT::allocationFailure );
    mx::math::detail::svdDeletionHooks<double>().operation = throwLengthError;
    REQUIRE( movedFromWorkspace.prepare( 2, 1, backendT::stableCore ) == statusT::allocationFailure );
    mx::math::detail::svdDeletionHooks<double>().operation = nullptr;
    REQUIRE( movedFromWorkspace.prepare( 2, 1, backendT::stableCore ) == statusT::success );

    REQUIRE( result.status() == statusT::notComputed );
    REQUIRE( result.prepare( 0, 0 ) == statusT::invalidInput );
    REQUIRE( result.prepare( 2, 3 ) == statusT::invalidInput );
    REQUIRE( workspace.prepare( 0, -1, backendT::stableCore ) == statusT::invalidInput );
    REQUIRE( workspace.prepare( 2, 1, static_cast<backendT>( 99 ) ) == statusT::invalidInput );

    failingOperation = mx::math::detail::svdDeletionTestOperation::prepareResult;
    mx::math::detail::svdDeletionHooks<double>().operation = throwAllocation;
    REQUIRE( result.prepare( 2, 2 ) == statusT::allocationFailure );

    failingOperation = mx::math::detail::svdDeletionTestOperation::prepareWorkspace;
    REQUIRE( workspace.prepare( 2, 1, backendT::stableCore ) == statusT::allocationFailure );

    failingOperation = mx::math::detail::svdDeletionTestOperation::validateFactor;
    REQUIRE( mx::math::validateSvdDeletionFactor( factor ) == statusT::allocationFailure );
    mx::math::detail::svdDeletionHooks<double>().operation = nullptr;

    failingOperation = mx::math::detail::svdDeletionTestOperation::prepareResult;
    mx::math::detail::svdDeletionHooks<double>().operation = throwLengthError;
    REQUIRE( result.prepare( 3, 3 ) == statusT::allocationFailure );
    failingOperation = mx::math::detail::svdDeletionTestOperation::prepareWorkspace;
    REQUIRE( workspace.prepare( 3, 1, backendT::stableCore ) == statusT::allocationFailure );
    failingOperation = mx::math::detail::svdDeletionTestOperation::validateFactor;
    REQUIRE( mx::math::validateSvdDeletionFactor( factor ) == statusT::allocationFailure );
    mx::math::detail::svdDeletionHooks<double>().operation = nullptr;

    matrixT noDeletedRows( 0, 2 );
    failingOperation = mx::math::detail::svdDeletionTestOperation::prepareResult;
    mx::math::detail::svdDeletionHooks<double>().operation = throwAllocation;
    mx::math::svdDeletionResult<realT> identityFailure;
    REQUIRE( mx::math::svdDeletionLeadingCore( identityFailure, singular, noDeletedRows, 2, workspace ) ==
             statusT::allocationFailure );

    matrixT oneDeletedRow( 1, 2 );
    oneDeletedRow << 0.25, 0.5;
    mx::math::svdDeletionResult<realT> coreResultFailure;
    REQUIRE( mx::math::svdDeletionLeadingCore( coreResultFailure, singular, oneDeletedRow, 2, workspace ) ==
             statusT::allocationFailure );

    mx::math::svdDeletionResult<realT> removeResultFailure;
    const std::vector<Eigen::Index> oneRow{ 0 };
    REQUIRE( mx::math::svdRemoveRows( removeResultFailure,
                                      singular,
                                      factor,
                                      oneRow,
                                      2,
                                      workspace,
                                      backendT::leadingCovariance ) == statusT::allocationFailure );
    mx::math::detail::svdDeletionHooks<double>().operation = nullptr;

    REQUIRE( workspace.prepare( 2, 1, backendT::stableCore ) == statusT::success );
    REQUIRE( workspace.prepare( 2, 1, backendT::stableCore ) == statusT::success );

    const std::vector<Eigen::Index> allRows{ 0, 1 };
    REQUIRE( mx::math::svdRemoveRows( result, singular, factor, allRows, 2, workspace, backendT::stableCore ) ==
             statusT::invalidInput );
    const std::vector<Eigen::Index> duplicate{ 0, 0 };
    REQUIRE( mx::math::svdRemoveRows( result, singular, factor, duplicate, 2, workspace, backendT::stableCore ) ==
             statusT::invalidInput );
    const std::vector<Eigen::Index> unsorted{ 1, 0 };
    REQUIRE( mx::math::svdRemoveRows( result, singular, factor, unsorted, 2, workspace, backendT::stableCore ) ==
             statusT::invalidInput );
    const std::vector<Eigen::Index> negative{ -1 };
    REQUIRE( mx::math::svdRemoveRows( result, singular, factor, negative, 2, workspace, backendT::stableCore ) ==
             statusT::invalidInput );
    const std::vector<Eigen::Index> outOfRange{ 2 };
    REQUIRE( mx::math::svdRemoveColumns( result, singular, factor, outOfRange, 2, workspace, backendT::stableCore ) ==
             statusT::invalidInput );

    matrixT nonfinite = factor;
    nonfinite( 0, 0 ) = std::numeric_limits<realT>::infinity();
    REQUIRE( mx::math::validateSvdDeletionFactor( nonfinite ) == statusT::invalidInput );
    REQUIRE( mx::math::validateSvdDeletionFactor( factor, realT( -1 ) ) == statusT::invalidInput );

    matrixT tooWide( 1, 2 );
    tooWide.setZero();
    REQUIRE( mx::math::validateSvdDeletionFactor( tooWide ) == statusT::invalidInput );
    const std::vector<Eigen::Index> noIndices;
    REQUIRE( mx::math::svdRemoveRows( result, singular, tooWide, noIndices, 2, workspace, backendT::stableCore ) ==
             statusT::invalidInput );
    matrixT emptyFactor( 2, 0 );
    REQUIRE( mx::math::validateSvdDeletionFactor( emptyFactor ) == statusT::invalidInput );
    matrixT overflowFactor( 2, 1 );
    overflowFactor.setConstant( std::numeric_limits<realT>::max() );
    REQUIRE( mx::math::validateSvdDeletionFactor( overflowFactor ) == statusT::invalidInput );

    vectorT invalidSingular = singular;
    invalidSingular( 1 ) = -1;
    const std::vector<Eigen::Index> one{ 0 };
    REQUIRE( mx::math::svdRemoveRows( result, invalidSingular, factor, one, 2, workspace, backendT::stableCore ) ==
             statusT::invalidInput );
    invalidSingular << 1, 2;
    REQUIRE( mx::math::svdRemoveRows( result, invalidSingular, factor, one, 2, workspace, backendT::stableCore ) ==
             statusT::invalidInput );
    invalidSingular = singular;
    invalidSingular( 0 ) = std::numeric_limits<realT>::infinity();
    REQUIRE( mx::math::svdRemoveRows( result, invalidSingular, factor, one, 2, workspace, backendT::stableCore ) ==
             statusT::invalidInput );
    REQUIRE( mx::math::svdRemoveRows( result, singular, factor, one, 0, workspace, backendT::stableCore ) ==
             statusT::invalidInput );

    matrixT wrongColumns( 0, 1 );
    REQUIRE( mx::math::svdDeletionCore( result, singular, wrongColumns, 2, workspace, backendT::stableCore ) ==
             statusT::invalidInput );
    REQUIRE( mx::math::svdDeletionLeadingCore( result, singular, wrongColumns, 2, workspace ) ==
             statusT::invalidInput );
    matrixT nonfiniteDeleted( 1, 2 );
    nonfiniteDeleted << 0, std::numeric_limits<realT>::infinity();
    REQUIRE( mx::math::svdDeletionCore( result, singular, nonfiniteDeleted, 2, workspace, backendT::stableCore ) ==
             statusT::invalidInput );
    REQUIRE(
        mx::math::svdDeletionCore( result, singular, factor.topRows( 1 ), 2, workspace, static_cast<backendT>( 99 ) ) ==
        statusT::invalidInput );

    matrixT selectedFinite = factor;
    selectedFinite( 1, 0 ) = std::numeric_limits<realT>::infinity();
    REQUIRE( mx::math::svdDeletionSucceeded(
        mx::math::svdRemoveRows( result, singular, selectedFinite, one, 2, workspace, backendT::leadingCovariance ) ) );
    const std::vector<Eigen::Index> selectNonfinite{ 1 };
    REQUIRE( mx::math::svdRemoveRows( result,
                                      singular,
                                      selectedFinite,
                                      selectNonfinite,
                                      2,
                                      workspace,
                                      backendT::leadingCovariance ) == statusT::invalidInput );

    const Eigen::Index huge = static_cast<Eigen::Index>( std::numeric_limits<MXLAPACK_INT>::max() );
    REQUIRE( workspace.prepare( huge, 0, backendT::stableCore ) == statusT::invalidInput );
    REQUIRE( workspace.prepare( huge, 1, backendT::stableCore ) == statusT::invalidInput );
}

/// SVD deletion reports workspace query failures
/** Verifies SYEVR and GESVD workspace-query failure reporting through svdDeletionWorkspace::prepare and
 * svdDeletionCore.
 *
 * \ingroup svdDowndate_unit_tests
 */
TEST_CASE( "SVD deletion reports workspace query failures", "[math::svdDowndate][errors][query]" )
{
    hookGuard guard;

    SECTION( "SYEVR query INFO" )
    {
        mx::math::svdDeletionWorkspace<realT> workspace;
        syevrMode = solverHookMode::queryFailure;
        mx::math::detail::svdDeletionHooks<double>().syevr = syevrHook;
        REQUIRE( workspace.prepare( 3, 1, backendT::leadingCovariance ) == statusT::workspaceQueryFailure );
        REQUIRE( workspace.lapackInfo() == 61 );
        REQUIRE( !workspace.prepared() );
    }

    SECTION( "SYEVR floating query size" )
    {
        mx::math::svdDeletionWorkspace<realT> workspace;
        syevrMode = solverHookMode::invalidQuery;
        mx::math::detail::svdDeletionHooks<double>().syevr = syevrHook;
        REQUIRE( workspace.prepare( 3, 1, backendT::leadingCovariance ) == statusT::workspaceQueryFailure );
        REQUIRE( workspace.lapackInfo() == 0 );
    }

    SECTION( "SYEVR integer query size" )
    {
        mx::math::svdDeletionWorkspace<realT> workspace;
        syevrMode = solverHookMode::invalidIntegerQuery;
        mx::math::detail::svdDeletionHooks<double>().syevr = syevrHook;
        REQUIRE( workspace.prepare( 3, 1, backendT::stableCore ) == statusT::workspaceQueryFailure );
    }

    SECTION( "GESVD query INFO" )
    {
        mx::math::svdDeletionWorkspace<realT> workspace;
        gesvdMode = solverHookMode::queryFailure;
        mx::math::detail::svdDeletionHooks<double>().gesvd = gesvdHook;
        REQUIRE( workspace.prepare( 3, 1, backendT::stableCore ) == statusT::workspaceQueryFailure );
        REQUIRE( workspace.lapackInfo() == 62 );
    }

    SECTION( "GESVD query size" )
    {
        mx::math::svdDeletionWorkspace<realT> workspace;
        gesvdMode = solverHookMode::invalidQuery;
        mx::math::detail::svdDeletionHooks<double>().gesvd = gesvdHook;
        REQUIRE( workspace.prepare( 3, 1, backendT::stableCore ) == statusT::workspaceQueryFailure );
        REQUIRE( workspace.lapackInfo() == 0 );
    }

    SECTION( "query failure propagated to a result" )
    {
        vectorT singular( 2 );
        singular << 4, 1;
        matrixT deleted( 1, 2 );
        deleted << 0.25, 0.5;
        mx::math::svdDeletionResult<realT> result;
        mx::math::svdDeletionWorkspace<realT> workspace;
        syevrMode = solverHookMode::queryFailure;
        mx::math::detail::svdDeletionHooks<double>().syevr = syevrHook;
        REQUIRE( mx::math::svdDeletionCore( result, singular, deleted, 2, workspace, backendT::leadingCovariance ) ==
                 statusT::workspaceQueryFailure );
        REQUIRE( result.status() == statusT::workspaceQueryFailure );
        REQUIRE( result.lapackInfo() == 61 );
        REQUIRE( result.minimumPSDValue() == Approx( 0.0 ) );
    }
}

/// SVD deletion reports numerical solver outcomes
/** Verifies solve failures, malformed solver output, PSD clamping, indefiniteness, and stable rescaling through
 * svdDeletionLeadingCore and svdDeletionStableCore.
 *
 * \ingroup svdDowndate_unit_tests
 */
TEST_CASE( "SVD deletion reports numerical solver outcomes", "[math::svdDowndate][errors][solver]" )
{
    hookGuard guard;
    vectorT singular( 3 );
    singular << 8, 4, 2;
    matrixT oneDeleted( 1, 3 );
    oneDeleted << 0.25, 0.2, 0.1;
    matrixT twoDeleted( 2, 3 );
    twoDeleted << 0.25, 0.2, 0.1, 0.1, 0.15, 0.2;

    SECTION( "leading SYEVR solve failure" )
    {
        mx::math::svdDeletionResult<realT> result;
        mx::math::svdDeletionWorkspace<realT> workspace;
        REQUIRE( workspace.prepare( 3, 1, backendT::leadingCovariance ) == statusT::success );
        syevrMode = solverHookMode::solveFailure;
        mx::math::detail::svdDeletionHooks<double>().syevr = syevrHook;
        REQUIRE( mx::math::svdDeletionLeadingCore( result, singular, oneDeleted, 3, workspace ) ==
                 statusT::solverFailure );
        REQUIRE( result.lapackInfo() == 71 );
    }

    SECTION( "leading SYEVR count mismatch" )
    {
        mx::math::svdDeletionResult<realT> result;
        mx::math::svdDeletionWorkspace<realT> workspace;
        REQUIRE( workspace.prepare( 3, 1, backendT::leadingCovariance ) == statusT::success );
        syevrMode = solverHookMode::countMismatch;
        mx::math::detail::svdDeletionHooks<double>().syevr = syevrHook;
        REQUIRE( mx::math::svdDeletionLeadingCore( result, singular, oneDeleted, 3, workspace ) ==
                 statusT::solverFailure );
        REQUIRE( result.lapackInfo() == 0 );
    }

    SECTION( "leading non-finite values and vectors" )
    {
        for( const solverHookMode mode : { solverHookMode::nonFiniteValue, solverHookMode::nonFiniteVector } )
        {
            mx::math::svdDeletionResult<realT> result;
            mx::math::svdDeletionWorkspace<realT> workspace;
            REQUIRE( workspace.prepare( 3, 1, backendT::leadingCovariance ) == statusT::success );
            syevrMode = mode;
            mx::math::detail::svdDeletionHooks<double>().syevr = syevrHook;
            REQUIRE( mx::math::svdDeletionLeadingCore( result, singular, oneDeleted, 3, workspace ) ==
                     statusT::nonFiniteOutput );
        }
    }

    SECTION( "leading invalid ordering" )
    {
        mx::math::svdDeletionResult<realT> result;
        mx::math::svdDeletionWorkspace<realT> workspace;
        REQUIRE( workspace.prepare( 3, 1, backendT::leadingCovariance ) == statusT::success );
        syevrMode = solverHookMode::invalidOrdering;
        mx::math::detail::svdDeletionHooks<double>().syevr = syevrHook;
        REQUIRE( mx::math::svdDeletionLeadingCore( result, singular, oneDeleted, 3, workspace ) ==
                 statusT::invalidSolverOutput );
    }

    SECTION( "leading clamping and indefiniteness" )
    {
        mx::math::svdDeletionResult<realT> result;
        mx::math::svdDeletionWorkspace<realT> workspace;
        REQUIRE( workspace.prepare( 3, 1, backendT::leadingCovariance ) == statusT::success );
        mx::math::detail::svdDeletionHooks<double>().syevr = syevrHook;
        syevrMode = solverHookMode::roundoffClamp;
        REQUIRE( mx::math::svdDeletionLeadingCore( result, singular, oneDeleted, 3, workspace ) ==
                 statusT::successWithClamping );
        REQUIRE( result.clampedEigenvalues() == 1 );
        REQUIRE( result.minimumPSDValue() == Approx( -std::numeric_limits<realT>::epsilon() ) );

        syevrMode = solverHookMode::indefinite;
        REQUIRE( mx::math::svdDeletionLeadingCore( result, singular, oneDeleted, 3, workspace ) ==
                 statusT::nonPositiveSemidefinite );
        REQUIRE( result.minimumPSDValue() == Approx( -0.25 ) );
    }

    SECTION( "stable complement solver outcomes" )
    {
        mx::math::svdDeletionResult<realT> result;
        mx::math::svdDeletionWorkspace<realT> workspace;
        REQUIRE( workspace.prepare( 3, 2, backendT::stableCore ) == statusT::success );
        mx::math::detail::svdDeletionHooks<double>().syevr = syevrHook;

        for( const solverHookMode mode : { solverHookMode::nonFiniteValue, solverHookMode::nonFiniteVector } )
        {
            syevrMode = mode;
            REQUIRE( mx::math::svdDeletionStableCore( result, singular, oneDeleted, 3, workspace ) ==
                     statusT::nonFiniteOutput );
        }
        syevrMode = solverHookMode::countMismatch;
        REQUIRE( mx::math::svdDeletionStableCore( result, singular, oneDeleted, 3, workspace ) ==
                 statusT::solverFailure );
        syevrMode = solverHookMode::invalidOrdering;
        REQUIRE( mx::math::svdDeletionStableCore( result, singular, twoDeleted, 3, workspace ) ==
                 statusT::invalidSolverOutput );
        syevrMode = solverHookMode::indefinite;
        REQUIRE( mx::math::svdDeletionStableCore( result, singular, oneDeleted, 3, workspace ) ==
                 statusT::nonPositiveSemidefinite );
        REQUIRE( result.minimumPSDValue() == Approx( -0.25 ) );
        syevrMode = solverHookMode::roundoffClamp;
        REQUIRE( mx::math::svdDeletionStableCore( result, singular, oneDeleted, 3, workspace ) ==
                 statusT::successWithClamping );
        REQUIRE( result.clampedEigenvalues() == 1 );
    }

    SECTION( "stable GESVD solver outcomes" )
    {
        for( const auto [mode, expected] : std::vector<std::pair<solverHookMode, statusT>>{
                 { solverHookMode::solveFailure, statusT::solverFailure },
                 { solverHookMode::nonFiniteValue, statusT::nonFiniteOutput },
                 { solverHookMode::nonFiniteVector, statusT::nonFiniteOutput },
                 { solverHookMode::invalidOrdering, statusT::invalidSolverOutput },
                 { solverHookMode::negativeSpectrum, statusT::invalidSolverOutput } } )
        {
            mx::math::svdDeletionResult<realT> result;
            mx::math::svdDeletionWorkspace<realT> workspace;
            REQUIRE( workspace.prepare( 3, 1, backendT::stableCore ) == statusT::success );
            gesvdMode = mode;
            mx::math::detail::svdDeletionHooks<double>().gesvd = gesvdHook;
            REQUIRE( mx::math::svdDeletionStableCore( result, singular, oneDeleted, 3, workspace ) == expected );
            if( mode == solverHookMode::solveFailure )
            {
                REQUIRE( result.lapackInfo() == 72 );
            }
            mx::math::detail::svdDeletionHooks<double>().gesvd = nullptr;
        }
    }

    SECTION( "stable rescaling retains a tiny normalized singular value" )
    {
        vectorT largeSingular( 1 );
        largeSingular << 1e200;
        matrixT deleted( 1, 1 );
        deleted << 0.5;
        mx::math::svdDeletionResult<realT> result;
        mx::math::svdDeletionWorkspace<realT> workspace;
        REQUIRE( workspace.prepare( 1, 1, backendT::stableCore ) == statusT::success );
        gesvdMode = solverHookMode::tinySpectrum;
        mx::math::detail::svdDeletionHooks<double>().gesvd = gesvdHook;
        REQUIRE( mx::math::svdDeletionStableCore( result, largeSingular, deleted, 1, workspace ) == statusT::success );
        REQUIRE( result.singularValues()( 0 ) == Approx( 1.0 ).epsilon( 1e-12 ) );
        REQUIRE( result.squaredSingularValues()( 0 ) == Approx( 1.0 ).epsilon( 1e-12 ) );
    }
}

/// SVD deletion reports rescaling overflow
/** Verifies explicit squared-result overflow reporting in unchanged, leading-core, and stable-core deletion
 * paths.
 *
 * \ingroup svdDowndate_unit_tests
 */
TEST_CASE( "SVD deletion reports rescaling overflow", "[math::svdDowndate][errors][overflow]" )
{
    vectorT singular( 1 );
    singular << std::numeric_limits<realT>::max() / 2;
    matrixT factor( 2, 1 );
    factor << std::sqrt( 0.5 ), std::sqrt( 0.5 );
    const std::vector<Eigen::Index> none;
    const std::vector<Eigen::Index> deletedIndex{ 0 };

    for( const backendT backend : { backendT::leadingCovariance, backendT::stableCore } )
    {
        CAPTURE( mx::math::svdDeletionBackendName( backend ) );
        mx::math::svdDeletionResult<realT> result;
        mx::math::svdDeletionWorkspace<realT> workspace;
        REQUIRE( mx::math::svdRemoveRows( result, singular, factor, none, 1, workspace, backend ) ==
                 statusT::rescalingOverflow );
        const statusT deletionStatus =
            mx::math::svdRemoveRows( result, singular, factor, deletedIndex, 1, workspace, backend );
        CAPTURE( result.singularValues()( 0 ), result.squaredSingularValues()( 0 ) );
        REQUIRE( deletionStatus == statusT::rescalingOverflow );
    }
}

} // namespace unitTest::math_svdDowndate_test
