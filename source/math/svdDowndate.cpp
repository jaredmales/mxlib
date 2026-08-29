/** \file
 * \brief Definitions of reusable thin-SVD row and column deletion updates.
 *
 * \ingroup gen_math_files
 */

//***********************************************************************//
// Copyright 2026 Jared R. Males (jaredmales@gmail.com)
//
// This file is part of mxlib.
//
// mxlib is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// mxlib is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with mxlib.  If not, see <http://www.gnu.org/licenses/>.
//***********************************************************************//

#include "math/svdDowndate.hpp"

#include "math/floatUtils.hpp"

#include <algorithm>
#include <cmath>
#include <limits>
#include <new>
#include <stdexcept>
#include <utility>
#include <vector>

namespace mx
{
namespace math
{

template <typename realT>
class svdDeletionResult<realT>::storage
{
  public:
    /// Complete descending updated singular spectrum.
    svdDeletionVector<realT> m_singularValues;

    /// Complete descending squared singular spectrum.
    svdDeletionVector<realT> m_squaredSingularValues;

    /// Preserved-side rotation storage; only the leading output-rank columns are active.
    svdDeletionMatrix<realT> m_rotation;

    /// Most recent operation status.
    svdDeletionStatus m_status{ svdDeletionStatus::notComputed };

    /// Backend that produced the current outputs.
    svdDeletionBackend m_backend{ svdDeletionBackend::stableCore };

    /// Prepared base rank.
    Eigen::Index m_baseRank{ 0 };

    /// Prepared output rank.
    Eigen::Index m_outputRank{ 0 };

    /// Allocated maximum output rank for the current base rank.
    Eigen::Index m_maximumOutputRank{ 0 };

    /// Number of clamped roundoff-scale negative eigenvalues.
    Eigen::Index m_clampedEigenvalues{ 0 };

    /// Smallest pre-clamp eigenvalue from the normalized backend PSD core.
    realT m_minimumPSDValue{ 0 };

    /// Underlying LAPACK failure code, or zero.
    MXLAPACK_INT m_lapackInfo{ 0 };
};

template <typename realT>
class svdDeletionWorkspace<realT>::storage
{
  public:
    /// Gathered deleted rows, with the active operation occupying the leading rows.
    svdDeletionMatrix<realT> m_deletedRows;

    /// Base singular values normalized by their maximum to prevent overflow.
    svdDeletionVector<realT> m_normalizedSingularValues;

    /// Scaled deleted-row transpose used by the symmetric core.
    svdDeletionMatrix<realT> m_scaledDeletedTranspose;

    /// Symmetric covariance or complement matrix overwritten by SYEVR.
    svdDeletionMatrix<realT> m_psdCore;

    /// Full eigenvectors returned by SYEVR.
    svdDeletionMatrix<realT> m_eigenvectors;

    /// Full ascending eigenvalues returned by SYEVR.
    svdDeletionVector<realT> m_eigenvalues;

    /// Complement square root used by the stable backend.
    svdDeletionMatrix<realT> m_complementRoot;

    /// Complement-preserving rectangular core overwritten by GESVD.
    svdDeletionMatrix<realT> m_stableCore;

    /// Right singular vectors returned row-wise by GESVD.
    svdDeletionMatrix<realT> m_rightTranspose;

    /// Singular values returned by GESVD.
    svdDeletionVector<realT> m_stableSingularValues;

    /// Eigenvector support indices required by SYEVR.
    std::vector<MXLAPACK_INT> m_support;

    /// Floating-point workspace required by SYEVR.
    std::vector<realT> m_syevrWork;

    /// Integer workspace required by SYEVR.
    std::vector<MXLAPACK_INT> m_syevrIWork;

    /// Floating-point workspace required by GESVD.
    std::vector<realT> m_gesvdWork;

    /// Prepared base rank.
    Eigen::Index m_baseRank{ 0 };

    /// Prepared maximum deletion count.
    Eigen::Index m_maximumDeleted{ 0 };

    /// Prepared numerical backend.
    svdDeletionBackend m_backend{ svdDeletionBackend::stableCore };

    /// True after successful dimensioning and LAPACK workspace queries.
    bool m_prepared{ false };

    /// Underlying LAPACK workspace-query failure code, or zero.
    MXLAPACK_INT m_lapackInfo{ 0 };
};

namespace detail
{

namespace
{

// Return true when an Eigen dimension can be passed to the configured LAPACK integer ABI.
bool validLapackDimension( Eigen::Index dimension /* [in] candidate dimension */ )
{
    return dimension >= 0 && static_cast<unsigned long long>( dimension ) <=
                                 static_cast<unsigned long long>( std::numeric_limits<MXLAPACK_INT>::max() );
}

// Return true when a dense array shape fits Eigen indexing and addressable storage.
template <typename valueT>
bool validArrayShape( Eigen::Index rows, /* [in] candidate row count */
                      Eigen::Index cols /* [in] candidate column count */ )
{
    if( rows < 0 || cols < 0 || ( rows > 0 && cols > std::numeric_limits<Eigen::Index>::max() / rows ) )
    {
        return false; // LCOV_EXCL_LINE -- requires dimensions wider than the configured LAPACK integer ABI.
    }
    const Eigen::Index elements = rows * cols;
    return static_cast<unsigned long long>( elements ) <=
           static_cast<unsigned long long>( std::numeric_limits<std::size_t>::max() / sizeof( valueT ) );
}

// Return true when an element count fits addressable `std::vector` storage.
template <typename valueT>
bool validVectorSize( Eigen::Index elements /* [in] candidate element count */ )
{
    return elements >= 0 &&
           static_cast<unsigned long long>( elements ) <=
               static_cast<unsigned long long>( std::numeric_limits<std::size_t>::max() / sizeof( valueT ) );
}

// Return true when every dense Eigen element is finite under mxlib's fast-math-safe predicate.
template <typename derivedT>
bool allFinite( const Eigen::DenseBase<derivedT> &values /* [in] dense object to inspect */ )
{
    for( Eigen::Index index = 0; index < values.size(); ++index )
    {
        if( !math::isFinite( values.derived().data()[index] ) )
        {
            return false;
        }
    }
    return true;
}

// Return true when the leading values are finite, nonnegative, and nonincreasing.
template <typename derivedT>
bool validDescendingSpectrum( const Eigen::DenseBase<derivedT> &values /* [in] spectrum to inspect */ )
{
    for( Eigen::Index index = 0; index < values.size(); ++index )
    {
        const auto value = values.derived()( index );
        if( !math::isFinite( value ) || value < 0 || ( index > 0 && value > values.derived()( index - 1 ) ) )
        {
            return false;
        }
    }
    return true;
}

// Return true when the leading values are finite and nondecreasing.
template <typename derivedT>
bool validAscendingSpectrum( const Eigen::DenseBase<derivedT> &values /* [in] spectrum to inspect */ )
{
    for( Eigen::Index index = 0; index < values.size(); ++index )
    {
        const auto value = values.derived()( index );
        if( !math::isFinite( value ) || ( index > 0 && value < values.derived()( index - 1 ) ) )
        {
            return false;
        }
    }
    return true;
}

// Invoke an optional deterministic operation hook.
template <typename realT>
void callOperationHook( svdDeletionTestOperation operation /* [in] operation about to execute */ )
{
    auto &hooks = svdDeletionHooks<realT>();
    if( hooks.operation )
    {
        hooks.operation( operation );
    }
}

// Invoke the injected or production SYEVR implementation.
template <typename realT>
MXLAPACK_INT callSyevr( char jobz,              /* [in] eigenvector request */
                        char range,             /* [in] eigenvalue selection mode */
                        char uplo,              /* [in] populated input triangle */
                        MXLAPACK_INT n,         /* [in] matrix order */
                        realT *matrix,          /* [in,out] symmetric input matrix */
                        MXLAPACK_INT lda,       /* [in] input leading dimension */
                        realT valueLower,       /* [in] lower value bound */
                        realT valueUpper,       /* [in] upper value bound */
                        MXLAPACK_INT indexLow,  /* [in] one-based lower index */
                        MXLAPACK_INT indexHigh, /* [in] one-based upper index */
                        realT tolerance,        /* [in] convergence tolerance */
                        MXLAPACK_INT *found,    /* [out] selected eigenvalue count */
                        realT *eigenvalues,     /* [out] ascending eigenvalues */
                        realT *eigenvectors,    /* [out] eigenvectors in columns */
                        MXLAPACK_INT ldz,       /* [in] eigenvector leading dimension */
                        MXLAPACK_INT *support,  /* [out] eigenvector support */
                        realT *work,            /* [in,out] floating workspace */
                        MXLAPACK_INT lwork,     /* [in] floating workspace size */
                        MXLAPACK_INT *iwork,    /* [in,out] integer workspace */
                        MXLAPACK_INT liwork /* [in] integer workspace size */ )
{
    auto &hooks = svdDeletionHooks<realT>();
    if( hooks.syevr )
    {
        return hooks.syevr( jobz,
                            range,
                            uplo,
                            n,
                            matrix,
                            lda,
                            valueLower,
                            valueUpper,
                            indexLow,
                            indexHigh,
                            tolerance,
                            found,
                            eigenvalues,
                            eigenvectors,
                            ldz,
                            support,
                            work,
                            lwork,
                            iwork,
                            liwork );
    }

    return math::syevr<realT>( jobz,
                               range,
                               uplo,
                               n,
                               matrix,
                               lda,
                               valueLower,
                               valueUpper,
                               indexLow,
                               indexHigh,
                               tolerance,
                               found,
                               eigenvalues,
                               eigenvectors,
                               ldz,
                               support,
                               work,
                               lwork,
                               iwork,
                               liwork );
}

// Invoke the injected or production GESVD implementation.
template <typename realT>
MXLAPACK_INT callGesvd( char jobu,         /* [in] left singular-vector request */
                        char jobvt,        /* [in] right singular-vector request */
                        MXLAPACK_INT rows, /* [in] matrix row count */
                        MXLAPACK_INT cols, /* [in] matrix column count */
                        realT *matrix,     /* [in,out] input matrix */
                        MXLAPACK_INT lda,  /* [in] input leading dimension */
                        realT *singular,   /* [out] descending singular values */
                        realT *left,       /* [out] left singular vectors when requested */
                        MXLAPACK_INT ldu,  /* [in] left-vector leading dimension */
                        realT *rightT,     /* [out] transposed right singular vectors */
                        MXLAPACK_INT ldvt, /* [in] right-vector leading dimension */
                        realT *work,       /* [in,out] floating workspace */
                        MXLAPACK_INT lwork /* [in] floating workspace size */ )
{
    auto &hooks = svdDeletionHooks<realT>();
    if( hooks.gesvd )
    {
        return hooks.gesvd( jobu, jobvt, rows, cols, matrix, lda, singular, left, ldu, rightT, ldvt, work, lwork );
    }

    return math::gesvd<realT>( jobu, jobvt, rows, cols, matrix, lda, singular, left, ldu, rightT, ldvt, work, lwork );
}

// Convert a finite LAPACK floating workspace query to an integer size.
template <typename realT>
bool querySize( MXLAPACK_INT &size, /* [out] validated integer workspace size */
                realT query /* [in] floating-point LAPACK query result */ )
{
    if( !math::isFinite( query ) || query < realT( 1 ) ||
        static_cast<long double>( query ) > static_cast<long double>( std::numeric_limits<MXLAPACK_INT>::max() ) )
    {
        return false;
    }

    size = static_cast<MXLAPACK_INT>( std::ceil( query ) );
    return size >= 1;
}

// Return the construction-scale tolerance for a theoretically PSD core.
template <typename realT>
realT psdTolerance( Eigen::Index dimension, /* [in] core dimension */
                    realT scale /* [in] nonnegative core construction scale */ )
{
    return realT( 64 ) * std::numeric_limits<realT>::epsilon() *
           static_cast<realT>( std::max<Eigen::Index>( 1, dimension ) ) * scale;
}

} // namespace

template <typename realT>
svdDeletionTestHooks<realT> &svdDeletionHooks()
{
    static svdDeletionTestHooks<realT> hooks;
    return hooks;
}

// Implementation access shared by the public deletion entry points.
template <typename realT>
struct svdDeletionImplementation
{
    using resultT = svdDeletionResult<realT>;
    using workspaceT = svdDeletionWorkspace<realT>;

    // Publish a failure and its optional underlying LAPACK status.
    static svdDeletionStatus fail( resultT &result,          /* [out] result status to invalidate */
                                   svdDeletionStatus status, /* [in] failure status */
                                   MXLAPACK_INT lapackInfo = 0 /* [in] underlying LAPACK status */ )
    {
        if( !result.ensureStorage() )
        {
            return svdDeletionStatus::allocationFailure;
        }
        result.m_storage->m_status = status;
        result.m_storage->m_lapackInfo = lapackInfo;
        result.m_storage->m_clampedEigenvalues = 0;
        if( status != svdDeletionStatus::nonPositiveSemidefinite )
        {
            result.m_storage->m_minimumPSDValue = realT( 0 );
        }
        return status;
    }

    // Return true for a supported numerical backend value.
    static bool validBackend( svdDeletionBackend backend /* [in] backend value to validate */ )
    {
        return backend == svdDeletionBackend::leadingCovariance || backend == svdDeletionBackend::stableCore;
    }

    // Validate singular values and requested rank without checking factor rows.
    static bool validSingularInputs( svdDeletionConstVectorRef<realT> singularValues, /* [in] base singular values */
                                     Eigen::Index outputRank /* [in] requested output rank */ )
    {
        const Eigen::Index rank = singularValues.size();
        if( rank <= 0 || !validLapackDimension( rank ) || outputRank <= 0 || outputRank > rank ||
            !allFinite( singularValues ) )
        {
            return false;
        }

        for( Eigen::Index index = 0; index < rank; ++index )
        {
            if( singularValues( index ) < realT( 0 ) ||
                ( index > 0 && singularValues( index ) > singularValues( index - 1 ) ) )
            {
                return false;
            }
        }
        return true;
    }

    // Validate common core inputs without checking global factor orthonormality.
    static bool validCoreInputs( svdDeletionConstVectorRef<realT> singularValues, /* [in] base singular values */
                                 svdDeletionConstMatrixRef<realT> deletedRows,    /* [in] deleted factor rows */
                                 Eigen::Index outputRank /* [in] requested output rank */ )
    {
        const Eigen::Index rank = singularValues.size();
        return validSingularInputs( singularValues, outputRank ) && deletedRows.cols() == rank &&
               validLapackDimension( deletedRows.rows() ) &&
               deletedRows.rows() <= std::numeric_limits<MXLAPACK_INT>::max() - rank && allFinite( deletedRows );
    }

    // Normalize base singular values and return their common scale.
    static realT normalizeSingularValues( workspaceT &workspace, /* [out] normalized values */
                                          svdDeletionConstVectorRef<realT> singularValues /* [in] base values */ )
    {
        const realT scale = singularValues( 0 );
        if( scale == realT( 0 ) )
        {
            workspace.m_storage->m_normalizedSingularValues.setZero();
            return scale;
        }

        int scaleExponent{ 0 };
        const realT scaleMantissa = std::frexp( scale, &scaleExponent );
        for( Eigen::Index index = 0; index < singularValues.size(); ++index )
        {
            int valueExponent{ 0 };
            const realT valueMantissa = std::frexp( singularValues( index ), &valueExponent );
            workspace.m_storage->m_normalizedSingularValues( index ) =
                std::ldexp( valueMantissa / scaleMantissa, valueExponent - scaleExponent );
        }
        return scale;
    }

    // Return a ratio without forming an underflowed reciprocal under fast-math builds.
    static realT normalizedRatio( realT value, /* [in] nonnegative numerator */
                                  realT scale /* [in] positive denominator */ )
    {
        int scaleExponent{ 0 };
        int valueExponent{ 0 };
        const realT scaleMantissa = std::frexp( scale, &scaleExponent );
        const realT valueMantissa = std::frexp( value, &valueExponent );
        return std::ldexp( valueMantissa / scaleMantissa, valueExponent - scaleExponent );
    }

    // Rescale one normalized singular value without squaring before restoring input scale.
    static bool rescaleSingularValue( realT &singularValue,        /* [out] rescaled singular value */
                                      realT &squaredSingularValue, /* [out] rescaled squared value */
                                      realT normalizedValue,       /* [in] normalized singular value */
                                      realT scale /* [in] base singular-value scale */ )
    {
        const realT maximum = std::numeric_limits<realT>::max();
        singularValue = normalizedValue * scale;
        if( !math::isFinite( singularValue ) || singularValue > std::sqrt( maximum ) )
        {
            return false;
        }
        squaredSingularValue = singularValue * singularValue;
        return math::isFinite( squaredSingularValue );
    }

    // Rescale one normalized squared singular value without overflowing either public representation.
    static bool rescaleValue( realT &singularValue,         /* [out] rescaled singular value */
                              realT &squaredSingularValue,  /* [out] rescaled squared singular value */
                              realT normalizedSquaredValue, /* [in] normalized squared value */
                              realT scale /* [in] base singular-value scale */ )
    {
        return rescaleSingularValue( singularValue, squaredSingularValue, std::sqrt( normalizedSquaredValue ), scale );
    }

    // Ensure result and worker storage can execute one core operation.
    static svdDeletionStatus prepare( resultT &result,           /* [in,out] output storage */
                                      workspaceT &workspace,     /* [in,out] reusable scratch */
                                      Eigen::Index rank,         /* [in] base rank */
                                      Eigen::Index deletedCount, /* [in] actual deleted-row count */
                                      Eigen::Index outputRank,   /* [in] requested output rank */
                                      svdDeletionBackend backend /* [in] selected backend */ )
    {
        svdDeletionStatus status = result.prepare( rank, outputRank );
        if( status != svdDeletionStatus::success )
        {
            return status;
        }

        if( !workspace.prepared() || workspace.baseRank() != rank || workspace.maximumDeleted() < deletedCount ||
            workspace.backend() != backend )
        {
            status = workspace.prepare( rank, deletedCount, backend );
            if( status != svdDeletionStatus::success )
            {
                return fail( result, status, workspace.lapackInfo() );
            }
        }

        return svdDeletionStatus::success;
    }

    // Copy deleted factor rows into the active portion of a prepared workspace.
    static void loadDeletedRows( workspaceT &workspace, /* [out] gathered workspace rows */
                                 svdDeletionConstMatrixRef<realT> deletedRows /* [in] rows to copy */ )
    {
        if( deletedRows.rows() > 0 )
        {
            workspace.m_storage->m_deletedRows.matrix().topRows( deletedRows.rows() ) = deletedRows.matrix();
        }
    }

    // Publish the unchanged represented singular system for an empty deletion.
    static svdDeletionStatus identity( resultT &result,                                 /* [out] identity update */
                                       svdDeletionConstVectorRef<realT> singularValues, /* [in] base singular values */
                                       Eigen::Index outputRank,                         /* [in] requested output rank */
                                       svdDeletionBackend backend /* [in] reported backend */ )
    {
        const svdDeletionStatus status = result.prepare( singularValues.size(), outputRank );
        if( status != svdDeletionStatus::success )
        {
            return status;
        }

        result.m_storage->m_rotation.setZero();
        for( Eigen::Index index = 0; index < singularValues.size(); ++index )
        {
            result.m_storage->m_singularValues( index ) = singularValues( index );
            if( singularValues( index ) > realT( 0 ) &&
                singularValues( index ) > std::sqrt( std::numeric_limits<realT>::max() ) )
            {
                return fail( result, svdDeletionStatus::rescalingOverflow );
            }
            result.m_storage->m_squaredSingularValues( index ) = singularValues( index ) * singularValues( index );
            if( index < outputRank )
            {
                result.m_storage->m_rotation( index, index ) = realT( 1 );
            }
        }
        result.m_storage->m_backend = backend;
        result.m_storage->m_status = svdDeletionStatus::success;
        result.m_storage->m_clampedEigenvalues = 0;
        if( backend == svdDeletionBackend::leadingCovariance )
        {
            const realT scale = singularValues( 0 );
            realT minimumNormalized{ 0 };
            if( scale != realT( 0 ) )
            {
                minimumNormalized = normalizedRatio( singularValues( singularValues.size() - 1 ), scale );
            }
            result.m_storage->m_minimumPSDValue = minimumNormalized * minimumNormalized;
        }
        else
        {
            result.m_storage->m_minimumPSDValue = realT( 1 );
        }
        result.m_storage->m_lapackInfo = 0;
        return result.m_storage->m_status;
    }

    // Publish an arbitrary orthonormal system for an identically zero represented matrix.
    static svdDeletionStatus zeroSystem( resultT &result,         /* [out] zero-valued result */
                                         Eigen::Index outputRank, /* [in] requested output rank */
                                         svdDeletionBackend backend /* [in] reported backend */ )
    {
        result.m_storage->m_rotation.setZero();
        result.m_storage->m_singularValues.setZero();
        result.m_storage->m_squaredSingularValues.setZero();
        for( Eigen::Index index = 0; index < outputRank; ++index )
        {
            result.m_storage->m_rotation( index, index ) = realT( 1 );
        }
        result.m_storage->m_backend = backend;
        result.m_storage->m_status = svdDeletionStatus::success;
        result.m_storage->m_clampedEigenvalues = 0;
        result.m_storage->m_minimumPSDValue = realT( 0 );
        result.m_storage->m_lapackInfo = 0;
        return result.m_storage->m_status;
    }

    // Execute a normalized symmetric covariance solve with deleted rows already loaded.
    static svdDeletionStatus leadingLoaded( resultT &result,                                 /* [out] updated system */
                                            svdDeletionConstVectorRef<realT> singularValues, /* [in] base values */
                                            Eigen::Index deletedCount,                       /* [in] active rows */
                                            Eigen::Index outputRank,                         /* [in] requested rank */
                                            workspaceT &workspace /* [in,out] scratch */ )
    {
        const Eigen::Index rank = singularValues.size();
        const realT scale = normalizeSingularValues( workspace, singularValues );
        if( scale == realT( 0 ) )
        {
            return zeroSystem( result, outputRank, svdDeletionBackend::leadingCovariance );
        }

        for( Eigen::Index column = 0; column < deletedCount; ++column )
        {
            for( Eigen::Index row = 0; row < rank; ++row )
            {
                workspace.m_storage->m_scaledDeletedTranspose( row, column ) =
                    workspace.m_storage->m_normalizedSingularValues( row ) *
                    workspace.m_storage->m_deletedRows( column, row );
            }
        }

        const auto scaledDeleted = workspace.m_storage->m_scaledDeletedTranspose.matrix().leftCols( deletedCount );
        workspace.m_storage->m_psdCore.matrix().noalias() = -scaledDeleted * scaledDeleted.transpose();
        for( Eigen::Index index = 0; index < rank; ++index )
        {
            workspace.m_storage->m_psdCore( index, index ) += workspace.m_storage->m_normalizedSingularValues( index ) *
                                                              workspace.m_storage->m_normalizedSingularValues( index );
        }
        MXLAPACK_INT found{ 0 };
        const MXLAPACK_INT lapackRank = static_cast<MXLAPACK_INT>( rank );
        const MXLAPACK_INT info =
            callSyevr<realT>( 'V',
                              'A',
                              'L',
                              lapackRank,
                              workspace.m_storage->m_psdCore.data(),
                              lapackRank,
                              realT( 0 ),
                              realT( 0 ),
                              1,
                              lapackRank,
                              math::lamch<realT>( 'S' ),
                              &found,
                              workspace.m_storage->m_eigenvalues.data(),
                              workspace.m_storage->m_eigenvectors.data(),
                              lapackRank,
                              workspace.m_storage->m_support.data(),
                              workspace.m_storage->m_syevrWork.data(),
                              static_cast<MXLAPACK_INT>( workspace.m_storage->m_syevrWork.size() ),
                              workspace.m_storage->m_syevrIWork.data(),
                              static_cast<MXLAPACK_INT>( workspace.m_storage->m_syevrIWork.size() ) );
        if( info != 0 || found != lapackRank )
        {
            return fail( result, svdDeletionStatus::solverFailure, info );
        }
        if( !allFinite( workspace.m_storage->m_eigenvalues ) || !allFinite( workspace.m_storage->m_eigenvectors ) )
        {
            return fail( result, svdDeletionStatus::nonFiniteOutput );
        }
        if( !validAscendingSpectrum( workspace.m_storage->m_eigenvalues ) )
        {
            return fail( result, svdDeletionStatus::invalidSolverOutput );
        }

        result.m_storage->m_minimumPSDValue = workspace.m_storage->m_eigenvalues( 0 );
        const realT tolerance =
            psdTolerance<realT>( rank, realT( 1 ) + static_cast<realT>( scaledDeleted.squaredNorm() ) );
        if( result.m_storage->m_minimumPSDValue < -tolerance )
        {
            return fail( result, svdDeletionStatus::nonPositiveSemidefinite );
        }

        result.m_storage->m_clampedEigenvalues = 0;
        for( Eigen::Index index = 0; index < rank; ++index )
        {
            if( workspace.m_storage->m_eigenvalues( index ) < realT( 0 ) )
            {
                workspace.m_storage->m_eigenvalues( index ) = realT( 0 );
                ++result.m_storage->m_clampedEigenvalues;
            }
        }

        for( Eigen::Index output = 0; output < rank; ++output )
        {
            const Eigen::Index source = rank - 1 - output;
            if( !rescaleValue( result.m_storage->m_singularValues( output ),
                               result.m_storage->m_squaredSingularValues( output ),
                               workspace.m_storage->m_eigenvalues( source ),
                               scale ) )
            {
                return fail( result, svdDeletionStatus::rescalingOverflow );
            }
            if( output < outputRank )
            {
                result.m_storage->m_rotation.matrix().col( output ) =
                    workspace.m_storage->m_eigenvectors.matrix().col( source );
            }
        }

        result.m_storage->m_backend = svdDeletionBackend::leadingCovariance;
        result.m_storage->m_lapackInfo = 0;
        result.m_storage->m_status = result.m_storage->m_clampedEigenvalues > 0 ? svdDeletionStatus::successWithClamping
                                                                                : svdDeletionStatus::success;
        return result.m_storage->m_status;
    }

    // Execute the symmetric covariance deletion core after workspace loading.
    static svdDeletionStatus leading( resultT &result,                                 /* [out] updated system */
                                      svdDeletionConstVectorRef<realT> singularValues, /* [in] base values */
                                      svdDeletionConstMatrixRef<realT> deletedRows,    /* [in] deleted rows */
                                      Eigen::Index outputRank,                         /* [in] requested rank */
                                      workspaceT &workspace /* [in,out] scratch */ )
    {
        if( !validCoreInputs( singularValues, deletedRows, outputRank ) )
        {
            return fail( result, svdDeletionStatus::invalidInput );
        }
        if( deletedRows.rows() == 0 )
        {
            return identity( result, singularValues, outputRank, svdDeletionBackend::leadingCovariance );
        }

        const Eigen::Index rank = singularValues.size();
        svdDeletionStatus status =
            prepare( result, workspace, rank, deletedRows.rows(), outputRank, svdDeletionBackend::leadingCovariance );
        if( status != svdDeletionStatus::success )
        {
            return status;
        }
        loadDeletedRows( workspace, deletedRows );
        return leadingLoaded( result, singularValues, deletedRows.rows(), outputRank, workspace );
    }

    // Execute a normalized complement-preserving SVD with deleted rows already loaded.
    static svdDeletionStatus stableLoaded( resultT &result,                                 /* [out] updated system */
                                           svdDeletionConstVectorRef<realT> singularValues, /* [in] base values */
                                           Eigen::Index deletedCount,                       /* [in] active rows */
                                           Eigen::Index outputRank,                         /* [in] requested rank */
                                           workspaceT &workspace /* [in,out] scratch */ )
    {
        const Eigen::Index rank = singularValues.size();
        const realT scale = normalizeSingularValues( workspace, singularValues );
        if( scale == realT( 0 ) )
        {
            return zeroSystem( result, outputRank, svdDeletionBackend::stableCore );
        }

        const auto factorRows = workspace.m_storage->m_deletedRows.matrix().topRows( deletedCount );
        auto complement = workspace.m_storage->m_psdCore.matrix().topLeftCorner( deletedCount, deletedCount );
        complement.setIdentity();
        complement.noalias() -= factorRows * factorRows.transpose();
        MXLAPACK_INT found{ 0 };
        const MXLAPACK_INT lapackComplement = static_cast<MXLAPACK_INT>( deletedCount );
        MXLAPACK_INT info = callSyevr<realT>( 'V',
                                              'A',
                                              'L',
                                              lapackComplement,
                                              workspace.m_storage->m_psdCore.data(),
                                              static_cast<MXLAPACK_INT>( workspace.m_storage->m_psdCore.rows() ),
                                              realT( 0 ),
                                              realT( 0 ),
                                              1,
                                              lapackComplement,
                                              math::lamch<realT>( 'S' ),
                                              &found,
                                              workspace.m_storage->m_eigenvalues.data(),
                                              workspace.m_storage->m_eigenvectors.data(),
                                              static_cast<MXLAPACK_INT>( workspace.m_storage->m_eigenvectors.rows() ),
                                              workspace.m_storage->m_support.data(),
                                              workspace.m_storage->m_syevrWork.data(),
                                              static_cast<MXLAPACK_INT>( workspace.m_storage->m_syevrWork.size() ),
                                              workspace.m_storage->m_syevrIWork.data(),
                                              static_cast<MXLAPACK_INT>( workspace.m_storage->m_syevrIWork.size() ) );
        if( info != 0 || found != lapackComplement )
        {
            return fail( result, svdDeletionStatus::solverFailure, info );
        }
        for( Eigen::Index index = 0; index < deletedCount; ++index )
        {
            if( !math::isFinite( workspace.m_storage->m_eigenvalues( index ) ) )
            {
                return fail( result, svdDeletionStatus::nonFiniteOutput );
            }
            for( Eigen::Index row = 0; row < deletedCount; ++row )
            {
                if( !math::isFinite( workspace.m_storage->m_eigenvectors( row, index ) ) )
                {
                    return fail( result, svdDeletionStatus::nonFiniteOutput );
                }
            }
        }
        if( !validAscendingSpectrum( workspace.m_storage->m_eigenvalues.head( deletedCount ) ) )
        {
            return fail( result, svdDeletionStatus::invalidSolverOutput );
        }

        result.m_storage->m_minimumPSDValue = workspace.m_storage->m_eigenvalues( 0 );
        const realT tolerance =
            psdTolerance<realT>( deletedCount, realT( 1 ) + static_cast<realT>( factorRows.squaredNorm() ) );
        if( result.m_storage->m_minimumPSDValue < -tolerance )
        {
            return fail( result, svdDeletionStatus::nonPositiveSemidefinite );
        }

        result.m_storage->m_clampedEigenvalues = 0;
        for( Eigen::Index index = 0; index < deletedCount; ++index )
        {
            if( workspace.m_storage->m_eigenvalues( index ) < realT( 0 ) )
            {
                workspace.m_storage->m_eigenvalues( index ) = realT( 0 );
                ++result.m_storage->m_clampedEigenvalues;
            }
        }

        auto complementRoot =
            workspace.m_storage->m_complementRoot.matrix().topLeftCorner( deletedCount, deletedCount );
        for( Eigen::Index row = 0; row < deletedCount; ++row )
        {
            const realT root = std::sqrt( workspace.m_storage->m_eigenvalues( row ) );
            for( Eigen::Index column = 0; column < deletedCount; ++column )
            {
                complementRoot( row, column ) = root * workspace.m_storage->m_eigenvectors( column, row );
            }
        }

        auto topCore = workspace.m_storage->m_stableCore.matrix().topRows( rank );
        topCore.setIdentity();
        topCore.noalias() -= factorRows.transpose() * factorRows;
        workspace.m_storage->m_stableCore.matrix().middleRows( rank, deletedCount ).noalias() =
            complementRoot * factorRows;
        for( Eigen::Index column = 0; column < rank; ++column )
        {
            workspace.m_storage->m_stableCore.matrix().col( column ).head( rank + deletedCount ) *=
                workspace.m_storage->m_normalizedSingularValues( column );
        }

        realT unusedLeft{ 0 };
        const MXLAPACK_INT coreRows = static_cast<MXLAPACK_INT>( rank + deletedCount );
        const MXLAPACK_INT lapackRank = static_cast<MXLAPACK_INT>( rank );
        info = callGesvd<realT>( 'N',
                                 'S',
                                 coreRows,
                                 lapackRank,
                                 workspace.m_storage->m_stableCore.data(),
                                 static_cast<MXLAPACK_INT>( workspace.m_storage->m_stableCore.rows() ),
                                 workspace.m_storage->m_stableSingularValues.data(),
                                 &unusedLeft,
                                 1,
                                 workspace.m_storage->m_rightTranspose.data(),
                                 lapackRank,
                                 workspace.m_storage->m_gesvdWork.data(),
                                 static_cast<MXLAPACK_INT>( workspace.m_storage->m_gesvdWork.size() ) );
        if( info != 0 )
        {
            return fail( result, svdDeletionStatus::solverFailure, info );
        }
        if( !allFinite( workspace.m_storage->m_stableSingularValues ) ||
            !allFinite( workspace.m_storage->m_rightTranspose ) )
        {
            return fail( result, svdDeletionStatus::nonFiniteOutput );
        }
        if( !validDescendingSpectrum( workspace.m_storage->m_stableSingularValues ) )
        {
            return fail( result, svdDeletionStatus::invalidSolverOutput );
        }

        for( Eigen::Index output = 0; output < rank; ++output )
        {
            const realT normalizedValue = workspace.m_storage->m_stableSingularValues( output );
            if( !rescaleSingularValue( result.m_storage->m_singularValues( output ),
                                       result.m_storage->m_squaredSingularValues( output ),
                                       normalizedValue,
                                       scale ) )
            {
                return fail( result, svdDeletionStatus::rescalingOverflow );
            }
            if( output < outputRank )
            {
                result.m_storage->m_rotation.matrix().col( output ) =
                    workspace.m_storage->m_rightTranspose.matrix().row( output ).transpose();
            }
        }

        result.m_storage->m_backend = svdDeletionBackend::stableCore;
        result.m_storage->m_lapackInfo = 0;
        result.m_storage->m_status = result.m_storage->m_clampedEigenvalues > 0 ? svdDeletionStatus::successWithClamping
                                                                                : svdDeletionStatus::success;
        return result.m_storage->m_status;
    }

    // Execute the complement-preserving small-SVD deletion core after workspace loading.
    static svdDeletionStatus stable( resultT &result,                                 /* [out] updated system */
                                     svdDeletionConstVectorRef<realT> singularValues, /* [in] base values */
                                     svdDeletionConstMatrixRef<realT> deletedRows,    /* [in] deleted rows */
                                     Eigen::Index outputRank,                         /* [in] requested rank */
                                     workspaceT &workspace /* [in,out] scratch */ )
    {
        if( !validCoreInputs( singularValues, deletedRows, outputRank ) )
        {
            return fail( result, svdDeletionStatus::invalidInput );
        }
        if( deletedRows.rows() == 0 )
        {
            return identity( result, singularValues, outputRank, svdDeletionBackend::stableCore );
        }

        const Eigen::Index rank = singularValues.size();
        svdDeletionStatus status =
            prepare( result, workspace, rank, deletedRows.rows(), outputRank, svdDeletionBackend::stableCore );
        if( status != svdDeletionStatus::success )
        {
            return status;
        }
        loadDeletedRows( workspace, deletedRows );
        return stableLoaded( result, singularValues, deletedRows.rows(), outputRank, workspace );
    }

    // Gather arbitrary factor rows and dispatch to the selected core.
    static svdDeletionStatus remove( resultT &result,                                 /* [out] updated system */
                                     svdDeletionConstVectorRef<realT> singularValues, /* [in] base values */
                                     svdDeletionConstMatrixRef<realT> factor,         /* [in] deleted-side factor */
                                     std::span<const Eigen::Index> indices,           /* [in] rows to remove */
                                     Eigen::Index outputRank,                         /* [in] requested rank */
                                     workspaceT &workspace,                           /* [in,out] scratch */
                                     svdDeletionBackend backend /* [in] selected backend */ )
    {
        const Eigen::Index rank = singularValues.size();
        if( !validSingularInputs( singularValues, outputRank ) || !validBackend( backend ) || factor.rows() <= 0 ||
            factor.rows() < rank || factor.cols() != rank || !validLapackDimension( factor.rows() ) ||
            indices.size() >= static_cast<std::size_t>( factor.rows() ) ||
            indices.size() > static_cast<std::size_t>( std::numeric_limits<MXLAPACK_INT>::max() - rank ) )
        {
            return fail( result, svdDeletionStatus::invalidInput );
        }

        Eigen::Index previous = -1;
        for( const Eigen::Index index : indices )
        {
            if( index < 0 || index >= factor.rows() || index <= previous )
            {
                return fail( result, svdDeletionStatus::invalidInput );
            }
            previous = index;
        }

        if( indices.empty() )
        {
            return identity( result, singularValues, outputRank, backend );
        }

        const Eigen::Index deletedCount = static_cast<Eigen::Index>( indices.size() );
        const svdDeletionStatus status = prepare( result, workspace, rank, deletedCount, outputRank, backend );
        if( status != svdDeletionStatus::success )
        {
            return status;
        }

        for( Eigen::Index row = 0; row < static_cast<Eigen::Index>( indices.size() ); ++row )
        {
            for( Eigen::Index column = 0; column < rank; ++column )
            {
                const realT value = factor( indices[row], column );
                if( !math::isFinite( value ) )
                {
                    return fail( result, svdDeletionStatus::invalidInput );
                }
                workspace.m_storage->m_deletedRows( row, column ) = value;
            }
        }

        return backend == svdDeletionBackend::leadingCovariance
                   ? leadingLoaded( result, singularValues, deletedCount, outputRank, workspace )
                   : stableLoaded( result, singularValues, deletedCount, outputRank, workspace );
    }
};

} // namespace detail

const char *svdDeletionBackendName( svdDeletionBackend backend )
{
    switch( backend )
    {
        case svdDeletionBackend::leadingCovariance:
            return "leadingCovariance";
        case svdDeletionBackend::stableCore:
            return "stableCore";
    }
    return "unknown";
}

const char *svdDeletionStatusName( svdDeletionStatus status )
{
    switch( status )
    {
        case svdDeletionStatus::notComputed:
            return "notComputed";
        case svdDeletionStatus::success:
            return "success";
        case svdDeletionStatus::successWithClamping:
            return "successWithClamping";
        case svdDeletionStatus::invalidInput:
            return "invalidInput";
        case svdDeletionStatus::allocationFailure:
            return "allocationFailure";
        case svdDeletionStatus::workspaceQueryFailure:
            return "workspaceQueryFailure";
        case svdDeletionStatus::solverFailure:
            return "solverFailure";
        case svdDeletionStatus::nonFiniteOutput:
            return "nonFiniteOutput";
        case svdDeletionStatus::invalidSolverOutput:
            return "invalidSolverOutput";
        case svdDeletionStatus::rescalingOverflow:
            return "rescalingOverflow";
        case svdDeletionStatus::nonPositiveSemidefinite:
            return "nonPositiveSemidefinite";
        case svdDeletionStatus::factorNotOrthonormal:
            return "factorNotOrthonormal";
    }
    return "unknown";
}

bool svdDeletionSucceeded( svdDeletionStatus status ) noexcept
{
    return status == svdDeletionStatus::success || status == svdDeletionStatus::successWithClamping;
}

template <typename realT>
svdDeletionResult<realT>::svdDeletionResult() : m_storage( std::make_unique<storage>() )
{
}

template <typename realT>
svdDeletionResult<realT>::~svdDeletionResult() = default;

template <typename realT>
svdDeletionResult<realT>::svdDeletionResult( svdDeletionResult &&other ) noexcept = default;

template <typename realT>
svdDeletionResult<realT> &svdDeletionResult<realT>::operator=( svdDeletionResult &&other ) noexcept = default;

template <typename realT>
bool svdDeletionResult<realT>::ensureStorage() noexcept
{
    if( m_storage )
    {
        return true;
    }

    try
    {
        detail::callOperationHook<realT>( detail::svdDeletionTestOperation::prepareResult );
        m_storage = std::make_unique<storage>();
    }
    catch( const std::bad_alloc & )
    {
        return false;
    }
    catch( const std::length_error & )
    {
        return false;
    }
    return true;
}

template <typename realT>
svdDeletionStatus svdDeletionResult<realT>::prepare( Eigen::Index baseRank, Eigen::Index outputRank )
{
    if( !ensureStorage() )
    {
        return svdDeletionStatus::allocationFailure;
    }
    auto &stored = *m_storage;

    if( baseRank <= 0 || outputRank <= 0 || outputRank > baseRank || !detail::validLapackDimension( baseRank ) ||
        !detail::validArrayShape<realT>( baseRank, outputRank ) || !detail::validArrayShape<realT>( outputRank, 1 ) )
    {
        stored.m_status = svdDeletionStatus::invalidInput;
        return stored.m_status;
    }

    try
    {
        if( stored.m_baseRank != baseRank || stored.m_maximumOutputRank < outputRank )
        {
            detail::callOperationHook<realT>( detail::svdDeletionTestOperation::prepareResult );
            stored.m_singularValues.resize( baseRank );
            stored.m_squaredSingularValues.resize( baseRank );
            stored.m_rotation.resize( baseRank, outputRank );
            stored.m_maximumOutputRank = outputRank;
        }
    }
    catch( const std::bad_alloc & )
    {
        stored.m_status = svdDeletionStatus::allocationFailure;
        return stored.m_status;
    }
    catch( const std::length_error & )
    {
        stored.m_status = svdDeletionStatus::allocationFailure;
        return stored.m_status;
    }

    stored.m_baseRank = baseRank;
    stored.m_outputRank = outputRank;
    stored.m_status = svdDeletionStatus::notComputed;
    stored.m_lapackInfo = 0;
    stored.m_clampedEigenvalues = 0;
    stored.m_minimumPSDValue = realT( 0 );
    return svdDeletionStatus::success;
}

template <typename realT>
svdDeletionStatus svdDeletionResult<realT>::status() const noexcept
{
    return m_storage ? m_storage->m_status : svdDeletionStatus::notComputed;
}

template <typename realT>
svdDeletionBackend svdDeletionResult<realT>::backend() const noexcept
{
    return m_storage ? m_storage->m_backend : svdDeletionBackend::stableCore;
}

template <typename realT>
svdDeletionConstVectorRef<realT> svdDeletionResult<realT>::singularValues() const noexcept
{
    static const svdDeletionVector<realT> empty;
    return m_storage ? m_storage->m_singularValues : empty;
}

template <typename realT>
svdDeletionConstVectorRef<realT> svdDeletionResult<realT>::squaredSingularValues() const noexcept
{
    static const svdDeletionVector<realT> empty;
    return m_storage ? m_storage->m_squaredSingularValues : empty;
}

template <typename realT>
svdDeletionConstMatrixRef<realT> svdDeletionResult<realT>::rotation() const noexcept
{
    static const svdDeletionMatrix<realT> empty;
    if( !m_storage )
    {
        return empty;
    }
    return m_storage->m_rotation.matrix().leftCols( m_storage->m_outputRank ).array();
}

template <typename realT>
Eigen::Index svdDeletionResult<realT>::baseRank() const noexcept
{
    return m_storage ? m_storage->m_baseRank : 0;
}

template <typename realT>
Eigen::Index svdDeletionResult<realT>::outputRank() const noexcept
{
    return m_storage ? m_storage->m_outputRank : 0;
}

template <typename realT>
Eigen::Index svdDeletionResult<realT>::maximumOutputRank() const noexcept
{
    return m_storage ? m_storage->m_maximumOutputRank : 0;
}

template <typename realT>
Eigen::Index svdDeletionResult<realT>::clampedEigenvalues() const noexcept
{
    return m_storage ? m_storage->m_clampedEigenvalues : 0;
}

template <typename realT>
realT svdDeletionResult<realT>::minimumPSDValue() const noexcept
{
    return m_storage ? m_storage->m_minimumPSDValue : realT( 0 );
}

template <typename realT>
MXLAPACK_INT svdDeletionResult<realT>::lapackInfo() const noexcept
{
    return m_storage ? m_storage->m_lapackInfo : 0;
}

template <typename realT>
svdDeletionWorkspace<realT>::svdDeletionWorkspace() : m_storage( std::make_unique<storage>() )
{
}

template <typename realT>
svdDeletionWorkspace<realT>::~svdDeletionWorkspace()
{
}

template <typename realT>
svdDeletionWorkspace<realT>::svdDeletionWorkspace( svdDeletionWorkspace &&other ) noexcept = default;

template <typename realT>
svdDeletionWorkspace<realT> &svdDeletionWorkspace<realT>::operator=( svdDeletionWorkspace &&other ) noexcept = default;

template <typename realT>
bool svdDeletionWorkspace<realT>::ensureStorage() noexcept
{
    if( m_storage )
    {
        return true;
    }

    try
    {
        detail::callOperationHook<realT>( detail::svdDeletionTestOperation::prepareWorkspace );
        m_storage = std::make_unique<storage>();
    }
    catch( const std::bad_alloc & )
    {
        return false;
    }
    catch( const std::length_error & )
    {
        return false;
    }
    return true;
}

template <typename realT>
svdDeletionStatus
svdDeletionWorkspace<realT>::prepare( Eigen::Index baseRank, Eigen::Index maximumDeleted, svdDeletionBackend backend )
{
    if( !ensureStorage() )
    {
        return svdDeletionStatus::allocationFailure;
    }
    auto &stored = *m_storage;
    stored.m_lapackInfo = 0;
    if( baseRank <= 0 || maximumDeleted < 0 ||
        ( backend != svdDeletionBackend::leadingCovariance && backend != svdDeletionBackend::stableCore ) ||
        !detail::validLapackDimension( baseRank ) || !detail::validLapackDimension( maximumDeleted ) ||
        maximumDeleted > std::numeric_limits<MXLAPACK_INT>::max() - baseRank )
    {
        return svdDeletionStatus::invalidInput;
    }
    if( stored.m_prepared && stored.m_baseRank == baseRank && stored.m_maximumDeleted >= maximumDeleted &&
        stored.m_backend == backend )
    {
        return svdDeletionStatus::success;
    }

    const Eigen::Index eigenSize = backend == svdDeletionBackend::leadingCovariance ? baseRank : maximumDeleted;
    if( eigenSize > std::numeric_limits<Eigen::Index>::max() / 2 )
    {
        return svdDeletionStatus::invalidInput; // LCOV_EXCL_LINE -- requires an Eigen index narrower than LAPACK.
    }
    const Eigen::Index supportSize = 2 * std::max<Eigen::Index>( 1, eigenSize );
    Eigen::Index minimumGesvdWork{ 0 };
    bool validShapes =
        detail::validArrayShape<realT>( maximumDeleted, baseRank ) && detail::validArrayShape<realT>( baseRank, 1 ) &&
        detail::validArrayShape<realT>( eigenSize, eigenSize ) && detail::validArrayShape<realT>( eigenSize, 1 ) &&
        detail::validVectorSize<MXLAPACK_INT>( supportSize );
    if( backend == svdDeletionBackend::leadingCovariance )
    {
        validShapes = validShapes && detail::validArrayShape<realT>( baseRank, maximumDeleted );
    }
    else
    {
        validShapes = validShapes && detail::validArrayShape<realT>( maximumDeleted, maximumDeleted ) &&
                      detail::validArrayShape<realT>( baseRank + maximumDeleted, baseRank ) &&
                      detail::validArrayShape<realT>( baseRank, baseRank );
        const unsigned long long rank = static_cast<unsigned long long>( baseRank );
        const unsigned long long deleted = static_cast<unsigned long long>( maximumDeleted );
        const unsigned long long lapackMaximum =
            static_cast<unsigned long long>( std::numeric_limits<MXLAPACK_INT>::max() );
        validShapes = validShapes && rank <= lapackMaximum / 5 && rank <= ( lapackMaximum - deleted ) / 4;
        if( validShapes )
        {
            minimumGesvdWork = std::max<Eigen::Index>( 1, std::max( 5 * baseRank, 4 * baseRank + maximumDeleted ) );
            validShapes = detail::validVectorSize<realT>( minimumGesvdWork );
        }
    }
    if( !validShapes )
    {
        return svdDeletionStatus::invalidInput;
    }

    try
    {
        detail::callOperationHook<realT>( detail::svdDeletionTestOperation::prepareWorkspace );

        stored.m_deletedRows.resize( maximumDeleted, baseRank );
        stored.m_normalizedSingularValues.resize( baseRank );
        stored.m_scaledDeletedTranspose.resize( 0, 0 );
        stored.m_complementRoot.resize( 0, 0 );
        stored.m_stableCore.resize( 0, 0 );
        stored.m_rightTranspose.resize( 0, 0 );
        stored.m_stableSingularValues.resize( 0 );

        if( backend == svdDeletionBackend::leadingCovariance )
        {
            stored.m_scaledDeletedTranspose.resize( baseRank, maximumDeleted );
            stored.m_psdCore.resize( baseRank, baseRank );
        }
        else
        {
            stored.m_psdCore.resize( maximumDeleted, maximumDeleted );
            stored.m_complementRoot.resize( maximumDeleted, maximumDeleted );
            stored.m_stableCore.resize( baseRank + maximumDeleted, baseRank );
            stored.m_rightTranspose.resize( baseRank, baseRank );
            stored.m_stableSingularValues.resize( baseRank );
        }

        stored.m_eigenvectors.resize( eigenSize, eigenSize );
        stored.m_eigenvalues.resize( eigenSize );
        stored.m_support.resize( static_cast<std::size_t>( supportSize ) );
        stored.m_syevrWork.clear();
        stored.m_syevrIWork.clear();
        stored.m_gesvdWork.clear();

        if( eigenSize > 0 )
        {
            stored.m_psdCore.matrix().setIdentity();
            stored.m_eigenvectors.setZero();
            stored.m_eigenvalues.setZero();
            realT workQuery{ 0 };
            MXLAPACK_INT integerWorkQuery{ 0 };
            MXLAPACK_INT found{ 0 };
            const MXLAPACK_INT lapackEigenSize = static_cast<MXLAPACK_INT>( eigenSize );
            const MXLAPACK_INT info = detail::callSyevr<realT>( 'V',
                                                                'A',
                                                                'L',
                                                                lapackEigenSize,
                                                                stored.m_psdCore.data(),
                                                                lapackEigenSize,
                                                                realT( 0 ),
                                                                realT( 0 ),
                                                                1,
                                                                lapackEigenSize,
                                                                math::lamch<realT>( 'S' ),
                                                                &found,
                                                                stored.m_eigenvalues.data(),
                                                                stored.m_eigenvectors.data(),
                                                                lapackEigenSize,
                                                                stored.m_support.data(),
                                                                &workQuery,
                                                                -1,
                                                                &integerWorkQuery,
                                                                -1 );
            MXLAPACK_INT realWorkSize{ 0 };
            if( info != 0 || !detail::querySize( realWorkSize, workQuery ) || integerWorkQuery < 1 ||
                !detail::validVectorSize<realT>( realWorkSize ) ||
                !detail::validVectorSize<MXLAPACK_INT>( integerWorkQuery ) )
            {
                clear();
                stored.m_lapackInfo = info;
                return svdDeletionStatus::workspaceQueryFailure;
            }
            stored.m_syevrWork.resize( static_cast<std::size_t>( realWorkSize ) );
            stored.m_syevrIWork.resize( static_cast<std::size_t>( integerWorkQuery ) );
        }

        if( backend == svdDeletionBackend::stableCore )
        {
            stored.m_stableCore.setZero();
            stored.m_rightTranspose.setZero();
            stored.m_stableSingularValues.setZero();
            realT unusedLeft{ 0 };
            realT workQuery{ 0 };
            const MXLAPACK_INT rows = static_cast<MXLAPACK_INT>( baseRank + maximumDeleted );
            const MXLAPACK_INT cols = static_cast<MXLAPACK_INT>( baseRank );
            const MXLAPACK_INT info = detail::callGesvd<realT>( 'N',
                                                                'S',
                                                                rows,
                                                                cols,
                                                                stored.m_stableCore.data(),
                                                                rows,
                                                                stored.m_stableSingularValues.data(),
                                                                &unusedLeft,
                                                                1,
                                                                stored.m_rightTranspose.data(),
                                                                cols,
                                                                &workQuery,
                                                                -1 );
            MXLAPACK_INT workSize{ 0 };
            if( info != 0 || !detail::querySize( workSize, workQuery ) )
            {
                clear();
                stored.m_lapackInfo = info;
                return svdDeletionStatus::workspaceQueryFailure;
            }
            workSize = std::max( workSize, static_cast<MXLAPACK_INT>( minimumGesvdWork ) );
            if( !detail::validVectorSize<realT>( workSize ) )
            {
                clear(); // LCOV_EXCL_LINE -- requires size_t narrower than the configured LAPACK integer ABI.
                return svdDeletionStatus::workspaceQueryFailure; // LCOV_EXCL_LINE
            }
            stored.m_gesvdWork.resize( static_cast<std::size_t>( workSize ) );
        }
    }
    catch( const std::bad_alloc & )
    {
        clear();
        return svdDeletionStatus::allocationFailure;
    }
    catch( const std::length_error & )
    {
        clear();
        return svdDeletionStatus::allocationFailure;
    }

    stored.m_baseRank = baseRank;
    stored.m_maximumDeleted = maximumDeleted;
    stored.m_backend = backend;
    stored.m_prepared = true;
    return svdDeletionStatus::success;
}

template <typename realT>
void svdDeletionWorkspace<realT>::clear() noexcept
{
    if( !m_storage )
    {
        return;
    }

    m_storage->m_deletedRows.resize( 0, 0 );
    m_storage->m_normalizedSingularValues.resize( 0 );
    m_storage->m_scaledDeletedTranspose.resize( 0, 0 );
    m_storage->m_psdCore.resize( 0, 0 );
    m_storage->m_eigenvectors.resize( 0, 0 );
    m_storage->m_eigenvalues.resize( 0 );
    m_storage->m_complementRoot.resize( 0, 0 );
    m_storage->m_stableCore.resize( 0, 0 );
    m_storage->m_rightTranspose.resize( 0, 0 );
    m_storage->m_stableSingularValues.resize( 0 );
    m_storage->m_support.clear();
    m_storage->m_syevrWork.clear();
    m_storage->m_syevrIWork.clear();
    m_storage->m_gesvdWork.clear();
    m_storage->m_baseRank = 0;
    m_storage->m_maximumDeleted = 0;
    m_storage->m_backend = svdDeletionBackend::stableCore;
    m_storage->m_prepared = false;
    m_storage->m_lapackInfo = 0;
}

template <typename realT>
bool svdDeletionWorkspace<realT>::prepared() const noexcept
{
    return m_storage && m_storage->m_prepared;
}

template <typename realT>
Eigen::Index svdDeletionWorkspace<realT>::baseRank() const noexcept
{
    return m_storage ? m_storage->m_baseRank : 0;
}

template <typename realT>
Eigen::Index svdDeletionWorkspace<realT>::maximumDeleted() const noexcept
{
    return m_storage ? m_storage->m_maximumDeleted : 0;
}

template <typename realT>
svdDeletionBackend svdDeletionWorkspace<realT>::backend() const noexcept
{
    return m_storage ? m_storage->m_backend : svdDeletionBackend::stableCore;
}

template <typename realT>
MXLAPACK_INT svdDeletionWorkspace<realT>::lapackInfo() const noexcept
{
    return m_storage ? m_storage->m_lapackInfo : 0;
}

template <typename realT>
svdDeletionStatus validateSvdDeletionFactorImpl( svdDeletionConstMatrixRef<realT> factor, realT tolerance )
{
    if( factor.rows() < factor.cols() || factor.cols() <= 0 || !detail::validLapackDimension( factor.rows() ) ||
        !detail::validLapackDimension( factor.cols() ) || !math::isFinite( tolerance ) || tolerance < realT( 0 ) ||
        !detail::validArrayShape<realT>( factor.cols(), factor.cols() ) || !detail::allFinite( factor ) )
    {
        return svdDeletionStatus::invalidInput;
    }

    try
    {
        detail::callOperationHook<realT>( detail::svdDeletionTestOperation::validateFactor );
        svdDeletionMatrix<realT> gram( factor.cols(), factor.cols() );
        gram.matrix().noalias() = factor.matrix().transpose() * factor.matrix();
        if( !detail::allFinite( gram ) )
        {
            return svdDeletionStatus::invalidInput;
        }
        for( Eigen::Index index = 0; index < gram.rows(); ++index )
        {
            gram( index, index ) -= realT( 1 );
        }

        if( tolerance == realT( 0 ) )
        {
            tolerance = realT( 64 ) * std::numeric_limits<realT>::epsilon() *
                        static_cast<realT>( std::max( factor.rows(), factor.cols() ) );
        }
        if( gram.abs().maxCoeff() > tolerance )
        {
            return svdDeletionStatus::factorNotOrthonormal;
        }
    }
    catch( const std::bad_alloc & )
    {
        return svdDeletionStatus::allocationFailure;
    }
    catch( const std::length_error & )
    {
        return svdDeletionStatus::allocationFailure;
    }

    return svdDeletionStatus::success;
}

svdDeletionStatus validateSvdDeletionFactor( svdDeletionConstMatrixRef<float> factor, float tolerance )
{
    return validateSvdDeletionFactorImpl<float>( factor, tolerance );
}

svdDeletionStatus validateSvdDeletionFactor( svdDeletionConstMatrixRef<double> factor, double tolerance )
{
    return validateSvdDeletionFactorImpl<double>( factor, tolerance );
}

template <typename realT>
svdDeletionStatus svdDeletionLeadingCore( svdDeletionResult<realT> &result,
                                          std::type_identity_t<svdDeletionConstVectorRef<realT>> singularValues,
                                          std::type_identity_t<svdDeletionConstMatrixRef<realT>> deletedRows,
                                          Eigen::Index outputRank,
                                          svdDeletionWorkspace<realT> &workspace )
{
    return detail::svdDeletionImplementation<realT>::leading( result,
                                                              singularValues,
                                                              deletedRows,
                                                              outputRank,
                                                              workspace );
}

template <typename realT>
svdDeletionStatus svdDeletionStableCore( svdDeletionResult<realT> &result,
                                         std::type_identity_t<svdDeletionConstVectorRef<realT>> singularValues,
                                         std::type_identity_t<svdDeletionConstMatrixRef<realT>> deletedRows,
                                         Eigen::Index outputRank,
                                         svdDeletionWorkspace<realT> &workspace )
{
    return detail::svdDeletionImplementation<realT>::stable( result,
                                                             singularValues,
                                                             deletedRows,
                                                             outputRank,
                                                             workspace );
}

template <typename realT>
svdDeletionStatus svdDeletionCore( svdDeletionResult<realT> &result,
                                   std::type_identity_t<svdDeletionConstVectorRef<realT>> singularValues,
                                   std::type_identity_t<svdDeletionConstMatrixRef<realT>> deletedRows,
                                   Eigen::Index outputRank,
                                   svdDeletionWorkspace<realT> &workspace,
                                   svdDeletionBackend backend )
{
    if( backend != svdDeletionBackend::leadingCovariance && backend != svdDeletionBackend::stableCore )
    {
        return detail::svdDeletionImplementation<realT>::fail( result, svdDeletionStatus::invalidInput );
    }
    return backend == svdDeletionBackend::leadingCovariance
               ? detail::svdDeletionImplementation<realT>::leading( result,
                                                                    singularValues,
                                                                    deletedRows,
                                                                    outputRank,
                                                                    workspace )
               : detail::svdDeletionImplementation<realT>::stable( result,
                                                                   singularValues,
                                                                   deletedRows,
                                                                   outputRank,
                                                                   workspace );
}

template <typename realT>
svdDeletionStatus svdRemoveRows( svdDeletionResult<realT> &result,
                                 std::type_identity_t<svdDeletionConstVectorRef<realT>> singularValues,
                                 std::type_identity_t<svdDeletionConstMatrixRef<realT>> leftFactor,
                                 std::span<const Eigen::Index> deletedIndices,
                                 Eigen::Index outputRank,
                                 svdDeletionWorkspace<realT> &workspace,
                                 svdDeletionBackend backend )
{
    return detail::svdDeletionImplementation<realT>::remove( result,
                                                             singularValues,
                                                             leftFactor,
                                                             deletedIndices,
                                                             outputRank,
                                                             workspace,
                                                             backend );
}

template <typename realT>
svdDeletionStatus svdRemoveColumns( svdDeletionResult<realT> &result,
                                    std::type_identity_t<svdDeletionConstVectorRef<realT>> singularValues,
                                    std::type_identity_t<svdDeletionConstMatrixRef<realT>> rightFactor,
                                    std::span<const Eigen::Index> deletedIndices,
                                    Eigen::Index outputRank,
                                    svdDeletionWorkspace<realT> &workspace,
                                    svdDeletionBackend backend )
{
    return detail::svdDeletionImplementation<realT>::remove( result,
                                                             singularValues,
                                                             rightFactor,
                                                             deletedIndices,
                                                             outputRank,
                                                             workspace,
                                                             backend );
}

template class svdDeletionResult<float>;
template class svdDeletionResult<double>;
template class svdDeletionWorkspace<float>;
template class svdDeletionWorkspace<double>;

template svdDeletionStatus svdDeletionLeadingCore<float>( svdDeletionResult<float> &,
                                                          std::type_identity_t<svdDeletionConstVectorRef<float>>,
                                                          std::type_identity_t<svdDeletionConstMatrixRef<float>>,
                                                          Eigen::Index,
                                                          svdDeletionWorkspace<float> & );
template svdDeletionStatus svdDeletionLeadingCore<double>( svdDeletionResult<double> &,
                                                           std::type_identity_t<svdDeletionConstVectorRef<double>>,
                                                           std::type_identity_t<svdDeletionConstMatrixRef<double>>,
                                                           Eigen::Index,
                                                           svdDeletionWorkspace<double> & );

template svdDeletionStatus svdDeletionStableCore<float>( svdDeletionResult<float> &,
                                                         std::type_identity_t<svdDeletionConstVectorRef<float>>,
                                                         std::type_identity_t<svdDeletionConstMatrixRef<float>>,
                                                         Eigen::Index,
                                                         svdDeletionWorkspace<float> & );
template svdDeletionStatus svdDeletionStableCore<double>( svdDeletionResult<double> &,
                                                          std::type_identity_t<svdDeletionConstVectorRef<double>>,
                                                          std::type_identity_t<svdDeletionConstMatrixRef<double>>,
                                                          Eigen::Index,
                                                          svdDeletionWorkspace<double> & );

template svdDeletionStatus svdDeletionCore<float>( svdDeletionResult<float> &,
                                                   std::type_identity_t<svdDeletionConstVectorRef<float>>,
                                                   std::type_identity_t<svdDeletionConstMatrixRef<float>>,
                                                   Eigen::Index,
                                                   svdDeletionWorkspace<float> &,
                                                   svdDeletionBackend );
template svdDeletionStatus svdDeletionCore<double>( svdDeletionResult<double> &,
                                                    std::type_identity_t<svdDeletionConstVectorRef<double>>,
                                                    std::type_identity_t<svdDeletionConstMatrixRef<double>>,
                                                    Eigen::Index,
                                                    svdDeletionWorkspace<double> &,
                                                    svdDeletionBackend );

template svdDeletionStatus svdRemoveRows<float>( svdDeletionResult<float> &,
                                                 std::type_identity_t<svdDeletionConstVectorRef<float>>,
                                                 std::type_identity_t<svdDeletionConstMatrixRef<float>>,
                                                 std::span<const Eigen::Index>,
                                                 Eigen::Index,
                                                 svdDeletionWorkspace<float> &,
                                                 svdDeletionBackend );
template svdDeletionStatus svdRemoveRows<double>( svdDeletionResult<double> &,
                                                  std::type_identity_t<svdDeletionConstVectorRef<double>>,
                                                  std::type_identity_t<svdDeletionConstMatrixRef<double>>,
                                                  std::span<const Eigen::Index>,
                                                  Eigen::Index,
                                                  svdDeletionWorkspace<double> &,
                                                  svdDeletionBackend );

template svdDeletionStatus svdRemoveColumns<float>( svdDeletionResult<float> &,
                                                    std::type_identity_t<svdDeletionConstVectorRef<float>>,
                                                    std::type_identity_t<svdDeletionConstMatrixRef<float>>,
                                                    std::span<const Eigen::Index>,
                                                    Eigen::Index,
                                                    svdDeletionWorkspace<float> &,
                                                    svdDeletionBackend );
template svdDeletionStatus svdRemoveColumns<double>( svdDeletionResult<double> &,
                                                     std::type_identity_t<svdDeletionConstVectorRef<double>>,
                                                     std::type_identity_t<svdDeletionConstMatrixRef<double>>,
                                                     std::span<const Eigen::Index>,
                                                     Eigen::Index,
                                                     svdDeletionWorkspace<double> &,
                                                     svdDeletionBackend );

template detail::svdDeletionTestHooks<float> &detail::svdDeletionHooks<float>();
template detail::svdDeletionTestHooks<double> &detail::svdDeletionHooks<double>();

} // namespace math
} // namespace mx
