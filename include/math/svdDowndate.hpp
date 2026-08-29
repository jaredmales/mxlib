/** \file
 * \brief Reusable row and column deletion updates for thin singular value decompositions.
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

#ifndef math_svdDowndate_hpp
#define math_svdDowndate_hpp

#include <Eigen/Dense>

#include <memory>
#include <span>
#include <type_traits>

#include "templateLapack.hpp"

namespace mx
{
namespace math
{

/** \addtogroup svd_downdate
 * @{
 */

/// Numerical backend used to delete rows or columns from thin-SVD factors.
enum class svdDeletionBackend
{
    leadingCovariance, ///< Symmetric leading-spectrum core; fastest when small singular values are not required.
    stableCore         ///< Complement-preserving small SVD; avoids squaring singular-value conditioning.
};

/// Completion status for an SVD deletion operation.
enum class svdDeletionStatus
{
    notComputed,             ///< No operation has published a result.
    success,                 ///< The operation completed without numerical clamping.
    successWithClamping,     ///< The operation completed after clamping roundoff-scale negative eigenvalues.
    invalidInput,            ///< Dimensions, values, indices, or requested output rank are invalid.
    allocationFailure,       ///< Result or workspace allocation failed.
    workspaceQueryFailure,   ///< LAPACK returned an invalid or failed workspace query.
    solverFailure,           ///< LAPACK failed during the numerical solve.
    nonFiniteOutput,         ///< LAPACK returned a non-finite singular system.
    invalidSolverOutput,     ///< LAPACK returned a finite spectrum with invalid ordering or sign.
    rescalingOverflow,       ///< A finite normalized result cannot be represented after restoring input scale.
    nonPositiveSemidefinite, ///< A theoretically PSD core has a materially negative eigenvalue.
    factorNotOrthonormal     ///< A requested singular-factor validation failed.
};

/// Return a stable text representation of an SVD deletion backend.
const char *svdDeletionBackendName( svdDeletionBackend backend /**< [in] backend to describe */ );

/// Return a stable text representation of an SVD deletion status.
const char *svdDeletionStatusName( svdDeletionStatus status /**< [in] status to describe */ );

/// Return true when a status represents usable numerical output.
bool svdDeletionSucceeded( svdDeletionStatus status /**< [in] status to classify */ ) noexcept;

/// Column-major dynamic matrix used by the SVD deletion API.
template <typename realT>
using svdDeletionMatrix = Eigen::Array<realT, Eigen::Dynamic, Eigen::Dynamic, Eigen::ColMajor>;

/// Dynamic column vector used by the SVD deletion API.
template <typename realT>
using svdDeletionVector = Eigen::Array<realT, Eigen::Dynamic, 1>;

/// Non-owning read-only reference to a compatible column-major SVD deletion matrix.
template <typename realT>
using svdDeletionConstMatrixRef = Eigen::Ref<const svdDeletionMatrix<realT>>;

/// Non-owning read-only reference to a compatible SVD deletion vector.
template <typename realT>
using svdDeletionConstVectorRef = Eigen::Ref<const svdDeletionVector<realT>>;

template <typename realT>
class svdDeletionResult;

template <typename realT>
class svdDeletionWorkspace;

namespace detail
{

template <typename realT>
struct svdDeletionImplementation;

/// \cond svdDeletion_test_detail

/// Internal stages exposed only for deterministic failure-path tests.
enum class svdDeletionTestOperation
{
    validateFactor,
    prepareResult,
    prepareWorkspace
};

/// Function signature for deterministic allocation-failure injection.
using svdDeletionOperationHookT = void ( * )( svdDeletionTestOperation );

/// Function signature matching the LAPACK SYEVR wrapper.
template <typename realT>
using svdDeletionSyevrHookT = MXLAPACK_INT ( * )( char,
                                                  char,
                                                  char,
                                                  MXLAPACK_INT,
                                                  realT *,
                                                  MXLAPACK_INT,
                                                  realT,
                                                  realT,
                                                  MXLAPACK_INT,
                                                  MXLAPACK_INT,
                                                  realT,
                                                  MXLAPACK_INT *,
                                                  realT *,
                                                  realT *,
                                                  MXLAPACK_INT,
                                                  MXLAPACK_INT *,
                                                  realT *,
                                                  MXLAPACK_INT,
                                                  MXLAPACK_INT *,
                                                  MXLAPACK_INT );

/// Function signature matching the LAPACK GESVD wrapper.
template <typename realT>
using svdDeletionGesvdHookT = MXLAPACK_INT ( * )( char,
                                                  char,
                                                  MXLAPACK_INT,
                                                  MXLAPACK_INT,
                                                  realT *,
                                                  MXLAPACK_INT,
                                                  realT *,
                                                  realT *,
                                                  MXLAPACK_INT,
                                                  realT *,
                                                  MXLAPACK_INT,
                                                  realT *,
                                                  MXLAPACK_INT );

/// Test hooks for one floating-point specialization.
template <typename realT>
struct svdDeletionTestHooks
{
    /// Optional operation hook; null performs the production operation.
    svdDeletionOperationHookT operation{ nullptr };

    /// Optional SYEVR hook; null calls `math::syevr<realT>`.
    svdDeletionSyevrHookT<realT> syevr{ nullptr };

    /// Optional GESVD hook; null calls `math::gesvd<realT>`.
    svdDeletionGesvdHookT<realT> gesvd{ nullptr };
};

/// Access the process-wide SVD deletion test hooks for one scalar type.
template <typename realT>
svdDeletionTestHooks<realT> &svdDeletionHooks();

/// \endcond

} // namespace detail

/// Result of deleting rows or columns from a represented thin SVD.
/** If \f$A=U\Sigma V^T\f$, the returned rotation applies to the factor on the side that was not deleted. For row
 * deletion the updated right factor is `V * rotation()`; for column deletion the updated left factor is
 * `U * rotation()`. Singular values are descending.
 *
 * Output arrays remain allocated across calls. Their numerical contents are valid only when status() reports a
 * successful operation.
 *
 * \tparam realT floating-point type; supported explicit instantiations are float and double.
 */
template <typename realT>
class svdDeletionResult
{
  public:
    /// Construct an empty result.
    svdDeletionResult();

    /// Release result storage using mxlib's Eigen allocation configuration.
    ~svdDeletionResult();

    /// Results cannot be copied across the mxlib ABI boundary.
    svdDeletionResult( const svdDeletionResult &other /**< [in] result that copying is forbidden from */ ) = delete;

    /// Results cannot be copy-assigned across the mxlib ABI boundary.
    svdDeletionResult &
    operator=( const svdDeletionResult &other /**< [in] result that copy assignment is forbidden from */ ) = delete;

    /// Move owned result storage from another result.
    svdDeletionResult( svdDeletionResult &&other /**< [in,out] result to move from */ ) noexcept;

    /// Replace this result by moving owned storage from another result.
    svdDeletionResult &operator=( svdDeletionResult &&other /**< [in,out] result to move from */ ) noexcept;

    /// Prepare or reuse output storage for a base and requested active output rank.
    /** Storage grows when required but does not shrink while `baseRank` is unchanged. Calling this before a hot loop
     * with the largest planned output rank reserves the rotation capacity for later smaller requests.
     */
    svdDeletionStatus prepare( Eigen::Index baseRank, /**< [in] number of supplied singular triplets */
                               Eigen::Index outputRank /**< [in] number of updated triplets to publish */ );

    /// Return the most recent operation status.
    svdDeletionStatus status() const noexcept;

    /// Return the backend that produced the current result.
    svdDeletionBackend backend() const noexcept;

    /// Return an unaligned borrowed view of all `baseRank()` descending updated singular values.
    /** The view remains valid until this result is prepared, assigned, moved, or destroyed. */
    svdDeletionConstVectorRef<realT> singularValues() const noexcept;

    /// Return an unaligned borrowed view of all corresponding descending squared singular values.
    /** The view remains valid until this result is prepared, assigned, moved, or destroyed. */
    svdDeletionConstVectorRef<realT> squaredSingularValues() const noexcept;

    /// Return the preserved-side rotation, with updated directions in columns.
    /** The unaligned borrowed view remains valid until this result is prepared, assigned, moved, or destroyed. */
    svdDeletionConstMatrixRef<realT> rotation() const noexcept;

    /// Return the base factorization rank for which storage is prepared.
    Eigen::Index baseRank() const noexcept;

    /// Return the requested published rank.
    Eigen::Index outputRank() const noexcept;

    /// Return the allocated rotation-column capacity for the current base rank.
    Eigen::Index maximumOutputRank() const noexcept;

    /// Return the number of roundoff-scale negative eigenvalues clamped to zero.
    Eigen::Index clampedEigenvalues() const noexcept;

    /// Return the smallest pre-clamp eigenvalue from the normalized backend PSD validation core.
    /** A non-positive-semidefinite failure preserves the offending value. An empty stable-core deletion uses one as
     * the vacuous minimum of its zero-dimensional complement core.
     */
    realT minimumPSDValue() const noexcept;

    /// Return the underlying LAPACK status from the most recent failed query or solve.
    MXLAPACK_INT lapackInfo() const noexcept;

  private:
    friend struct detail::svdDeletionImplementation<realT>;

    /// Allocate private result storage after a move, translating allocation failure to false.
    bool ensureStorage() noexcept;

    class storage;

    /// Opaque result storage allocated and released with mxlib's Eigen build configuration.
    std::unique_ptr<storage> m_storage;
};

/// Reusable, non-shared storage for SVD deletion operations.
/** Call prepare() before a hot loop to allocate and query LAPACK once. A later operation reuses the storage whenever
 * its dimensions and backend fit the prepared capacity. The workspace is noncopyable but movable, allowing callers
 * to own one independent workspace per worker without sharing in-flight numerical state.
 *
 * \tparam realT floating-point type; supported explicit instantiations are float and double.
 */
template <typename realT>
class svdDeletionWorkspace
{
  public:
    /// Construct an empty workspace.
    svdDeletionWorkspace();

    /// Release owned scratch storage.
    ~svdDeletionWorkspace();

    /// Workspaces cannot be copied.
    svdDeletionWorkspace( const svdDeletionWorkspace &other /**< [in] workspace that copying is forbidden from */ ) =
        delete;

    /// Workspaces cannot be copy-assigned.
    svdDeletionWorkspace &operator=(
        const svdDeletionWorkspace &other /**< [in] workspace that copy assignment is forbidden from */ ) = delete;

    /// Move owned storage and preparation state from another workspace.
    svdDeletionWorkspace( svdDeletionWorkspace &&other /**< [in,out] workspace to move from */ ) noexcept;

    /// Replace this workspace by moving owned storage and preparation state.
    svdDeletionWorkspace &operator=( svdDeletionWorkspace &&other /**< [in,out] workspace to move from */ ) noexcept;

    /// Prepare reusable storage and LAPACK work arrays.
    svdDeletionStatus prepare( Eigen::Index baseRank,       /**< [in] maximum supplied singular rank */
                               Eigen::Index maximumDeleted, /**< [in] maximum rows of the deleted-side factor */
                               svdDeletionBackend backend /**< [in] numerical backend to prepare */ );

    /// Release all prepared storage and reset dimensions.
    void clear() noexcept;

    /// Report whether this workspace has completed preparation.
    bool prepared() const noexcept;

    /// Return the prepared base rank.
    Eigen::Index baseRank() const noexcept;

    /// Return the prepared maximum deletion count.
    Eigen::Index maximumDeleted() const noexcept;

    /// Return the prepared numerical backend.
    svdDeletionBackend backend() const noexcept;

    /// Return the underlying LAPACK status from the most recent failed workspace query.
    MXLAPACK_INT lapackInfo() const noexcept;

  private:
    friend struct detail::svdDeletionImplementation<realT>;

    /// Allocate private workspace storage after a move, translating allocation failure to false.
    bool ensureStorage() noexcept;

    class storage;

    /// Opaque scratch storage allocated and released with mxlib's Eigen build configuration.
    std::unique_ptr<storage> m_storage;
};

/// Validate that a supplied thin singular-vector factor has orthonormal columns.
/** This is an optional one-time base-factor check. Hot deletion calls assume the factor contract and do not repeat
 * this `O(n q^2)` operation. A zero tolerance selects a dimension-scaled default.
 */
svdDeletionStatus
validateSvdDeletionFactor( svdDeletionConstMatrixRef<float> factor, /**< [in] thin singular-vector factor to validate */
                           float tolerance = 0 /**< [in] maximum absolute Gram-matrix error, or zero for automatic */ );

/// Validate that a supplied double-precision thin singular-vector factor has orthonormal columns.
svdDeletionStatus validateSvdDeletionFactor(
    svdDeletionConstMatrixRef<double> factor, /**< [in] thin singular-vector factor to validate */
    double tolerance = 0 /**< [in] maximum absolute Gram-matrix error, or zero for automatic */ );

/// Delete supplied singular-factor rows with the full-spectrum symmetric covariance core.
/** Given deleted-side rows `F`, this solves
 *
 * \f[
 * H = \Sigma (I-F^T F) \Sigma = W \Lambda W^T.
 * \f]
 *
 * The result rotation is `W` and its singular values are `sqrt(diag(Lambda))`. The complete spectrum is evaluated
 * for PSD validation even when only `outputRank` leading directions are published. This backend forms normal
 * equations and therefore does not promise high relative accuracy for the smallest singular values.
 * Classical structured eigensystem downdates are described by \cite bunch_nielsen_1978 and \cite gu_eisenstat_1995;
 * they are candidates for a future non-dense backend behind the same interface.
 */
template <typename realT>
svdDeletionStatus svdDeletionLeadingCore(
    svdDeletionResult<realT> &result, /**< [out] updated spectrum, rotation, and diagnostics */
    std::type_identity_t<svdDeletionConstVectorRef<realT>> singularValues, /**< [in] descending base singular values */
    std::type_identity_t<svdDeletionConstMatrixRef<realT>> deletedRows,    /**< [in] deleted-side factor rows */
    Eigen::Index outputRank, /**< [in] number of leading updated directions to publish */
    svdDeletionWorkspace<realT> &workspace /**< [in,out] reusable, worker-private scratch storage */ );

/// Delete supplied singular-factor rows with the complement-preserving small-SVD core.
/** Let `F` contain deleted-side singular-factor rows and choose `B` so that
 * `B^T B = I - F F^T`. The at-most `(q+c) x q` core
 *
 * \f[
 * K = \begin{bmatrix}(I-F^TF)\Sigma \\ BF\Sigma\end{bmatrix}
 * \f]
 *
 * obeys `K^T K = Sigma (I-F^T F) Sigma`. Its right singular vectors are the preserved-side rotation. This is the
 * default generic backend because it avoids explicitly squaring the represented singular values. See
 * \cite brand_2006 and \cite long_males_2021.
 *
 * Both current backends use dense LAPACK eigensolvers or SVDs and therefore have cubic asymptotic cost in the base
 * rank. The reusable interface intentionally leaves room for a later structured secular-equation backend.
 */
template <typename realT>
svdDeletionStatus svdDeletionStableCore(
    svdDeletionResult<realT> &result, /**< [out] updated spectrum, rotation, and diagnostics */
    std::type_identity_t<svdDeletionConstVectorRef<realT>> singularValues, /**< [in] descending base singular values */
    std::type_identity_t<svdDeletionConstMatrixRef<realT>> deletedRows,    /**< [in] deleted-side factor rows */
    Eigen::Index outputRank, /**< [in] number of leading updated directions to publish */
    svdDeletionWorkspace<realT> &workspace /**< [in,out] reusable, worker-private scratch storage */ );

/// Delete supplied singular-factor rows with an explicitly selected backend.
template <typename realT>
svdDeletionStatus svdDeletionCore(
    svdDeletionResult<realT> &result, /**< [out] updated spectrum, rotation, and diagnostics */
    std::type_identity_t<svdDeletionConstVectorRef<realT>> singularValues, /**< [in] descending base singular values */
    std::type_identity_t<svdDeletionConstMatrixRef<realT>> deletedRows,    /**< [in] deleted-side factor rows */
    Eigen::Index outputRank,                /**< [in] number of leading updated directions to publish */
    svdDeletionWorkspace<realT> &workspace, /**< [in,out] reusable, worker-private scratch storage */
    svdDeletionBackend backend = svdDeletionBackend::stableCore /**< [in] numerical backend */ );

/// Delete physical rows from the matrix represented by a thin SVD.
/** For `A=U Sigma V^T`, this gathers `U[deletedIndices,:]`. The returned rotation applies to `V`. Supplying complete
 * factors makes the deletion identical to a direct SVD of the physically retained matrix. Supplying truncated
 * factors deletes rows exactly from that represented low-rank matrix.
 */
template <typename realT>
svdDeletionStatus svdRemoveRows(
    svdDeletionResult<realT> &result, /**< [out] updated spectrum, rotation, and diagnostics */
    std::type_identity_t<svdDeletionConstVectorRef<realT>> singularValues, /**< [in] descending base singular values */
    std::type_identity_t<svdDeletionConstMatrixRef<realT>> leftFactor,     /**< [in] thin left factor `U` */
    std::span<const Eigen::Index> deletedIndices, /**< [in] sorted, unique row indices to delete */
    Eigen::Index outputRank,                      /**< [in] number of leading updated directions */
    svdDeletionWorkspace<realT> &workspace,       /**< [in,out] reusable, worker-private scratch */
    svdDeletionBackend backend = svdDeletionBackend::stableCore /**< [in] numerical backend */ );

/// Delete physical columns from the matrix represented by a thin SVD.
/** For `A=U Sigma V^T`, this gathers `V[deletedIndices,:]`. The returned rotation applies to `U`. Supplying complete
 * factors makes the deletion identical to a direct SVD of the physically retained matrix. Supplying truncated
 * factors deletes columns exactly from that represented low-rank matrix.
 */
template <typename realT>
svdDeletionStatus svdRemoveColumns(
    svdDeletionResult<realT> &result, /**< [out] updated spectrum, rotation, and diagnostics */
    std::type_identity_t<svdDeletionConstVectorRef<realT>> singularValues, /**< [in] descending base singular values */
    std::type_identity_t<svdDeletionConstMatrixRef<realT>> rightFactor,    /**< [in] thin right factor `V` */
    std::span<const Eigen::Index> deletedIndices, /**< [in] sorted, unique column indices to delete */
    Eigen::Index outputRank,                      /**< [in] number of leading updated directions */
    svdDeletionWorkspace<realT> &workspace,       /**< [in,out] reusable, worker-private scratch */
    svdDeletionBackend backend = svdDeletionBackend::stableCore /**< [in] numerical backend */ );

extern template class svdDeletionResult<float>;
extern template class svdDeletionResult<double>;
extern template class svdDeletionWorkspace<float>;
extern template class svdDeletionWorkspace<double>;

extern template svdDeletionStatus svdDeletionLeadingCore<float>( svdDeletionResult<float> &,
                                                                 std::type_identity_t<svdDeletionConstVectorRef<float>>,
                                                                 std::type_identity_t<svdDeletionConstMatrixRef<float>>,
                                                                 Eigen::Index,
                                                                 svdDeletionWorkspace<float> & );
extern template svdDeletionStatus
svdDeletionLeadingCore<double>( svdDeletionResult<double> &,
                                std::type_identity_t<svdDeletionConstVectorRef<double>>,
                                std::type_identity_t<svdDeletionConstMatrixRef<double>>,
                                Eigen::Index,
                                svdDeletionWorkspace<double> & );

extern template svdDeletionStatus svdDeletionStableCore<float>( svdDeletionResult<float> &,
                                                                std::type_identity_t<svdDeletionConstVectorRef<float>>,
                                                                std::type_identity_t<svdDeletionConstMatrixRef<float>>,
                                                                Eigen::Index,
                                                                svdDeletionWorkspace<float> & );
extern template svdDeletionStatus
svdDeletionStableCore<double>( svdDeletionResult<double> &,
                               std::type_identity_t<svdDeletionConstVectorRef<double>>,
                               std::type_identity_t<svdDeletionConstMatrixRef<double>>,
                               Eigen::Index,
                               svdDeletionWorkspace<double> & );

extern template svdDeletionStatus svdDeletionCore<float>( svdDeletionResult<float> &,
                                                          std::type_identity_t<svdDeletionConstVectorRef<float>>,
                                                          std::type_identity_t<svdDeletionConstMatrixRef<float>>,
                                                          Eigen::Index,
                                                          svdDeletionWorkspace<float> &,
                                                          svdDeletionBackend );
extern template svdDeletionStatus svdDeletionCore<double>( svdDeletionResult<double> &,
                                                           std::type_identity_t<svdDeletionConstVectorRef<double>>,
                                                           std::type_identity_t<svdDeletionConstMatrixRef<double>>,
                                                           Eigen::Index,
                                                           svdDeletionWorkspace<double> &,
                                                           svdDeletionBackend );

extern template svdDeletionStatus svdRemoveRows<float>( svdDeletionResult<float> &,
                                                        std::type_identity_t<svdDeletionConstVectorRef<float>>,
                                                        std::type_identity_t<svdDeletionConstMatrixRef<float>>,
                                                        std::span<const Eigen::Index>,
                                                        Eigen::Index,
                                                        svdDeletionWorkspace<float> &,
                                                        svdDeletionBackend );
extern template svdDeletionStatus svdRemoveRows<double>( svdDeletionResult<double> &,
                                                         std::type_identity_t<svdDeletionConstVectorRef<double>>,
                                                         std::type_identity_t<svdDeletionConstMatrixRef<double>>,
                                                         std::span<const Eigen::Index>,
                                                         Eigen::Index,
                                                         svdDeletionWorkspace<double> &,
                                                         svdDeletionBackend );

extern template svdDeletionStatus svdRemoveColumns<float>( svdDeletionResult<float> &,
                                                           std::type_identity_t<svdDeletionConstVectorRef<float>>,
                                                           std::type_identity_t<svdDeletionConstMatrixRef<float>>,
                                                           std::span<const Eigen::Index>,
                                                           Eigen::Index,
                                                           svdDeletionWorkspace<float> &,
                                                           svdDeletionBackend );
extern template svdDeletionStatus svdRemoveColumns<double>( svdDeletionResult<double> &,
                                                            std::type_identity_t<svdDeletionConstVectorRef<double>>,
                                                            std::type_identity_t<svdDeletionConstMatrixRef<double>>,
                                                            std::span<const Eigen::Index>,
                                                            Eigen::Index,
                                                            svdDeletionWorkspace<double> &,
                                                            svdDeletionBackend );

/** @} */

} // namespace math
} // namespace mx

#endif // math_svdDowndate_hpp
