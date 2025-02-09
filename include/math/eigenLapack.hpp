/** \file eigenLapack.hpp
 * \brief Interfaces to Lapack and BLAS for Eigen-like arrays.
 *
 * \author Jared R. Males (jaredmales@gmail.com)
 *
 * \ingroup gen_math_files
 *
 */

//***********************************************************************//
// Copyright 2015, 2016, 2017 Jared R. Males (jaredmales@gmail.com)
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

#ifndef math_eigenLapack_hpp
#define math_eigenLapack_hpp

#pragma GCC system_header
#include <Eigen/Dense>

#include <cmath>

#include "templateBLAS.hpp"
#include "templateLapack.hpp"

#include "../sys/timeUtils.hpp"

// #include "vectorUtils.hpp"

#include "../math/plot/gnuPlot.hpp"

namespace mx
{
namespace math
{

/// Calculates the lower triangular part of the covariance matrix of ims.
/** Uses cblas_ssyrk.  cv is resized to ims.cols() X ims.cols().
 * Calculates \f$ cv = A^T*A \f$.
 *
 *
 * \tparam eigenT1 is the eigen matrix/array type of cv.
 * \tparam eigenT2 is the eigen matrix/array type of ims
 *
 * \ingroup eigen_lapack
 */
template <typename eigenT1, typename eigenT2>
void eigenSYRK( eigenT1 &cv,       ///< [out] is the eigen matrix/array where to store the result
                const eigenT2 &ims ///< [in] is the eigen matrix/array (images as columns) to
                                   ///< calculate the covariance of
)
{
    cv.resize( ims.cols(), ims.cols() );

    math::syrk<typename eigenT1::Scalar>( /*const enum CBLAS_ORDER Order*/ CblasColMajor,
                                          /*const enum CBLAS_UPLO Uplo*/ CblasLower,
                                          /*const enum CBLAS_TRANSPOSE Trans*/ CblasTrans,
                                          /*const MXLAPACK_INT N*/ ims.cols(),
                                          /*const MXLAPACK_INT K*/ ims.rows(),
                                          /*const float alpha*/ 1.0,
                                          /*const float *A*/ ims.data(),
                                          /*const MXLAPACK_INT lda*/ ims.rows(),
                                          /*const float beta*/ 0.,
                                          /*float *C*/ cv.data(),
                                          /*const MXLAPACK_INT ldc*/ cv.rows() );
}

/// A struct to hold the working memory for eigenSYEVR and maintain it between calls if desired.
/** \todo this should have the working memory for the first exploratory call to ?syevr as well.
 */
template <typename floatT>
struct syevrMem
{
    char UPLO{ 'L' };
    MXLAPACK_INT n{ 0 };
    char RANGE{ 'A' };
    MXLAPACK_INT numeig{ 0 };
    MXLAPACK_INT IL{ 0 };
    MXLAPACK_INT IU{ 0 };

    MXLAPACK_INT sizeISuppZ{ 0 };
    MXLAPACK_INT sizeMinWork{ 0 };
    MXLAPACK_INT sizeWork{ 0 };
    MXLAPACK_INT sizeMinIWork{ 0 };
    MXLAPACK_INT sizeIWork{ 0 };

    MXLAPACK_INT *iSuppZ{ nullptr };

    floatT *minWork{ nullptr };
    floatT *work{ nullptr };

    MXLAPACK_INT *minIWork{ nullptr };
    MXLAPACK_INT *iWork{ nullptr };

    // Working memory for the higher level ev calc routines
    Eigen::Array<floatT, Eigen::Dynamic, Eigen::Dynamic> cvd, evecsd, evalsd;

    syevrMem()
    {
    }

    ~syevrMem()
    {
        if( iSuppZ )
        {
            ::free( iSuppZ );
        }

        if( minWork )
        {
            ::free( minWork );
        }

        if( work )
        {
            ::free( work );
        }

        if( minIWork )
        {
            ::free( minIWork );
        }

        if( iWork )
        {
            ::free( iWork );
        }
    }

    void free()
    {
        if( iSuppZ )
        {
            ::free( iSuppZ );
        }
        sizeISuppZ = 0;

        if( minWork )
        {
            ::free( minWork );
        }
        sizeMinWork = 0;

        if( work )
        {
            ::free( work );
        }
        sizeWork = 0;

        if( minIWork )
        {
            ::free( minIWork );
        }
        sizeMinIWork = 0;

        if( iWork )
        {
            ::free( iWork );
        }
        sizeIWork = 0;
    }
};

/// Calculate select eigenvalues and eigenvectors of an Eigen Array
/** Uses the templateLapack wrapper for syevr.
 *
 * \tparam arrT is the eigen-like type containing the data
 *
 * \returns -1000 on an malloc allocation error.
 * \returns the return code from syevr (info) otherwise.
 *
 * \ingroup eigen_lapack
 */
template <typename arrT>
MXLAPACK_INT eigenSYEVR( arrT &eigvec,    /**< [out] will contain the eigenvectors as columns*/
                         arrT &eigval,    /**< [out] will contain the eigenvalues*/
                         arrT &X,         /**< [in]  is a square matrix which is either upper
                                                     or lower (default) triangular*/
                         int ev0 = 0,     /**< [in] [opt] is the first desired eigenvalue
                                                          (default = 0)*/
                         int ev1 = -1,    /**< [in] [opt] if >= ev0 then this is the last desired
                                                          eigenvalue. If -1 all eigenvalues are
                                                          returned.*/
                         char UPLO = 'L', /**< [in] [opt] specifies whether X is upper ('U')
                                                          or lower ('L') triangular.
                                                          Default is ('L').*/
                         syevrMem<typename arrT::Scalar> *mem = 0 /**< [in] [opt] holds the working memory arrays,
                                                                                  can be re-passed to avoid unnecessary
                                                                                  re-allocations*/
)
{
    typedef typename arrT::Scalar calcT;

    MXLAPACK_INT numeig, info;
    char RANGE = 'A';

    MXLAPACK_INT localMem = 0;

    if( mem == 0 )
    {
        mem = new syevrMem<calcT>;
        localMem = 1;
    }

    MXLAPACK_INT n = X.rows();

    MXLAPACK_INT IL = 1;
    MXLAPACK_INT IU = n;
    if( ev0 >= 0 && ev1 >= ev0 )
    {
        RANGE = 'I';
        IL = ev0 + 1; // This is FORTRAN, after all
        IU = ev1;
    }

    eigvec.resize( n, IU - IL + 1 );
    eigval.resize( n, 1 );

    if( UPLO != mem->UPLO || n != mem->n || RANGE != mem->RANGE || mem->numeig != numeig || mem->IL != IL ||
        mem->IU != IU )
    {
        if( mem->sizeISuppZ < 2 * n )
        {
            if( mem->iSuppZ )
            {
                free( mem->iSuppZ );
            }

            mem->sizeISuppZ = 2 * n;
            mem->iSuppZ = (MXLAPACK_INT *)malloc( mem->sizeISuppZ * sizeof( MXLAPACK_INT ) );

            if( mem->iSuppZ == NULL )
            {
                mxError( "eigenSYEVR", MXE_ALLOCERR, "malloc failed in eigenSYEVR." );
                if( localMem )
                {
                    delete mem;
                }
                return -1000;
            }
        }

        //  Allocate minimum allowed sizes for workspace
        if( mem->sizeMinWork < 26 * n )
        {
            if( mem->minWork )
            {
                free( mem->minWork );
            }
            mem->sizeMinWork = 26 * n;
            mem->minWork = (calcT *)malloc( mem->sizeMinWork * sizeof( calcT ) );
        }

        // MXLAPACK_INT sizeIWork = 10 * n;
        if( mem->sizeMinIWork < 10 * n )
        {
            if( mem->minIWork )
            {
                free( mem->minIWork );
            }
            mem->sizeMinIWork = 10 * n;
            mem->minIWork = (MXLAPACK_INT *)malloc( mem->sizeMinIWork * sizeof( MXLAPACK_INT ) );
        }

        //  Query for optimum sizes for workspace
        info = math::syevr<calcT>( 'V',
                                   RANGE,
                                   UPLO,
                                   n,
                                   X.data(),
                                   n,
                                   0,
                                   0,
                                   IL,
                                   IU,
                                   math::lamch<calcT>( 'S' ),
                                   &numeig,
                                   eigval.data(),
                                   eigvec.data(),
                                   n,
                                   mem->iSuppZ,
                                   mem->minWork,
                                   -1,
                                   mem->minIWork,
                                   -1 );

        if( info != 0 )
        {
            mxError( "eigenSYEVR", MXE_LAPACKERR, "error from SYEVR" );
            if( localMem )
            {
                delete mem;
            }

            return info;
        }

        // Now allocate optimum sizes
        /* -- tested increasing by x10, didn't improve performance at all
         */
        if( mem->sizeWork < ( (MXLAPACK_INT)mem->minWork[0] ) * ( 1 ) )
        {
            if( mem->work )
            {
                free( mem->work );
            }

            mem->sizeWork = ( (MXLAPACK_INT)mem->minWork[0] ) * 1;
            mem->work = (calcT *)malloc( ( mem->sizeWork ) * sizeof( calcT ) );
        }

        if( mem->sizeIWork < mem->minIWork[0] * 1 )
        {
            if( mem->iWork )
            {
                free( mem->iWork );
            }

            mem->sizeIWork = mem->minIWork[0] * 1;
            mem->iWork = (MXLAPACK_INT *)malloc( ( mem->sizeIWork ) * sizeof( MXLAPACK_INT ) );
        }

        if( ( mem->work == NULL ) || ( mem->iWork == NULL ) )
        {
            mxError( "eigenSYEVR", MXE_ALLOCERR, "malloc failed in eigenSYEVR." );
            if( localMem )
            {
                delete mem;
            }
            return -1000;
        }

        mem->UPLO = UPLO;
        mem->n = n;
        mem->RANGE = RANGE;
        mem->numeig = numeig;
        mem->IL = IL;
        mem->IU = IU;
    }

    // Now actually do the calculationg
    info = math::syevr<calcT>( 'V',
                               RANGE,
                               UPLO,
                               n,
                               X.data(),
                               n,
                               0,
                               0,
                               IL,
                               IU,
                               math::lamch<calcT>( 'S' ),
                               &numeig,
                               eigval.data(),
                               eigvec.data(),
                               n,
                               mem->iSuppZ,
                               mem->work,
                               mem->sizeWork,
                               mem->iWork,
                               mem->sizeIWork );

    /*  Cleanup and exit  */

    if( localMem )
    {
        delete mem;
    }

    return info;
}

/// Calculate the eigenvectors and eigenvalues given a triangular matrix
/** Eigen-decomposition of the matrix is performed using \ref eigenSYEVR().
 *
 * \tparam evCalcT is the type in which to perform eigen-decomposition, which may be different
 *         from the input array and output arrays (which must be the same).
 * \tparam eigenT is a 2D Eigen-like type
 *
 * \ingroup eigen_lapack
 */
template <typename _evCalcT = double, typename eigenT>
MXLAPACK_INT calcEigenVecs( eigenT &evecs,               /**< [out] on exit contains the eigen vectors*/
                            eigenT &evals,               /**< [out] on exit contains the eigen vectors*/
                            eigenT &cv,                  /**< [in] a lower-triangle (in the Lapack sense) square
                                                                   covariance matrix.*/
                            int nVecs = 0,               /**< [in] [opt] The maximum number of modes to solve
                                                                              for. If 0 all modes are solved for.*/
                            bool normalize = false,      /**< [in] [opt] flag specifying whether or not to
                                                                              normalize the eigenvectors.*/
                            bool check = false,          /**< [in] [opt] flag specifying whether or not to
                                                                              check the eigenvalues/vectors for
                                                                              validity. Requires normalize=true.*/
                            syevrMem<_evCalcT> *mem = 0, /**< [in] [opt] A memory structure which can be
                                                                              re-used by SYEVR for efficiency.*/
                            double *t_eigenv = nullptr   /**< [out] [opt] if not null, will be filled in
                                                                               with the time taken to calculate
                                                                               eigenvalues.*/ )
{
    typedef _evCalcT evCalcT;
    typedef typename eigenT::Scalar realT;

    bool localMem = false;

    if( mem == 0 )
    {
        mem = new syevrMem<evCalcT>;
        localMem = true;
    }

    if( cv.rows() != cv.cols() )
    {
        std::cerr << "calcEigenVecs: Non-square covariance matrix input to calcKLModes\n";
        if( localMem )
        {
            delete mem;
        }
        return -1;
    }

    MXLAPACK_INT tNims = cv.rows();

    if( nVecs <= 0 || nVecs > tNims )
    {
        nVecs = tNims;
    }

    if( t_eigenv )
    {
        *t_eigenv = sys::get_curr_time();
    }

    mem->cvd = cv.template cast<evCalcT>();

    // Calculate eigenvectors and eigenvalues
    /* SYEVR sorts eigenvalues in ascending order, so we specifiy the top n_modes
     */
    MXLAPACK_INT info = eigenSYEVR( mem->evecsd, mem->evalsd, mem->cvd, tNims - nVecs, tNims, 'L', mem );

    if( info != 0 )
    {
        std::cerr << "calcEigenVecs: eigenSYEVR returned an error (info = " << info << ")\n";
        if( localMem )
        {
            delete mem;
        }

        return -1;
    }

    evecs = mem->evecsd.template cast<realT>();
    evals = mem->evalsd.template cast<realT>();

    if( normalize )
    {
        // Normalize the eigenvectors
        if( !check )
        {
            for( MXLAPACK_INT i = 0; i < nVecs; ++i )
            {
                evecs.col( i ) = evecs.col( i ) / sqrt( evals( i ) );
            }
        }
        else // here we check for invalid results and 0 things out
        {
            for( MXLAPACK_INT i = 0; i < nVecs; ++i )
            {
                if( evals( i ) == 0 )
                {
                    std::cerr << "got 0 eigenvalue (# " << i << ")\n";
                    evecs.col( i ) *= 0;
                }
                else if( evals( i ) < 0 )
                {
                    std::cerr << "got < 0 eigenvalue (# " << i << ")\n";
                    evecs.col( i ) *= 0;
                }
                else if( !std::isnormal( evals( i ) ) )
                {
                    std::cerr << "got not-normal eigenvalue (# " << i << ")\n";
                    evecs.col( i ) *= 0;
                }
                else
                {
                    evecs.col( i ) = evecs.col( i ) / sqrt( evals( i ) );
                }

                for( int r = 0; r < evecs.rows(); ++r )
                {
                    if( !std::isnormal( evecs.col( i )( r ) ) )
                    {
                        std::cerr << "got not-normal eigenvector entry (# " << i << "," << r << ")\n";
                        evecs.col( i ) *= 0;
                        continue;
                    }
                }
            }
        }
    }

    if( t_eigenv )
    {
        *t_eigenv = sys::get_curr_time() - *t_eigenv;
    }

    if( localMem )
    {
        delete mem;
    }

    return 0;

} // calcKLModes

/// Calculate the K-L modes, or principle components, given a covariance matrix.
/** Eigen-decomposition of the covariance matrix is performed using \ref eigenSYEVR().
 *
 * \tparam evCalcT is the type in which to perform eigen-decomposition.
 * \tparam eigenT is a 2D Eigen-like type
 * \tparam eigenT1 is a 2D Eigen-like type.
 *
 * \ingroup eigen_lapack
 */
template <typename _evCalcT = double, typename eigenT, typename eigenT1>
MXLAPACK_INT calcKLModes( eigenT &klModes,             /**< [out] on exit contains the K-L modes (or P.C.s) */
                          eigenT &cv,                  /**< [in] a lower-triangle (in the Lapack sense) square
                                                                 covariance matrix.*/
                          const eigenT1 &Rims,         /**< [in] The reference data.  cv.rows() == Rims.cols().*/
                          int n_modes = 0,             /**< [in] [opt] Tbe maximum number of modes to solve for.
                                                                       If 0 all modes are solved for.*/
                          syevrMem<_evCalcT> *mem = 0, /**< [in] [opt] A memory structure which can be re-used
                                                                       by SYEVR for efficiency.*/
                          double *t_eigenv = nullptr,  /**< [out] [opt] if not null, will be filled in with the time
                                                                        taken to calculate eigenvalues.*/
                          double *t_klim = nullptr     /**< [out] [opt] if not null, will be filled in with the time
                                                                        taken to calculate the KL modes.*/)
{
    typedef _evCalcT evCalcT;
    typedef typename eigenT::Scalar realT;

    bool localMem = false;

    if( mem == 0 )
    {
        mem = new syevrMem<evCalcT>;
        localMem = true;
    }

    if( cv.rows() != Rims.cols() )
    {
        std::cerr << "Covariance matrix - reference image size mismatch in calcKLModes\n";
        return -1;
    }

    eigenT evecs, evals;

    MXLAPACK_INT tNims = cv.rows();
    MXLAPACK_INT tNpix = Rims.rows();

    if( n_modes <= 0 || n_modes > tNims )
    {
        n_modes = tNims;
    }

    // Eigen::Array<evCalcT, Eigen::Dynamic, Eigen::Dynamic> evecsd, evalsd;
    mem->cvd = cv.template cast<evCalcT>();
    MXLAPACK_INT info = calcEigenVecs( mem->evecsd, mem->evalsd, mem->cvd, n_modes, true, true, mem, t_eigenv );

    if( info != 0 )
    {
        std::cerr << "calckKLModes: eigenSYEVR returned an error (info = " << info << ")\n";
        return -1;
    }

    evecs = mem->evecsd.template cast<realT>();
    evals = mem->evalsd.template cast<realT>();

    klModes.resize( n_modes, tNpix );

    if( t_klim )
    {
        *t_klim = sys::get_curr_time();
    }

    // Now calculate KL images
    /*
     *  KL = E^T * R  ==> C = A^T * B
     */
    gemm<realT>( CblasColMajor,
                 CblasTrans,
                 CblasTrans,
                 n_modes,
                 tNpix,
                 tNims,
                 1.,
                 evecs.data(),
                 cv.rows(),
                 Rims.data(),
                 Rims.rows(),
                 0.,
                 klModes.data(),
                 klModes.rows() );

    if( t_klim )
    {
        *t_klim = sys::get_curr_time() - *t_klim;
    }

    if( localMem )
    {
        delete mem;
    }

    return 0;

} // calcKLModes

/// Compute the SVD of an Eigen::Array using LAPACK's xgesdd
/** Computes the SVD of A, \f$ A = U S V^T \f$.
 *
 * \returns 0 on success
 * \returns -i on error in ith parameter (from LAPACK xgesdd)
 * \returns >0 did not converge (from LAPACK xgesdd)
 *
 * \tparam dataT is either float or double.
 *
 * \ingroup eigen_lapack
 */
template <typename dataT>
MXLAPACK_INT eigenGESDD( Eigen::Array<dataT, -1, -1> &U,  ///< [out] the A.rows() x A.rows() left matrix
                         Eigen::Array<dataT, -1, -1> &S,  ///< [out] the A.cols() x 1 matrix of singular values
                         Eigen::Array<dataT, -1, -1> &VT, /**< [out] the A.cols() x A.cols() right matrix, note this
                                                                     is the transpose. */
                         Eigen::Array<dataT, -1, -1> &A   ///< [in]  the input matrix to be decomposed
)
{
    char JOBZ = 'A';
    MXLAPACK_INT M = A.rows();
    MXLAPACK_INT N = A.cols();
    MXLAPACK_INT LDA = M;
    S.resize( N, 1 );
    U.resize( M, M );
    MXLAPACK_INT LDU = M;
    VT.resize( N, N );
    MXLAPACK_INT LDVT = N;

    dataT wkOpt;
    MXLAPACK_INT LWORK = -1;

    MXLAPACK_INT *IWORK = new MXLAPACK_INT[8 * M];
    MXLAPACK_INT INFO;

    math::gesdd<
        dataT>( JOBZ, M, N, A.data(), LDA, S.data(), U.data(), LDU, VT.data(), LDVT, &wkOpt, LWORK, IWORK, INFO );

    LWORK = wkOpt;
    // delete WORK;
    dataT *WORK = new dataT[LWORK];

    INFO = math::gesdd<
        dataT>( JOBZ, M, N, A.data(), LDA, S.data(), U.data(), LDU, VT.data(), LDVT, WORK, LWORK, IWORK, INFO );

    delete[] WORK;
    delete[] IWORK;

    return INFO;
}

#define MX_PINV_NO_INTERACT 0
#define MX_PINV_PLOT 1
#define MX_PINV_ASK 2
#define MX_PINV_ASK_NMODES 4

/// Calculate the pseudo-inverse of a patrix given its SVD
/** Given the SVD of A, \f$ A = U S V^T \f$, as calculated by eigenGESDD the psuedo-inverse is
 * calculated as \f$ A^+ = V S^+ U^T\f$.
 *
 * The parameter \p interact is intepreted as a bitmask.  The values can be
 * - \ref MX_PINV_PLOT which will cause a plot to be displayed of the singular values
 * - \ref MX_PINV_ASK which will ask the user for a max. condition number using stdin
 * - \ref MX_PINV_ASK_NMODES which will ask the user for a max number of modes to include using
 * stdin.  Overrides MX_PINV_ASK. If \p interact is 0 then no interaction is used and \p
 * maxCondition controls the inversion.
 *
 * \tparam dataT is either float or double.
 *
 * \ingroup eigen_lapack
 */
template <typename dataT>
int eigenPseudoInverse( Eigen::Array<dataT, -1, -1> &PInv, ///< [out] The pseudo-inverse of A
                        dataT &condition,                  ///< [out] The final condition number.
                        int &nRejected,                    ///< [out] The number of eigenvectors rejected
                        Eigen::Array<dataT, -1, -1> &U,    ///< [in] the A.rows() x A.rows() left matrix
                        Eigen::Array<dataT, -1, -1> &S,    ///< [in] the A.cols() x 1 matrix of singular values
                        Eigen::Array<dataT, -1, -1> &VT,   /**< [in] the A.cols() x A.cols() right matrix, note this
                                                                      is the transpose. */
                        int minMN,                         ///< [in]  The minimum size of the matrix to invert.
                        dataT &maxCondition,               /**< [in]  If \> 0, the maximum condition number desired.
                                                                      If \<0 the number of modes to keep. Used to
                                                                      threshold the singular values.  Set to 0 to
                                                                      include all eigenvalues/vectors. Ignored if
                                                                      interactive. */
                        dataT alpha = 0,                   /**< [in] [opt] the Tikhonov regularization value, as
                                                                     a fraction of the highest singular value.
                                                                     If alpha < 0, then it is treated as a (positive)
                                                                     floor (as a fraction of highest singular value)
                                                                     for the singular values, which is not the same as
                                                                    Tikhonov (alpha > 0).*/
                        int interact = MX_PINV_NO_INTERACT /**< [in] [opt] a bitmask controlling interaction.
                                                                     See above.*/
)
{

    dataT Smax = S.maxCoeff();

    if( alpha > 0 )
    {
        for( MXLAPACK_INT i = 0; i < S.rows(); ++i )
        {
            S( i ) = ( pow( S( i ), 2 ) + pow( alpha * Smax, 2 ) ) / S( i );
        }
    }

    if( alpha < 0 )
    {
        for( MXLAPACK_INT i = 0; i < S.rows(); ++i )
        {
            S( i ) = S( i ) + -alpha * Smax;
        }
    }

    int modesToReject = 0;
    if( maxCondition < 0 ) // Rejecting mode numbers
    {
        modesToReject = -maxCondition;

        if( modesToReject - 1 < S.rows() )
        {
            maxCondition = Smax / S( modesToReject - 1, 0 );
        }
    }

    if( interact & MX_PINV_PLOT )
    {
        gnuPlot gp;
        gp.command( "set title \"SVD Singular Values\"" );
        gp.logy();
        gp.plot( S.data(), S.rows(), " w lp", "singular values" );
    }

    if( interact & MX_PINV_ASK && !( interact & MX_PINV_ASK_NMODES ) )
    {
        dataT mine;
        std::cout << "Maximum singular value: " << Smax << "\n";
        std::cout << "Minimum singular value: " << S.minCoeff() << "\n";
        std::cout << "Enter singular value threshold: ";
        std::cin >> mine;

        if( mine > 0 )
        {
            maxCondition = Smax / mine;
        }
        else
        {
            maxCondition = Smax / S( S.rows() - 1, 0 );
        }
    }
    else if( interact & MX_PINV_ASK_NMODES )
    {
        unsigned mine;
        std::cout << "Maximum singular value: " << Smax << "\n";
        std::cout << "Minimum singular value: " << S.minCoeff() << "\n";
        std::cout << "Enter number of modes to keep: ";
        std::cin >> mine;
        modesToReject = S.rows() - mine;

        if( modesToReject <= 0 || modesToReject > S.rows() )
        {
            modesToReject = 0;
        }

        maxCondition = -modesToReject;
    }

    Eigen::Array<dataT, -1, -1> sigma;
    sigma.resize( S.rows(), S.rows() );
    sigma.setZero();

    nRejected = 0;

    if( maxCondition > 0 )
    {
        dataT threshold = 0;

        if( maxCondition > 0 )
        {
            threshold = Smax / maxCondition;

            condition = 1;

            for( MXLAPACK_INT i = 0; i < S.rows(); ++i )
            {
                if( S( i ) >= threshold )
                {
                    sigma( i, i ) = 1. / S( i );
                    if( Smax / S( i ) > condition )
                    {
                        condition = Smax / S( i );
                    }
                }
                else
                {
                    sigma( i, i ) = 0;
                    ++nRejected;
                }
            }
        }
    }
    else // rejecting modes
    {
        std::cerr << "rejecting based on modes\n";
        std::cerr << " modes to reject: " << modesToReject << "\n";
        for( MXLAPACK_INT i = 0; i < S.rows(); ++i )
        {
            if( i < S.rows() - modesToReject )
            {
                sigma( i, i ) = 1. / S( i );
                if( Smax / S( i ) > condition )
                    condition = Smax / S( i );
            }
            else
            {
                sigma( i, i ) = 0;
                ++nRejected;
            }
        }
    }

    if( interact & MX_PINV_PLOT )
    {
        std::vector<dataT> vsig( sigma.rows() );
        for( int rr = 0; rr < sigma.rows(); ++rr )
        {
            vsig[rr] = sigma( rr, rr );
        }

        gnuPlot gp;
        gp.command( "set title \"Inverted Singular Values\"" );
        gp.logy();
        gp.plot( vsig.data(), vsig.size(), " w lp", "inverted singular values" );
    }

    if( interact & MX_PINV_ASK || interact & MX_PINV_ASK_NMODES )
    {
        dataT mine;
        std::cout << "Modes Rejected: " << nRejected << "\n";
        std::cout << "Condition Number: " << condition << "\n";
    }

    PInv = ( VT.matrix().transpose() * sigma.matrix().transpose() ) *
           U.block( 0, 0, U.rows(), minMN ).matrix().transpose();

    return 0;
}

/// Calculate the pseudo-inverse of a patrix using the SVD
/** First computes the SVD of A, \f$ A = U S V^T \f$, using eigenGESDD.  Then the psuedo-inverse is
 * calculated as \f$ A^+ = V S^+ U^T\f$.
 *
 * The parameter \p interact is intepreted as a bitmask.  The values can be
 * - \ref MX_PINV_PLOT which will cause a plot to be displayed of the singular values
 * - \ref MX_PINV_ASK which will ask the user for a max. condition number using stdin
 * - \ref MX_PINV_ASK_NMODES which will ask the user for a max number of modes to include using
 * stdin.  Overrides MX_PINV_ASK. If \p interact is 0 then no interaction is used and \p
 * maxCondition controls the inversion.
 *
 * \tparam dataT is either float or double.
 *
 * \ingroup eigen_lapack
 */
template <typename dataT>
int eigenPseudoInverse( Eigen::Array<dataT, -1, -1> &PInv, ///< [out] The pseudo-inverse of A
                        dataT &condition,                  ///< [out] The final condition number.
                        int &nRejected,                    ///< [out] The number of eigenvectors rejected
                        Eigen::Array<dataT, -1, -1> &U,    ///< [out] the A.rows() x A.rows() left matrix
                        Eigen::Array<dataT, -1, -1> &S,    ///< [out] the A.cols() x 1 matrix of singular values
                        Eigen::Array<dataT, -1, -1> &VT,   /**< [out] the A.cols() x A.cols() right matrix, note this
                                                                      is the transpose. */
                        Eigen::Array<dataT, -1, -1> &A,    ///< [in]  The matrix to invert.  This will be modified!
                        dataT &maxCondition,               /**< [in]  If \> 0, the maximum condition number desired.
                                                                      If \<0 the number of modes to keep. Used to
                                                                      threshold the singular values.  Set to 0 to
                                                                      include all eigenvalues/vectors. Ignored if
                                                                      interactive. */
                        dataT alpha = 0,                   /**< [in] [opt] the Tikhonov regularization value, as
                                                                     a fraction of the highest singular value.
                                                                     If alpha < 0, then it is treated as a (positive)
                                                                     floor (as a fraction of highest singular value)
                                                                     for the singular values, which is not the same as
                                                                    Tikhonov (alpha > 0).*/
                        int interact = MX_PINV_NO_INTERACT /**< [in] [opt] a bitmask controlling interaction.
                                                                     See above.*/
)
{

    int minMN = std::min( A.rows(), A.cols() );

    MXLAPACK_INT info;
    info = eigenGESDD( U, S, VT, A );

    if( info != 0 )
    {
        std::cerr << "eigenPseudoInverse: eigenGESDD failed with info = " << info << "\n";
        return info;
    }

    return eigenPseudoInverse( PInv, condition, nRejected, U, S, VT, minMN, maxCondition, alpha, interact );
}

/// Calculate the pseudo-inverse of a matrix using the SVD
/** First computes the SVD of A, \f$ A = U S V^T \f$, using eigenGESDD.  Then the psuedo-inverse is
 * calculated as \f$ A^+ = V S^+ U^T\f$. This interface does not provide access to U, S and VT.
 *
 * The parameter \p interact is intepreted as a bitmask.  The values can be
 * - \ref MX_PINV_PLOT which will cause a plot to be displayed of the singular values
 * - \ref MX_PINV_ASK which will ask the user for a max. condition number using stdin
 * - \ref MX_PINV_ASK_NMODES which will ask the user for a max number of modes to include using
 * stdin.  Overrides MX_PINV_ASK. If \p interact is 0 then no interaction is used and maxCondition
 * controls the inversion.
 *   *
 * \tparam dataT is either float or double.
 *
 * \overload
 *
 * \ingroup eigen_lapack
 */
template <typename dataT>
int eigenPseudoInverse( Eigen::Array<dataT, -1, -1> &PInv, ///< [out] The pseudo-inverse of A
                        dataT &condition,                  ///< [out] The final condition number.
                        int &nRejected,                    /**< [out] The number of eigenvectors
                                                                      rejected*/
                        Eigen::Array<dataT, -1, -1> &A,    /**< [in]  The matrix to invert, will be
                                                                       altered!*/
                        dataT &maxCondition,               /**< [in]  If \> 0, the maximum condition number desired.
                                                                      If \<0 the number of modes to keep. Used to
                                                                      threshold the singular values. Set to 0 to
                                                                      include all eigenvalues/vectors. Ignored if
                                                                      interactive.*/
                        dataT alpha = 0,                   /**< [in] [opt] the Tikhonov regularization value,
                                                                     as a fraction of the highest singular value. If
                                                                     alpha < 0, then it is treated as a (positive)
                                                                     floor (as a fraction of highest singular value)
                                                                     for the singular values, which is not the same
                                                                     as Tikhonov (alpha > 0).*/
                        int interact = MX_PINV_NO_INTERACT /**< [in] [opt] a bitmask controlling interaction.
                                                                     See above.*/
)
{
    Eigen::Array<dataT, -1, -1> S, U, VT;

    return eigenPseudoInverse( PInv, condition, nRejected, U, S, VT, A, maxCondition, alpha, interact );
}

} // namespace math
} // namespace mx

#endif // math_eigenLapack_hpp
