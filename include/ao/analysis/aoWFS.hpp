/** \file aoWFS.hpp
 * \author Jared R. Males (jaredmales@gmail.com)
 * \brief Definitions of various analytic wavefront sensors
 * \ingroup mxAO_files
 *
 */

#ifndef aoWFS_hpp
#define aoWFS_hpp

#include "../../math/constants.hpp"
#include "../../mxError.hpp"
#include "../../mxException.hpp"
#include "../../improc/eigenImage.hpp"
#include "../../ioutils/fits/fitsFile.hpp"

namespace mx
{
namespace AO
{
namespace analysis
{

/// The ideal wavefront sensor sensitivity function.
/** Provides the \f$ \beta_p \f$ parameter of Guyon, 2005 \cite guyon_2005
 * for the ideal WFS.
 *
 * This is the base class for all WFS.
 *
 * \tparam realT is the floating point type used for calculations
 * \tparam iosT is an output stream type with operator \<\< defined (default is std::ostream)
 *
 * \ingroup mxAOAnalytic
 */
template <typename realT, typename iosT = std::ostream>
struct wfs
{
    std::string _id;

    /// Constructor
    /** Only sets the value of _id.
     */
    wfs()
    {
        _id = "Ideal WFS";
    }

    /// Destructor
    /** Declared virtual so more complicated derived types can be created.
     */
    virtual ~wfs()
    {
        return;
    }

    /// Get the photon noise sensitivity at a spatial frequency.
    /** The sensitivity of the ideal WFS is 1 at all k \cite guyon_2005.
     *
     * \returns the sensitivity to photon noise parameter
     */
    virtual realT beta_p( int m,   ///< [in] the spatial frequency index for u (not used by this WFS)
                          int n,   ///< [in] the spatial frequency index for v (not used by this WFS)
                          realT D, ///< [in] the telescope diameter (not used by this WFS)
                          realT d, ///< [in] the sub-ap spacing (not used by this WFS)
                          realT r0 ///< [in] Fried's parameter (not used by this WFS)
    )
    {
        // Stuff a sock in the compiler's mouth:
        static_cast<void>( m );
        static_cast<void>( n );
        static_cast<void>( D );
        static_cast<void>( d );
        static_cast<void>( r0 );

        return static_cast<realT>( 1 );
    }

    /// Get the read noise sensitivity at a spatial frequency.
    /** Here we assume beta_r is the same as beta_p.
     *
     * \returns the sensitivity to read noise parameter
     */
    virtual realT beta_r( int m,   ///< [in] the spatial frequency index for u (not used by this WFS)
                          int n,   ///< [in] the spatial frequency index for v (not used by this WFS)
                          realT D, ///< [in] the telescope diameter (not used by this WFS)
                          realT d, ///< [in] the sub-ap spacing (not used by this WFS)
                          realT r0 ///< [in] Fried's parameter (not used by this WFS)
    )
    {
        return beta_p( m, n, D, d, r0 );
    }

    /// Dump the details of the WFS to an io stream.
    /** Is virtual so that derived types can add parameters.
     */
    virtual iosT &dumpWFS( iosT &ios )
    {
        ios << "# WFS Parameters:\n";
        ios << "#    ID = " << _id << '\n';

        return ios;
    }
};

/// The unmodulated pyramid wavefront sensor sensitivity function.
/** Provides the \f$ \beta_p \f$ parameter of Guyon, 2005 \cite guyon_2005
 * for the unmodulated PyWFS.
 *
 * \tparam realT is the floating point type used for calculations
 * \tparam iosT is an output stream type with operator \<\< defined (default is std::ostream)
 *
 * \ingroup mxAOAnalytic
 */
template <typename realT, typename iosT = std::ostream>
struct pywfsUnmod : public wfs<realT, iosT>
{
    pywfsUnmod()
    {
        this->_id = "Unmodulated Pyramid";
    }

    /// Get the photon noise sensitivity at a spatial frequency.
    /** The sensitivity of the unmodulated PyWFS is \f$ \sqrt{2} \f$ at all k.
     *
     * \returns the sensitivity to photon noise parameter
     */
    virtual realT beta_p( int m,   ///< [in] the spatial frequency index for u (not used by this WFS)
                          int n,   ///< [in] the spatial frequency index for v (not used by this WFS)
                          realT D, ///< [in] the telescope diameter (not used by this WFS)
                          realT d, ///< [in] the sub-ap spacing (not used by this WFS)
                          realT r0 ///< [in] Fried's parameter (not used by this WFS)
    )
    {
        // Stuff a sock in the compiler's mouth:
        static_cast<void>( m );
        static_cast<void>( n );
        static_cast<void>( D );
        static_cast<void>( d );
        static_cast<void>( r0 );

        return math::root_two<realT>();
    }

    /// Get the read noise sensitivity at a spatial frequency.
    /** Here we assume that beta_r is the same as beta_p.
     *
     * \returns the sensitivity to read noise parameter
     */
    virtual realT beta_r( int m,   ///< [in] the spatial frequency index for u (not used by this WFS)
                          int n,   ///< [in] the spatial frequency index for v (not used by this WFS)
                          realT D, ///< [in] the telescope diameter (not used by this WFS)
                          realT d, ///< [in] the sub-ap spacing (not used by this WFS)
                          realT r0 ///< [in] Fried's parameter (not used by this WFS)
    )
    {
        return beta_p( m, n, D, d, r0 );
    }
};

/// The asymptotic modulated pyramid wavefront sensor sensitivity function.
/** Provides the \f$ \beta_p \f$ parameter of Guyon, 2005 \cite guyon_2005
 * for the modulated PyWFS in the asymptotic limit.
 *
 * \tparam realT is the floating point type used for calculations
 * \tparam iosT is an output stream type with operator \<\< defined (default is std::ostream)
 *
 * \ingroup mxAOAnalytic
 */
template <typename realT, typename iosT = std::ostream>
struct pywfsModAsymptotic : public wfs<realT, iosT>
{
    pywfsModAsymptotic()
    {
        this->_id = "Asymptotic Modulated Pyramid";
    }

    /// Get the photon sensitivity at a spatial frequency.
    /** The photon noise sensitivity of the asymptotic modulated PyWFS is \f$ 2 \sqrt{2} \f$ at all k.
     *
     * \returns the sensitivity to photon noise parameter
     */
    virtual realT beta_p( int m,   ///< [in] the spatial frequency index for u (not used by this WFS)
                          int n,   ///< [in] the spatial frequency index for v (not used by this WFS)
                          realT D, ///< [in] the telescope diameter (not used by this WFS)
                          realT d, ///< [in] the sub-ap spacing (not used by this WFS)
                          realT r0 ///< [in] Fried's parameter (not used by this WFS)
    )
    {
        // Stuff a sock in the compiler's mouth:
        static_cast<void>( m );
        static_cast<void>( n );
        static_cast<void>( D );
        static_cast<void>( d );
        static_cast<void>( r0 );

        return static_cast<realT>( 2 ) * math::root_two<realT>();
    }

    /// Get the read noise sensitivity at a spatial frequency.
    /** Here we assume beta_r is the same as beta_p
     *
     * \returns the sensitivity to read noise parameter
     */
    virtual realT beta_r( int m,   ///< [in] the spatial frequency index for u (not used by this WFS)
                          int n,   ///< [in] the spatial frequency index for v (not used by this WFS)
                          realT D, ///< [in] the telescope diameter (not used by this WFS)
                          realT d, ///< [in] the sub-ap spacing (not used by this WFS)
                          realT r0 ///< [in] Fried's parameter (not used by this WFS)
    )
    {
        return beta_p( m, m, D, d, r0 );
    }
};

/// The shack hartmann wavefront sensor sensitivity function.
/** Provides the \f$ \beta_p \f$ parameter of Guyon, 2005 \cite guyon_2005
 * for the shack hartmann WFS.
 *
 * \tparam realT is the floating point type used for calculations
 * \tparam iosT is an output stream type with operator \<\< defined (default is std::ostream)
 *
 * \ingroup mxAOAnalytic
 */
template <typename realT, typename iosT = std::ostream>
struct shwfs : public wfs<realT, iosT>
{
    shwfs()
    {
        this->_id = "Shack Hartmann";
    }

    /// Get the photon noise sensitivity at a spatial frequency.
    /** The photon noise sensitivity of the shack hartmann WFS
     *
     * \returns the sensitivity to photon noise parameter
     */
    virtual realT beta_p( int m,   ///< [in] the spatial frequency index for u
                          int n,   ///< [in] the spatial frequency index for v
                          realT D, ///< [in] the telescope diameter
                          realT d, ///< [in] the sub-ap spacing
                          realT r0 ///< [in] Fried's parameter
    )
    {
        realT k = sqrt( m * m + n * n ) / D;

        return 1.48 / ( k * d ) * sqrt( 1 + pow( d / r0, 2 ) );
    }

    /// Get the read noise sensitivity at a spatial frequency.
    /** Here we assume beta_r = beta_p
     *
     * \returns the sensitivity to read noise parameter
     */
    virtual realT beta_r( int m,   ///< [in] the spatial frequency index for u
                          int n,   ///< [in] the spatial frequency index for v
                          realT D, ///< [in] the telescope diameter
                          realT d, ///< [in] the sub-ap spacing
                          realT r0 ///< [in] Fried's parameter
    )
    {
        return beta_p( m, n, D, d, r0 );
    }
};

/// The calculated WFS uses sensitivities provided by FITS files
/** Provides the \f$ \beta_p \f$ and \f$ \beta_r \f$ parameters
 * from FITS files.
 *
 * \tparam realT is the floating point type used for calculations
 * \tparam iosT is an output stream type with operator \<\< defined (default is std::ostream)
 *
 * \ingroup mxAOAnalytic
 */
template <typename realT, typename iosT = std::ostream>
struct calculatedWFS : public wfs<realT, iosT>
{
    std::string m_beta_p_file;
    std::string m_beta_r_file;
    bool m_sensitivity{ false };

    improc::eigenImage<realT> m_beta_p;
    improc::eigenImage<realT> m_beta_r;

    calculatedWFS()
    {
        this->_id = "Calculated WFS";
    }

    /// Get the photon noise sensitivity at a spatial frequency.
    /** The photon noise sensitivity from the FITS file is returned.
     *
     * \returns the sensitivity to photon noise parameter
     */
    virtual realT beta_p( int m,   ///< [in] the spatial frequency index for u
                          int n,   ///< [in] the spatial frequency index for v
                          realT D, ///< [in] the telescope diameter
                          realT d, ///< [in] the sub-ap spacing
                          realT r0 ///< [in] Fried's parameter
    )
    {

        if( m_beta_p.rows() == 0 )
        {
#pragma omp critical // Make sure we don't race in m/t
            {
                if( m_beta_p.rows() == 0 ) // Check again after critical locks
                {
                    fits::fitsFile<realT> ff;
                    ff.read( m_beta_p, m_beta_p_file );

                    if( m_sensitivity )
                    {
                        m_beta_p = 1.0 / m_beta_p;
                    }
                }
            }
        }

        int midx = 0.5 * ( m_beta_p.rows() - 1.0 ) + m;
        int nidx = 0.5 * ( m_beta_p.cols() - 1.0 ) + n;

        if( midx > m_beta_p.rows() - 1 || midx < 0 )
        {
#pragma omp critical
            {
                std::cerr << "calculatedWFS::beta_p: m index out of range. Got m = " << m
                          << "  / beta_p.rows = " << m_beta_p.rows() << "\n";
            }
            // Just return a huge number to make this spatial frequency unsenseable
            return std::numeric_limits<realT>::max();
        }

        if( nidx > m_beta_p.cols() - 1 || nidx < 0 )
        {
#pragma omp critical
            {
                std::cerr << "calculatedWFS::beta_p: n index out of range. Got n = " << n
                          << "  / beta_p.cols = " << m_beta_p.cols() << "\n";
            }
            // Just return a huge number to make this spatial frequency unsenseable
            return std::numeric_limits<realT>::max();
        }

        return m_beta_p( midx, nidx );
    }

    /// Get the read noise sensitivity at a spatial frequency.
    /** The read noise sensitivity from the FITS file is returned.
     *
     * \returns the sensitivity to read noise parameter
     */
    virtual realT beta_r( int m,   ///< [in] the spatial frequency index for u
                          int n,   ///< [in] the spatial frequency index for v
                          realT D, ///< [in] the telescope diameter
                          realT d, ///< [in] the sub-ap spacing
                          realT r0 ///< [in] Fried's parameter
    )
    {
        if( m_beta_r.rows() == 0 )
        {
            fits::fitsFile<realT> ff;
            ff.read( m_beta_r, m_beta_r_file );

            if( m_sensitivity )
                m_beta_r = 1.0 / m_beta_r;
        }

        int midx = 0.5 * ( m_beta_r.rows() - 1.0 ) + m;
        int nidx = 0.5 * ( m_beta_r.cols() - 1.0 ) + n;

        if( midx > m_beta_r.rows() - 1 || midx < 0 )
        {
            mxThrowException( err::sizeerr, "calculatedWFS::beta_r", "m index out of range" );
        }

        if( nidx > m_beta_r.cols() - 1 || nidx < 0 )
        {
            mxThrowException( err::sizeerr, "calculatedWFS::beta_r", "n index out of range" );
        }

        return m_beta_r( midx, nidx );
    }

    /// Dump the details of the WFS to an io stream.
    /** Is virtual so that derived types can add parameters.
     */
    virtual iosT &dumpWFS( iosT &ios )
    {
        wfs<realT, iosT>::dumpWFS( ios );

        ios << "#    beta_p = " << m_beta_p_file << '\n';
        ios << "#    beta_r = " << m_beta_r_file << '\n';
        ios << "#    sensitivity = " << std::boolalpha << m_sensitivity << "\n";

        return ios;
    }
};

} // namespace analysis
} // namespace AO
} // namespace mx

#endif // aoWFS_hpp
