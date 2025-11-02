/** \file pyramidSensor.hpp
 * \author Jared R. Males (jaredmales@gmail.com)
 * \brief Declaration and definition of a standard 4 quadrant pyramid WFS.
 * \ingroup mxAO_sim_files
 *
 */

#ifndef mx_AO_sim_pyramidSensor_hpp
#define mx_AO_sim_pyramidSensor_hpp

#include <omp.h>

#include "../../mxlib.hpp"

#include "../../wfp/imagingUtils.hpp"
#include "../../wfp/fraunhoferPropagator.hpp"
#include "../../sys/timeUtils.hpp"

#include "../../improc/eigenImage.hpp"
#include "../../improc/imageMasks.hpp"

#include "../../math/constants.hpp"
#include "../../math/geo.hpp"

#ifdef MXLIB_CUDA
#include "../../math/cuda/cublasHandle.hpp"
#endif

#include "wavefront.hpp"

#ifdef DEBUG
#define BREAD_CRUMB std::cout << "DEBUG: " << __FILE__ << " " << __LINE__ << "\n";
#else
#define BREAD_CRUMB
#endif

namespace mx
{
namespace AO
{
namespace sim
{

template <typename _realT>
struct wfsImageT
{
    typedef _realT realT;

    unsigned iterNo;

    /// The wavefront sensor detector image type
    typedef Eigen::Array<realT, Eigen::Dynamic, Eigen::Dynamic> imageT;

    imageT image;

    imageT tipImage;
};

/// A Pyramid Sensor Simulation
/**
 *
 * \tparam _realT is the real floating point type used for calculations
 * \tparam _detectorT is the detector used to record the PyWFS image.  Must conform to the mx::AO::sim detector
 * interface specifications.
 */
template <typename _realT, typename _detectorT, int _cudaGPU = 0>
class pyramidSensor
{
  public:
    static constexpr int cudaGPU = _cudaGPU;

    /// The real floating point type used for calculations
    typedef _realT realT;

    /// The complex floating point type used for calculations
    typedef std::complex<realT> complexT;

    /// The wavefront data type
    typedef wavefront<realT> wavefrontT;

    typedef improc::eigenImage<realT> realImageT;

    /// The wavefront complex field type
    typedef improc::eigenImage<std::complex<realT>> complexFieldT;

    /// The real array data type
    typedef typename wfp::fraunhoferPropagatorArrayT<realImageT, cudaGPU>::arrayT realArrayT;

    /// The complex array data type
    typedef typename wfp::fraunhoferPropagatorArrayT<complexFieldT, cudaGPU>::arrayT complexArrayT;

    /// The wavefront sensor detector image type
    typedef _detectorT detectorT;

  public:
    /// Default constructor
    pyramidSensor();

    /// Destructor
    ~pyramidSensor();

    /** \name Standard WFS Interface
     *
     * @{
     */

  protected:
    /* Standard WFS Interface: */
    uint32_t m_wfSz{ 0 };    ///< Size of the wavefront in pixels

    uint32_t m_detRows{ 0 }; ///< The number of rows of the WFS detector.  After forming the image the WFS detector
                             ///< plane is re-binned to this.

    uint32_t m_detCols{ 0 }; ///< The number of columns of the WFS detector.  After forming the image the WFS detector
                             ///< plane is re-binned to this.

    realT m_lambda{ 0 };     ///< Central wavelength, in meters

    /// \todo when the filter should be set with astrospectrum, and should re-calculate the central wavelength.
    /// \todo need to verify that the wavefront propagation is appropriately chromatic
    std::vector<realT> m_wavelengths;      ///< Vector of wavelengths in the WFS bandpass
    std::vector<realT> _wavelengthWeights; ///< The relative weights of the wavelengths

    int m_iTime{ 1 };                      ///< Integration time in loop steps

    int m_roTime{ 1 };                     ///< Readout time in loop steps

    realT m_simStep{ 0.001 };              ///< The simulation stepsize in seconds.

  public:
    /// Get the wavefront size in pixels
    /**
     * \returns the wavefront size in pixels
     */
    int wfSz();

    /// Set the wavefront size in pixels.
    /**
     */
    void wfSz( const uint32_t &sz /**< the new size*/ );

    /// Get the detector rows  in pixels
    /**
     * \returns m_detRows
     */
    uint32_t detRows();

    /// Get the detector columns  in pixels
    /**
     * \returns m_detCols
     */
    uint32_t detCols();

    /// Set the detector columns in pixels.
    /**
     */
    void detSize( const uint32_t &nrows, ///< The number of rows
                  const uint32_t &ncols  ///< The number of columns
    );

    /// Get the PyWFS central wavelength
    /**
     * \returns the central wavelength in meters
     */
    realT lambda();

    /// Set the PyWFS central wavelength
    /**
     */
    void lambda( const realT &l /**< The central wavelength, in meters*/ );

    /// Get the PyWFS integration time, in time steps
    int iTime();

    /// Set the PyWFS integration time, in time steps
    void iTime( const uint32_t &it /**<  the new integration time*/ );

    /// Get the PyWFS detector readout time, in time steps
    int roTime();

    /// Set the PyWFS detector readout time, in time steps
    void roTime( const uint32_t &rt /**< the new readout time*/ );

    /// Get the simulation step-size, in seconds.
    realT simStep();

    /// Set the simulation step-size, in seconds.
    void simStep( const realT &st /**< the new simulation step size*/ );

    /// Link this WFS to an AO system simulation
    template <typename AOSysT>
    void linkSystem( AOSysT &AOSys /**< The AO system simulation to link to*/ );

    /// Sense the wavefront aberrations
    /**  \returns true if a new wavefront measurement is ready.
     *  \returns false if still integrating.
     */
    bool senseWavefront( wavefrontT &pupilPlane /**< The input wavefront to be sensed*/ );

    /// Sense the wavefront aberrations in calibration mode
    /** Allows for faster calibrations.
     */
    bool senseWavefrontCal( wavefrontT &pupilPlane /**< The input wavefront to be sensed*/ );

  public:
    /// The WFS detector
    detectorT detector;

    /// The image on the detector, resized from m_wfsImage
    wfsImageT<realT> detectorImage;

    /// @}

    /** \name Pyramid Sensor Interface
     *
     * @{
     */
  protected:
    uint32_t m_nSides{ 4 }; ///< Number of sides in the pyramid

    /// The size of the pupil in wavefront pixels.
    /** This is the maximum diameter of the pupil in wavefront pixels.
     *
     */
    uint32_t m_pupilSz{ 0 };

    /// The separation of the pupil images in fraction of a pupil.  0 <= m_pupilSep, default 1.
    /** This sets the center-to-center separation of the pupils images in the focal plane wavefront.
     * Note that the separation in detector pixels then depends on the scaling between wavefront pixels
     * (m_wfSz) and detector pixels (m_detRows and m_detCols).
     *
     * This sets the size of the region in the pre-detection image that each pupil image
     * takes up, and therefore the size of the pre-detection image.
     * If the pupil (as defined in the input wavefront) is 60 pixels across
     * and m_pupilSep is set to 1.06667, then there will be a 2 pixel pad around each pupil image,
     * resulting in 4 pixels between each geometric pupil image.
     *
     * For a standard 4-sided pyramid, the pre-detection image will be
     * 2*m_pupilSep*m_pupilSz across.  For other n-sided pyramids, m_pupilSep still specifies the size of the pupil
     * image region, but the total image size will be a function of the resultant pupil positions.
     *
     * If m_pupilSep is less than 1, this will produce the "flattened pyramid", with overlap between
     * the pupil images.  In this case, image size will also be set by pupilSz to ensure that there are enough
     * pixels included to show all pupils.
     */
    realT m_pupilSep{ 1 };

    /// The angle by which to offset the pupils, in degrees. Default is 0.
    /** If this is 0, then a 4-sided pyramid makes a square as usual.  If this is set
     * to 45 degrees, then a 4-sided pyramid makes a diamond.
     *
     */
    realT m_angleOffset{ 0 };

    /// The size of the resulting PyWFS image in wavefront pixels.
    /** If \ref m_imageSzAuto is true, this is determined by number of sides (\ref m_nSides), the pupil size (\ref
     * m_pupilSz), and the pupil separation (\ref m_pupilSep).  For a 4 sided pyramid this will be the larger of
     * 2*m_pupilSep*m_pupilSz and 2*m_pupilSz.
     *
     * If , then this is used regardless of the optimum size.
     */
    uint32_t m_imageSz{ 0 };

    bool m_imageSzAuto{ true }; ///< Flag to track if \ref m_imageSz should be set to 0.

    realT m_wfPS{ 0 };          ///< Wavefront pixel scale, in meters/pixel

    realT m_D{ 0 };             ///< Telescope diameter, in meters

    uint32_t m_modSteps{ 20 };  ///< Number of modulation steps in one integration.  Can be set explicitly, but will be
                                ///< calculated if \ref m_perStep is set.

    realT m_perStep{ 1 };       /**< The minimum number of lamba/D per step in the modulation.
                                      Smaller will result in more steps.*/

    realT m_modRadius{ 3.0 };   ///< Radius of the modulation in pixels

  public:
    /// Get the number of pyramid sides
    /**
     * \returns the number of sides on the pyramid
     */
    int nSides();

    /// Set the number of sides on the pyramid
    /**
     */
    void nSides( const uint32_t &ns /**< The new number of sides on the pyramid*/ );

    /// Get the minimum number of modulation steps
    /**
     * \returns m_perStep;
     */
    realT perStep();

    /// Set the minimum number of modulation steps
    /**
     * \param mSt is the new number of modulation steps
     */
    void perStep( const realT &prStp /**< The minimum number of lamba/D per step to take*/ );

    /// Get the number of modulation steps
    /**
     * \returns m_modSteps, which is defined by perStep.
     */
    int modSteps();

    /// Get the radius of modulation
    /**
     * \returns m_modRadius;
     */
    realT modRadius();

    /// Set the modulation radius
    /**
     * \param mR is the new modulation radius in lambda/D
     */
    void modRadius( const realT &mR /**< [in] the new value of modulation radius */ );

    /// Get the wavefront pixel scale in meters per pixel
    /**
     * \returns the wavefront pixel scale in meters/pixel
     */
    realT wfPS();

    /// Get the telescope diameter
    /**
     * \returns the telescope diameter in meters
     */
    realT D();

    /// Set the telescope diameter
    /**
     * \param d is the new size in meters
     */
    void D( const realT &d /**< */ );

    /// Get the pupil size in pixels
    /** This is the pupil size in un-binned wavefront space
     *
     * \returns m_pupilSz
     * \returns m_pupilSz
     */
    uint32_t pupilSz();

    /// Set the pupil size in pixels.
    /** This is the size of the pupils in un-binned wavefront space.
     * See \ref m_pupilSz.
     *
     */
    void pupilSz( const uint32_t &sz /**< the new pupil size.*/ );

    /// Get the pupil separation as a fraction of pupil size
    /// Get the pupil separation as a fraction of pupil size
    /** This is the pupil separation in un-binned wavefront space
     *
     * \returns m_pupilSep
     */
    realT pupilSep();

    /// Set the pupil separation as a fraction of pupil size
    /// Set the pupil separation as a fraction of pupil size
    /** This is the separation of the pupils in un-binned wavefront space.
     * See \ref m_pupilSep.
     *
     */
    void pupilSep( const realT &sz /**< the new pupil separation.*/ );

    /// Get the angle offset
    /** See \ref m_angleOffset
     *
     * \returns m_angleOffset
     */
    realT angleOffset();

    /// Set the angle offset
    /** See \ref m_angleOffset.
     *
     */
    void angleOffset( const realT &ao /**< the new angle offset.*/ );

    /// Get the image size in wavefront pixels
    /** This is the size of the image in un-binned wavefront space
     *
     * \returns m_imageSz
     */
    uint32_t imageSz();

    /// Set the image size in wavefront pixels
    /** This is the size of the image in un-binned wavefront space
     * Setting a non-zero value also sets m_imageSizeAuto to false.
     * Setting 0 also sets m_imageSizeAuto to true.
     */
    void imageSz( const uint32_t &is );

    /// Get the value of the image size auto flag
    /** This controls whether image size is set automatically
     *
     * \returns m_imageSz
     */
    bool imageSzAuto();

    /// Set the value of the image size auto flag
    /** This controls whether image size is set automatically
     *
     */
    void imageSzAuto( const bool &ia );
    ///@}

    wfp::fraunhoferPropagator<complexFieldT, cudaGPU> m_frProp;

    bool m_opdMaskMade{ false };
    complexArrayT m_opdMask;

    bool m_tiltsMade{ false };
    std::vector<complexArrayT> m_tilts;

    bool m_preAllocated{ false };
    complexFieldT m_pupilPlaneCF;

    complexArrayT m_pupilPlaneCF_gpu;

    // Pre-allocated working memory:

    std::vector<complexArrayT> m_th_tiltedPlane; ///< Thread-local modulated wavefront

    std::vector<complexArrayT> m_th_focalPlane;  ///< Thread-local tip wavefront, used for FFT tilting

    std::vector<realArrayT> m_th_focalImage;     ///< Thread-local tip image

    std::vector<complexArrayT> m_th_sensorPlane; ///< Thread-local sensor-pupil-plane wavefront

    std::vector<realArrayT> m_th_sensorImage;    /**< Thread-local sensor-pupil-plane
                                                                             intensity image*/

    int m_iTime_counter{ 0 };

    int m_reading{ 0 };

    int m_roTime_counter{ 0 };

    std::vector<wavefrontT> _wavefronts;

    int m_lastWavefront{ 0 };

    realArrayT m_wfsTipImageAccum;
    realArrayT m_wfsImageAccum;

    /// The image formed by the WFS
    wfsImageT<realT> m_wfsImage;

  public:
    wfsImageT<realT> m_wfsTipImage;

  protected:
    void makeOpdMask();

    // These are called upload tilt but they are also used for opdMask. Maybe uploadField?
    template <int ccudaGPU = cudaGPU>
    error_t uploadTilt( complexArrayT &tilt, complexFieldT &ltilt, typename std::enable_if<ccudaGPU == 0>::type * = 0 );

    template <int ccudaGPU = cudaGPU>
    error_t uploadTilt( complexArrayT &tilt, complexFieldT &ltilt, typename std::enable_if<ccudaGPU == 1>::type * = 0 );

    void makeTilts();

    void allocThreadMem();

    void preAllocate();

    /// Convert a cpu complex field to a pointer to its CPU memory
    /** This just returns a pointer to the input array
     *
     */
    template <int ccudaGPU = cudaGPU>
    complexArrayT *uploadPupilPlaneCF( complexFieldT &cf /**< [in] the CPU complex field  */,
                                       typename std::enable_if<ccudaGPU == 0>::type * = 0 );

    /// Convert a cpu complex field to a pointer to its GPU memory
    /** This uploads the input array to the device and returns a cudaPtr.
     *
     */
    template <int ccudaGPU = cudaGPU>
    complexArrayT *uploadPupilPlaneCF( complexFieldT &cf /**< [in] the CPU complex field  */,
                                       typename std::enable_if<ccudaGPU == 1>::type * = 0 );

    /// Accumulate an image with a weight applied on the CPU
    /**
     * aim += im*w
     */
    template <int ccudaGPU = cudaGPU>
    void accumWeightedImage( realArrayT &aim, /**< [out] the image in which to accumulate the results*/
                             realArrayT &im,  /**< [in] the image to weight and then add to output*/
                             realT w,         /**< [in] the weight*/
                             typename std::enable_if<ccudaGPU == 0>::type * = 0 );

    /// Accumulate an image with a weight applied on the GPU
    /**
     * aim += im*w
     */
    template <int ccudaGPU = cudaGPU>
    void accumWeightedImage( realArrayT &aim, /**< [out] the image in which to accumulate the results*/
                             realArrayT &im,  /**< [in] the image to weight and then add to output*/
                             realT w,         /**< [in] the weight*/
                             typename std::enable_if<ccudaGPU == 1>::type * = 0 );

    /// Scale the accumulated image by number of accumulations and assign it to the output image. CPU version.
    /**
     * im = aim / naccums
     */
    template <int ccudaGPU = cudaGPU>
    void downloadAccumImage( realImageT &im,
                             realArrayT &aim,
                             uint32_t naccums,
                             typename std::enable_if<ccudaGPU == 0>::type * = 0 );

    /// Scale the accumulated image by number of accumulations and assign it to the output image. GPU version.
    /**
     * im = aim / naccums
     */
    template <int ccudaGPU = cudaGPU>
    void downloadAccumImage( realImageT &im,
                             realArrayT &aim,
                             uint32_t naccums,
                             typename std::enable_if<ccudaGPU == 1>::type * = 0 );

    void doSenseWavefront();
    void doSenseWavefront( wavefrontT & /**< */ );
    void doSenseWavefrontNoMod( wavefrontT & /**< */ );

    bool m_firstRun{ true };

  protected:
    // clang-format off
    #ifdef MXLIB_CUDA
    std::vector<mx::cuda::cublasHandle *> m_cublasHandle;
    #endif
    // clang-format on
};

template <typename realT, typename detectorT, int cudaGPU>
pyramidSensor<realT, detectorT, cudaGPU>::pyramidSensor()
{
    iTime( m_iTime );

    m_frProp.wholePixel( 0 );
}

template <typename realT, typename detectorT, int cudaGPU>
pyramidSensor<realT, detectorT, cudaGPU>::~pyramidSensor()
{
    // clang-format off
    #ifdef MXLIB_CUDA
    if(cudaGPU)
    {
        for(size_t nTh = 0; nTh < m_cublasHandle.size(); ++nTh)
        {
            delete m_cublasHandle[nTh];
        }
    }
    #endif // clang-format on
}

template <typename realT, typename detectorT, int cudaGPU>
int pyramidSensor<realT, detectorT, cudaGPU>::wfSz()
{
    return m_wfSz;
}

template <typename realT, typename detectorT, int cudaGPU>
void pyramidSensor<realT, detectorT, cudaGPU>::wfSz( const uint32_t &sz )
{
    if( m_wfSz == sz )
    {
        return;
    }

    m_wfSz = sz;

    m_tiltsMade = false;
    m_opdMaskMade = false;
    m_preAllocated = false;
}

template <typename realT, typename detectorT, int cudaGPU>
uint32_t pyramidSensor<realT, detectorT, cudaGPU>::detRows()
{
    return m_detRows;
}

template <typename realT, typename detectorT, int cudaGPU>
uint32_t pyramidSensor<realT, detectorT, cudaGPU>::detCols()
{
    return m_detCols;
}

template <typename realT, typename detectorT, int cudaGPU>
void pyramidSensor<realT, detectorT, cudaGPU>::detSize( const uint32_t &nrows, const uint32_t &ncols )
{
    if( m_detRows == nrows && m_detCols == ncols )
    {
        return;
    }

    m_detRows = nrows;
    m_detCols = ncols;

    detector.setSize( m_detRows, m_detCols );
    detectorImage.image.resize( m_detRows, m_detCols );

    m_opdMaskMade = false; // make sure size check is run on current settings
}

template <typename realT, typename detectorT, int cudaGPU>
realT pyramidSensor<realT, detectorT, cudaGPU>::lambda()
{
    return m_lambda;
}

template <typename realT, typename detectorT, int cudaGPU>
void pyramidSensor<realT, detectorT, cudaGPU>::lambda( const realT &l )
{
    m_lambda = l;

    //---------------------------------------
    //  Check if wavelength vector is filled
    //---------------------------------------
    if( m_wavelengths.size() == 0 )
    {
        m_wavelengths.resize( 1, m_lambda );
        _wavelengthWeights.resize( 1, 1.0 );
    }
}

template <typename realT, typename detectorT, int cudaGPU>
int pyramidSensor<realT, detectorT, cudaGPU>::iTime()
{
    return m_iTime;
}

template <typename realT, typename detectorT, int cudaGPU>
void pyramidSensor<realT, detectorT, cudaGPU>::iTime( const uint32_t &it )
{
    if( it < 1 )
    {
        throw mx::exception( error_t::invalidconfig, "iTime must be >= 1" );
    }

    m_iTime = it;

    _wavefronts.resize( m_iTime + 2 );
    m_lastWavefront = -1;

    detector.expTime( m_simStep * m_iTime );
}

template <typename realT, typename detectorT, int cudaGPU>
int pyramidSensor<realT, detectorT, cudaGPU>::roTime()
{
    return roTime;
}

template <typename realT, typename detectorT, int cudaGPU>
void pyramidSensor<realT, detectorT, cudaGPU>::roTime( const uint32_t &rt )
{
    if( rt < 1 )
    {
        throw mx::exception( error_t::invalidconfig, "roTime must be >= 1" );
    }

    m_roTime = rt;
}

template <typename realT, typename detectorT, int cudaGPU>
realT pyramidSensor<realT, detectorT, cudaGPU>::simStep()
{
    return m_simStep;
}

template <typename realT, typename detectorT, int cudaGPU>
void pyramidSensor<realT, detectorT, cudaGPU>::simStep( const realT &st )
{

    m_simStep = st;

    detector.expTime( m_simStep * m_iTime );
}

template <typename realT, typename detectorT, int cudaGPU>
template <typename AOSysT>
void pyramidSensor<realT, detectorT, cudaGPU>::linkSystem( AOSysT &AOSys )
{
    AOSys.wfs.wfPS( AOSys.m_wfPS );
    AOSys.wfs.D( AOSys.m_D );
}

template <typename realT, typename detectorT, int cudaGPU>
bool pyramidSensor<realT, detectorT, cudaGPU>::senseWavefront( wavefrontT &pupilPlane )
{

    ++m_lastWavefront;
    if( m_lastWavefront >= _wavefronts.size() )
        m_lastWavefront = 0;
    _wavefronts[m_lastWavefront].amplitude = pupilPlane.amplitude;
    _wavefronts[m_lastWavefront].phase = pupilPlane.phase;
    _wavefronts[m_lastWavefront].iterNo = pupilPlane.iterNo;

    // Always skip the first one for averaging to center of iTime.
    if( m_firstRun )
    {
        m_firstRun = false;
        return false;
    }

    ++m_iTime_counter;

    bool rv = false;

    if( m_reading )
    {
        ++m_roTime_counter;

        if( m_roTime_counter >= m_roTime )
        {
            detector.exposeImage( detectorImage.image, m_wfsImage.image );

            detectorImage.tipImage = m_wfsTipImage.image;
            detectorImage.iterNo = m_wfsImage.iterNo;

            m_roTime_counter = 0;
            m_reading = 0;
            rv = true;
        }
    }

    if( m_iTime_counter >= m_iTime )
    {
        doSenseWavefront();
        m_iTime_counter = 0;

        m_reading = 1;
        m_roTime_counter = 0;
    }

    return rv;
}

template <typename realT, typename detectorT, int cudaGPU>
bool pyramidSensor<realT, detectorT, cudaGPU>::senseWavefrontCal( wavefrontT &pupilPlane )
{

    BREAD_CRUMB;

    doSenseWavefront( pupilPlane );

    BREAD_CRUMB;

    detector.exposeImage( detectorImage.image, m_wfsImage.image );

    BREAD_CRUMB;

    detectorImage.tipImage = m_wfsTipImage.image;

    BREAD_CRUMB;

    return true;
}

/* Pyramid Specifics */

template <typename realT, typename detectorT, int cudaGPU>
int pyramidSensor<realT, detectorT, cudaGPU>::nSides()
{
    return m_nSides;
}

template <typename realT, typename detectorT, int cudaGPU>
void pyramidSensor<realT, detectorT, cudaGPU>::nSides( const uint32_t &ns )
{
    m_nSides = ns;
    m_opdMaskMade = false;
}

template <typename realT, typename detectorT, int cudaGPU>
realT pyramidSensor<realT, detectorT, cudaGPU>::wfPS()
{
    return m_wfPS;
}

template <typename realT, typename detectorT, int cudaGPU>
realT pyramidSensor<realT, detectorT, cudaGPU>::D()
{
    return m_D;
}

template <typename realT, typename detectorT, int cudaGPU>
void pyramidSensor<realT, detectorT, cudaGPU>::D( const realT &d )
{
    m_D = d;

    if( m_pupilSz >
        0 ) // Avoid making inf or nan so m_wfPS remains unset.  note that fast-math make detecting inf and nan hard.
    {
        m_wfPS = m_D / m_pupilSz;
    }
    else
    {
        m_wfPS = 0;
    }

    m_tiltsMade = false;
    m_opdMaskMade = false;
    m_preAllocated = false;
    m_opdMaskMade = false;
    m_preAllocated = false;
}

template <typename realT, typename detectorT, int cudaGPU>
realT pyramidSensor<realT, detectorT, cudaGPU>::perStep()
{
    return m_perStep;
}

template <typename realT, typename detectorT, int cudaGPU>
void pyramidSensor<realT, detectorT, cudaGPU>::perStep( const realT &prStp )
{
    m_perStep = prStp;

    if( m_modRadius <= 0 )
    {
        m_modSteps = 0;
        return;
    }

    realT radPerStep = m_perStep / m_modRadius;

    // Get the minimum number of steps to meet m_perStep while having symmetry for the quadrants
    m_modSteps = 1;
    while( math::half_pi<realT>() / m_modSteps > radPerStep )
    {
        ++m_modSteps;
    }

    m_modSteps *= 4;

    m_tiltsMade = false;
}

template <typename realT, typename detectorT, int cudaGPU>
int pyramidSensor<realT, detectorT, cudaGPU>::modSteps()
{
    return m_modSteps;
}

template <typename realT, typename detectorT, int cudaGPU>
realT pyramidSensor<realT, detectorT, cudaGPU>::modRadius()
{
    return m_modRadius;
}

template <typename realT, typename detectorT, int cudaGPU>
void pyramidSensor<realT, detectorT, cudaGPU>::modRadius( const realT &mR )
{
    m_modRadius = mR;
    perStep( m_perStep ); // to calculate m_modSteps;

    m_tiltsMade = false;
}

template <typename realT, typename detectorT, int cudaGPU>
uint32_t pyramidSensor<realT, detectorT, cudaGPU>::pupilSz()
{
    return m_pupilSz;
}

template <typename realT, typename detectorT, int cudaGPU>
void pyramidSensor<realT, detectorT, cudaGPU>::pupilSz( const uint32_t &sz )
{
    if( m_pupilSz == sz )
    {
        return;
    }

    m_pupilSz = sz;

    if( m_pupilSz >
        0 ) // Avoid making inf or nan so m_wfPS remains unset.  note that fast-math make detecting inf and nan hard.
    {
        m_wfPS = m_D / m_pupilSz;
    }
    else
    {
        m_wfPS = 0;
    }

    if( m_pupilSz >
        0 ) // Avoid making inf or nan so m_wfPS remains unset.  note that fast-math make detecting inf and nan hard.
    {
        m_wfPS = m_D / m_pupilSz;
    }
    else
    {
        m_wfPS = 0;
    }

    m_tiltsMade = false;
    m_opdMaskMade = false;
    m_preAllocated = false;
}

template <typename realT, typename detectorT, int cudaGPU>
realT pyramidSensor<realT, detectorT, cudaGPU>::pupilSep()
{
    return m_pupilSep;
}

template <typename realT, typename detectorT, int cudaGPU>
void pyramidSensor<realT, detectorT, cudaGPU>::pupilSep( const realT &sz )
{
    if( m_pupilSep == sz )
    {
        return;
    }

    m_pupilSep = sz;

    m_opdMaskMade = false;
    m_preAllocated = false;
}

template <typename realT, typename detectorT, int cudaGPU>
realT pyramidSensor<realT, detectorT, cudaGPU>::angleOffset()
{
    return m_angleOffset;
}

template <typename realT, typename detectorT, int cudaGPU>
void pyramidSensor<realT, detectorT, cudaGPU>::angleOffset( const realT &ao )
{
    if( m_angleOffset == ao )
    {
        return;
    }

    m_angleOffset = ao;

    m_opdMaskMade = false;
    m_preAllocated = false;
}

template <typename realT, typename detectorT, int cudaGPU>
uint32_t pyramidSensor<realT, detectorT, cudaGPU>::imageSz()
{
    if( !m_opdMaskMade )
    {
        makeOpdMask();
    }

    return m_imageSz;
}

template <typename realT, typename detectorT, int cudaGPU>
void pyramidSensor<realT, detectorT, cudaGPU>::imageSz( const uint32_t &sz )
{
    if( m_imageSz == sz )
    {
        return;
    }

    m_imageSz = sz;

    if( m_imageSz == 0 )
    {
        m_imageSzAuto = true;
    }
    else
    {
        m_imageSzAuto = false;
    }

    m_opdMaskMade = false;
    m_preAllocated = false;
}

template <typename realT, typename detectorT, int cudaGPU>
bool pyramidSensor<realT, detectorT, cudaGPU>::imageSzAuto()
{
    return m_imageSzAuto;
}

template <typename realT, typename detectorT, int cudaGPU>
void pyramidSensor<realT, detectorT, cudaGPU>::imageSzAuto( const bool &ia )
{
    m_imageSzAuto = ia;

    m_opdMaskMade = false;

    m_preAllocated = false;
}

template <typename realT, typename detectorT, int cudaGPU>
void pyramidSensor<realT, detectorT, cudaGPU>::makeOpdMask()
{

    if( m_wfPS == 0 )
    {
        throw mx::exception( error_t::invalidconfig,
                             "wavefront platescale (m_wfPS) is 0. Must set pupilSz and D first." );
    }

    if( !std::isfinite( m_wfPS ) || !std::isnormal( m_wfPS ) )
    {
        throw mx::exception( error_t::invalidconfig,
                             "wavefront platescale (m_wfPS) is infinite. Must set pupilSz and D first." );
    }

    std::cerr << m_wfPS << " " << m_D << "\n";

    if( m_D == 0 )
    {
        throw mx::exception( error_t::invalidconfig, "pupil diameter is 0. Must set D > 0 first." );
    }

    // Setup the Fraunhoffer Propagator
    m_frProp.setWavefrontSizePixels( m_wfSz );

    m_opdMask.resize( m_wfSz, m_wfSz );

    complexFieldT opdMask;
    opdMask.resize( m_wfSz, m_wfSz );

    complexFieldT opdMaskQ;
    opdMaskQ.resize( m_wfSz, m_wfSz );

    mx::improc::eigenImage<realT> mask;

    mask.resize( m_opdMask.rows(), m_opdMask.cols() );
    realT dang = mx::math::two_pi<realT>() / m_nSides;

    realT minx = 0;
    realT maxx = 0;
    realT miny = 0;
    realT maxy = 0;

    realT pupilRad = m_pupilSep * m_pupilSz / ( 2 * sin( dang / 2.0 ) );

    for( int n = 0; n < m_nSides; ++n )
    {
        realT ang = m_angleOffset * math::degreesT<realT>::radians + 0.5 * dang + n * dang;

        realT dx = pupilRad * cos( ang );

        if( dx < minx )
            minx = dx;
        if( dx > maxx )
            maxx = dx;

        realT dy = pupilRad * sin( ang );

        if( dy < miny )
            miny = dy;
        if( dy > maxy )
            maxy = dy;

        opdMaskQ.setConstant( std::complex<realT>( 0, 1 ) );
        wfp::tiltWavefront( opdMaskQ, dx, dy );
        mask.setZero();
        improc::maskWedge( mask,
                           0.5 * ( mask.rows() - 1 ),
                           0.5 * ( mask.cols() - 1 ),
                           math::rtod( ang ),
                           0.5 * math::rtod( dang ),
                           1 );
        wfp::extractMaskedPixels( opdMask, opdMaskQ, mask );
    }

    uploadTilt( m_opdMask, opdMask );

    int xsz =
        2 * std::max( { fabs( maxx ), fabs( minx ) } ) + 2 * std::max( { ( pupilRad / 2 ), ( (realT)m_pupilSz / 2 ) } );

    int ysz =
        2 * std::max( { fabs( maxy ), fabs( miny ) } ) + 2 * std::max( { ( pupilRad / 2 ), ( (realT)m_pupilSz / 2 ) } );

    if( m_imageSzAuto )
    {
        m_imageSz = std::max( xsz, ysz );

        if( m_pupilSep > 1 )
        {
            m_imageSz += ( m_pupilSep - 1.0 ) * m_pupilSz;
        }
    }

    if( m_imageSz > m_wfSz )
    {
        std::string msg = "image size (m_imageSz = " + std::to_string( m_imageSz ) + ") ";
        msg += "> wavefront size (m_wfSz = " + std::to_string( m_wfSz ) + "). ";
        msg += "Decrease number of sides (m_nSides = " + std::to_string( m_nSides ) + ") or increase wavefront size. ";
        throw mx::exception( error_t::invalidconfig, msg );
    }

    m_wfsImage.image.resize( m_imageSz, m_imageSz );
    m_wfsTipImage.image.resize( m_wfSz, m_wfSz );

    if( m_detRows == 0 || m_detCols == 0 )
    {
        detSize( m_imageSz, m_imageSz );
    }
    else if( m_detRows > m_imageSz || m_detCols > m_imageSz )
    {
        std::string msg = "image size (m_imageSz = " + std::to_string( m_imageSz ) + ") ";
        msg += "< detector size size (m_detRows = " + std::to_string( m_detRows );
        msg += " m_detCols = " + std::to_string( m_detCols ) + "). ";
        throw mx::exception( error_t::invalidconfig, msg );
    }

    m_opdMaskMade = true;
}

template <typename realT, typename detectorT, int cudaGPU>
template <int ccudaGPU>
error_t pyramidSensor<realT, detectorT, cudaGPU>::uploadTilt( complexArrayT &tilt,
                                                              complexFieldT &ltilt,
                                                              typename std::enable_if<ccudaGPU == 0>::type * )
{
    tilt = ltilt;
    return error_t::noerror;
}

template <typename realT, typename detectorT, int cudaGPU>
template <int ccudaGPU>
error_t pyramidSensor<realT, detectorT, cudaGPU>::uploadTilt( complexArrayT &tilt,
                                                              complexFieldT &ltilt,
                                                              typename std::enable_if<ccudaGPU == 1>::type * )
{
    return tilt.upload( ltilt.data(), ltilt.size() );
}

template <typename realT, typename detectorT, int cudaGPU>
void pyramidSensor<realT, detectorT, cudaGPU>::makeTilts()
{
    constexpr realT pi = math::pi<realT>();

    if( m_modSteps == 0 )
    {
        throw mx::exception( error_t::invalidconfig, "number of modulation steps (m_modSteps) has not been set." );
    }

    if( m_wfPS == 0 )
    {
        throw mx::exception( error_t::invalidconfig,
                             "wavefront platescale (m_wfPS) is 0. "
                             "Must set pupilSz and D first." );
    }

    if( !std::isfinite( m_wfPS ) )
    {
        throw mx::exception( error_t::invalidconfig,
                             "wavefront platescale (m_wfPS) is infinite. "
                             "Must set pupilSz and D first." );
    }

    if( m_D == 0 )
    {
        throw mx::exception( error_t::invalidconfig, "pupil diameter is 0. Must set D > 0 first." );
    }

    realT dang = 2 * pi / ( m_modSteps );
    realT dx, dy;

    m_tilts.resize( m_modSteps );

    std::cout << "WF Size: " << m_wfSz << "\n";
    std::cout << "WF PS:   " << m_wfPS << "\n";
    std::cout << "Lambda:  " << m_lambda << "\n";
    std::cout << "Pyr. PS: " << wfp::fftPlateScale<realT>( m_wfSz, m_wfPS, m_lambda ) * 206265. << " (mas/pix)\n";
    std::cout << "Mod. steps: " << m_modSteps << "\n";
    std::cout << "Mod rad: " << m_modRadius * ( m_lambda / m_D ) / wfp::fftPlateScale<realT>( m_wfSz, m_wfPS, m_lambda )
              << " (pixels)\n";

    for( int i = 0; i < m_modSteps; ++i )
    {
        dx = m_modRadius * ( m_lambda / m_D ) / wfp::fftPlateScale<realT>( m_wfSz, m_wfPS, m_lambda ) *
             cos( 0.0 * dang + dang * i );
        dy = m_modRadius * ( m_lambda / m_D ) / wfp::fftPlateScale<realT>( m_wfSz, m_wfPS, m_lambda ) *
             sin( 0.0 * dang + dang * i );

        complexFieldT tilt;
        tilt.resize( m_wfSz, m_wfSz );
        tilt.setConstant( std::complex<realT>( 0, 1 ) );

        wfp::tiltWavefront( tilt, dx, dy );

        uploadTilt( m_tilts[i], tilt );
    }

    m_tiltsMade = true;
}

template <typename realT, typename detectorT, int cudaGPU>
void pyramidSensor<realT, detectorT, cudaGPU>::allocThreadMem()
{
    if( !m_opdMaskMade )
    {
        makeOpdMask(); // Needed for m_imageSz
    }

    m_pupilPlaneCF.resize( m_wfSz, m_wfSz );

    if( cudaGPU )
    {
        m_pupilPlaneCF_gpu.resize( m_wfSz, m_wfSz );
    }

    int maxTh = omp_get_max_threads();
    m_th_tiltedPlane.resize( maxTh );

    m_th_focalPlane.resize( maxTh );

    m_th_focalImage.resize( maxTh );

    m_th_sensorPlane.resize( maxTh );

    m_th_sensorImage.resize( maxTh );

    // clang-format off
    #ifdef MXLIB_CUDA
    if(cudaGPU)
    {
        m_cublasHandle.resize(maxTh);
    }
    #endif // clang-format on

    for( int nTh = 0; nTh < maxTh; ++nTh )
    {
        m_th_tiltedPlane[nTh].resize( m_wfSz, m_wfSz );

        m_th_focalPlane[nTh].resize( m_wfSz, m_wfSz );

        m_th_focalImage[nTh].resize( m_wfSz, m_wfSz );

        m_th_sensorPlane[nTh].resize( m_wfSz, m_wfSz );

        m_th_sensorImage[nTh].resize( m_imageSz, m_imageSz );

        // clang-format off
        #ifdef MXLIB_CUDA
        if(cudaGPU)
        {
            m_cublasHandle[nTh] = new mx::cuda::cublasHandle(true);
        }
        #endif // clang-format on
    }

    m_wfsTipImageAccum.resize( m_wfSz, m_wfSz );
    m_wfsImageAccum.resize( m_imageSz, m_imageSz );

    m_preAllocated = true;
}

template <typename realT, typename detectorT, int cudaGPU>
void pyramidSensor<realT, detectorT, cudaGPU>::preAllocate()
{
    if( !m_tiltsMade )
    {
        makeTilts();
    }

    if( !m_opdMaskMade )
    {
        makeOpdMask();
    }

    if( !m_preAllocated )
    {
        allocThreadMem();
    }
}

template <typename realT, typename detectorT, int cudaGPU>
void pyramidSensor<realT, detectorT, cudaGPU>::doSenseWavefront()
{
    BREAD_CRUMB;

    wavefrontT pupilPlane;

    /* Here make average wavefront for now */
    int _firstWavefront = m_lastWavefront - m_iTime;
    if( _firstWavefront < 0 )
        _firstWavefront += _wavefronts.size();

    pupilPlane.amplitude = _wavefronts[_firstWavefront].amplitude;
    pupilPlane.phase = _wavefronts[_firstWavefront].phase;

    realT avgIt = _wavefronts[_firstWavefront].iterNo;

    BREAD_CRUMB;

    for( int i = 0; i < m_iTime; ++i )
    {
        ++_firstWavefront;
        if( (size_t)_firstWavefront >= _wavefronts.size() )
            _firstWavefront = 0;

        pupilPlane.amplitude += _wavefronts[_firstWavefront].amplitude;
        pupilPlane.phase += _wavefronts[_firstWavefront].phase;
        avgIt += _wavefronts[_firstWavefront].iterNo;
    }

    BREAD_CRUMB;

    pupilPlane.amplitude /= ( m_iTime + 1 );
    pupilPlane.phase /= ( m_iTime + 1 );

    avgIt /= ( m_iTime + 1.0 );

    /*=====================================*/
    doSenseWavefront( pupilPlane );

    m_wfsImage.iterNo = avgIt;
}

template <typename complexT>
void elWiseProduct_cpu( complexT *out, complexT *in1, complexT *in2, size_t nelem )
{
    for( int ii = 0; ii < nelem; ++ii )
    {
        out[ii] = in1[ii] * in2[ii];
    }
}

template <typename complexT>
void elWiseProduct_gpu( complexT *out, complexT *in1, complexT *in2, size_t nelem )
{
    BREAD_CRUMB;
    cuda::elementwiseXxY( out, in1, in2, nelem );
    BREAD_CRUMB;

    cudaError_t ce = cudaDeviceSynchronize();

    if( ce != cudaSuccess )
    {
        std::cerr << "cudaError " << cudaGetErrorName( ce ) << ": " << cudaGetErrorString( ce ) << '\n';
        return;
    }
    BREAD_CRUMB;
}

template <typename complexT, int cudaGPU>
void elWiseProduct( complexT *out, complexT *in1, complexT *in2, size_t nelem );

template <>
void elWiseProduct<std::complex<float>, 0>( std::complex<float> *out,
                                            std::complex<float> *in1,
                                            std::complex<float> *in2,
                                            size_t nelem )
{
    elWiseProduct_cpu( out, in1, in2, nelem );
}

template <>
void elWiseProduct<std::complex<double>, 0>( std::complex<double> *out,
                                             std::complex<double> *in1,
                                             std::complex<double> *in2,
                                             size_t nelem )
{
    elWiseProduct_cpu( out, in1, in2, nelem );
}

template <>
void elWiseProduct<std::complex<float>, 1>( std::complex<float> *out,
                                            std::complex<float> *in1,
                                            std::complex<float> *in2,
                                            size_t nelem )
{
    BREAD_CRUMB;
    elWiseProduct_gpu( reinterpret_cast<cuComplex *>( out ),
                       reinterpret_cast<cuComplex *>( in1 ),
                       reinterpret_cast<cuComplex *>( in2 ),
                       nelem );

    cudaError_t ce = cudaDeviceSynchronize();

    if( ce != cudaSuccess )
    {
        std::cerr << "cudaError " << cudaGetErrorName( ce ) << ": " << cudaGetErrorString( ce ) << '\n';
        return;
    }
    BREAD_CRUMB;
}

template <>
void elWiseProduct<std::complex<double>, 1>( std::complex<double> *out,
                                             std::complex<double> *in1,
                                             std::complex<double> *in2,
                                             size_t nelem )
{
    elWiseProduct_gpu( reinterpret_cast<cuDoubleComplex *>( out ),
                       reinterpret_cast<cuDoubleComplex *>( in1 ),
                       reinterpret_cast<cuDoubleComplex *>( in2 ),
                       nelem );
}

template <typename realT, typename detectorT, int cudaGPU>
template <int ccudaGPU>
pyramidSensor<realT, detectorT, cudaGPU>::complexArrayT *
pyramidSensor<realT, detectorT, cudaGPU>::uploadPupilPlaneCF( complexFieldT &cf,
                                                              typename std::enable_if<ccudaGPU == 0>::type * )
{
    return &cf;
}

template <typename realT, typename detectorT, int cudaGPU>
template <int ccudaGPU>
pyramidSensor<realT, detectorT, cudaGPU>::complexArrayT *
pyramidSensor<realT, detectorT, cudaGPU>::uploadPupilPlaneCF( complexFieldT &cf,
                                                              typename std::enable_if<ccudaGPU == 1>::type * )
{
    m_pupilPlaneCF_gpu.upload( cf.data() );
    return &m_pupilPlaneCF_gpu;
}

template <typename realT, typename detectorT, int cudaGPU>
template <int ccudaGPU>
void pyramidSensor<realT, detectorT, cudaGPU>::accumWeightedImage( realArrayT &aim,
                                                                   realArrayT &im,
                                                                   realT w,
                                                                   typename std::enable_if<ccudaGPU == 0>::type * )
{
    aim += im * w;
}

template <typename realT, typename detectorT, int cudaGPU>
template <int ccudaGPU>
void pyramidSensor<realT, detectorT, cudaGPU>::accumWeightedImage( realArrayT &aim,
                                                                   realArrayT &im,
                                                                   realT w,
                                                                   typename std::enable_if<ccudaGPU == 1>::type * )
{

    // clang-format off
    #ifdef MXLIB_CUDA

    mx::cuda::cublasTaxpy<realT>( *m_cublasHandle[omp_get_thread_num()], aim.size(), &w, im.data(), 1, aim.data(), 1 );

    #endif
    // clang-format on
}

template <typename realT, typename detectorT, int cudaGPU>
template <int ccudaGPU>
void pyramidSensor<realT, detectorT, cudaGPU>::downloadAccumImage( realImageT &im,
                                                                   realArrayT &aim,
                                                                   uint32_t naccums,
                                                                   typename std::enable_if<ccudaGPU == 0>::type * )
{
    im = aim / naccums;
}

template <typename realT, typename detectorT, int cudaGPU>
template <int ccudaGPU>
void pyramidSensor<realT, detectorT, cudaGPU>::downloadAccumImage( realImageT &im,
                                                                   realArrayT &aim,
                                                                   uint32_t naccums,
                                                                   typename std::enable_if<ccudaGPU == 1>::type * )
{
    aim.download( im.data() );
    im /= naccums;
}

template <typename realT, typename detectorT, int cudaGPU>
void pyramidSensor<realT, detectorT, cudaGPU>::doSenseWavefront( wavefrontT &pupilPlane )
{
    if( m_modRadius == 0 )
    {
        return doSenseWavefrontNoMod( pupilPlane );
    }

    if( !m_preAllocated )
    {
        preAllocate();
    }

    BREAD_CRUMB;

    BREAD_CRUMB;

    m_wfsTipImageAccum.setZero();
    m_wfsImageAccum.setZero();

    BREAD_CRUMB;

    for( size_t l = 0; l < m_wavelengths.size(); ++l )
    {
        pupilPlane.lambda = m_lambda;
        pupilPlane.getWavefront( m_pupilPlaneCF, m_wavelengths[l], m_wfSz );

        complexArrayT *pupilPlaneCF = uploadPupilPlaneCF( m_pupilPlaneCF );

        // clang-format off
        #pragma omp parallel// clang-format on
        {
            int nTh = omp_get_thread_num();

            // CUDA: these are cudaPtrs
            m_th_sensorImage[nTh].setZero();
            m_th_focalImage[nTh].setZero();

            BREAD_CRUMB;

            complexT *ppm_Data;
            complexT *tpm_Data;
            complexT *fpm_Data;
            complexT *opdm_Data;
            complexT *tim_Data;

            ppm_Data = reinterpret_cast<complexT *>( pupilPlaneCF->data() );
            tpm_Data = reinterpret_cast<complexT *>( m_th_tiltedPlane[nTh].data() );
            fpm_Data = reinterpret_cast<complexT *>( m_th_focalPlane[nTh].data() );
            opdm_Data = reinterpret_cast<complexT *>( m_opdMask.data() );

            int nelem = m_wfSz * m_wfSz;

            // clang-format off
            #pragma omp for // clang-format on
            for( int i = 0; i < m_modSteps; ++i )
            {

                tim_Data = reinterpret_cast<complexT*>(m_tilts[i].data());

                BREAD_CRUMB;

                //---------------------------------------------
                // Apply the modulating tip
                //---------------------------------------------
                elWiseProduct<complexT, cudaGPU>( tpm_Data, ppm_Data, tim_Data, nelem );

                BREAD_CRUMB;

                //---------------------------------------------
                // Propagate to Pyramid tip
                //---------------------------------------------
                m_frProp.propagatePupilToFocal( m_th_focalPlane[nTh], m_th_tiltedPlane[nTh], true );

                BREAD_CRUMB;

                //---------------------------------------------
                // Extract the tip image.
                //---------------------------------------------
                wfp::extractIntensityImageAccum<realArrayT, complexArrayT, cudaGPU>( m_th_focalImage[nTh],
                                                                                     0,
                                                                                     m_wfSz,
                                                                                     0,
                                                                                     m_wfSz,
                                                                                     m_th_focalPlane[nTh],
                                                                                     0,
                                                                                     0 );

                BREAD_CRUMB;

                //---------------------------------------------
                // Now apply the pyramid OPD
                //---------------------------------------------
                elWiseProduct<complexT, cudaGPU>( fpm_Data, fpm_Data, opdm_Data, nelem );

                BREAD_CRUMB;

                //---------------------------------------------
                // Propagate to sensor plane
                //---------------------------------------------
                m_frProp.propagateFocalToPupil( m_th_sensorPlane[nTh], m_th_focalPlane[nTh], true );

                BREAD_CRUMB;

                //---------------------------------------------
                // Extract the image.
                //---------------------------------------------

                wfp::extractIntensityImageAccum<realArrayT, complexArrayT, cudaGPU>( m_th_sensorImage[nTh],
                                                                                     0,
                                                                                     m_imageSz,
                                                                                     0,
                                                                                     m_imageSz,
                                                                                     m_th_sensorPlane[nTh],
                                                                                     0.5 * m_wfSz - m_imageSz / 2,
                                                                                     0.5 * m_wfSz - m_imageSz / 2 );



            } // for

            BREAD_CRUMB;

            // clang-format off
            #pragma omp critical // clang-format on
            {
                accumWeightedImage( m_wfsTipImageAccum, m_th_focalImage[nTh], _wavelengthWeights[l] );

                accumWeightedImage( m_wfsImageAccum, m_th_sensorImage[nTh], _wavelengthWeights[l] );
            }

        } // #pragma omp parallel

    } // l for wavelength

    BREAD_CRUMB;

    downloadAccumImage( m_wfsTipImage.image, m_wfsTipImageAccum, m_modSteps );
    downloadAccumImage( m_wfsImage.image, m_wfsImageAccum, m_modSteps );

    BREAD_CRUMB;

}

template <typename T>
void memCopy_cpu( T *out, T *in, size_t nelem )
{
    memcpy( out, in, nelem * sizeof( T ) );
}

template <typename T>
void memCopy_gpu( T *out, T *in, size_t nelem )
{
    cudaMemcpy( out, in, nelem * sizeof( T ), cudaMemcpyDeviceToDevice );
}

template <typename T, int cudaGPU>
void memCopy( T *out, T *in, size_t nelem );

template <>
void memCopy<std::complex<float>, 0>( std::complex<float> *out, std::complex<float> *in, size_t nelem )
{
    memCopy_cpu( out, in, nelem );
}

template <>
void memCopy<std::complex<double>, 0>( std::complex<double> *out, std::complex<double> *in, size_t nelem )
{
    memCopy_cpu( out, in, nelem );
}

template <>
void memCopy<std::complex<float>, 1>( std::complex<float> *out, std::complex<float> *in, size_t nelem )
{
    memCopy_gpu( out, in, nelem );
}

template <>
void memCopy<std::complex<double>, 1>( std::complex<double> *out, std::complex<double> *in, size_t nelem )
{
    memCopy_gpu( out, in, nelem );
}

template <typename realT, typename detectorT, int cudaGPU>
void pyramidSensor<realT, detectorT, cudaGPU>::doSenseWavefrontNoMod( wavefrontT &pupilPlane )
{
#ifndef MXLIB_CUDA
    BREAD_CRUMB;

    if( !m_opdMaskMade )
    {
        makeOpdMask();
    }

    m_wfsImage.image.resize( m_imageSz, m_imageSz );
    m_wfsImage.image.setZero();

    m_wfsTipImage.image.resize( m_wfSz, m_wfSz );
    m_wfsTipImage.image.setZero();

    pupilPlane.getWavefront( m_pupilPlaneCF, m_wfSz );

    complexArrayT *pupilPlaneCF = uploadPupilPlaneCF( m_pupilPlaneCF );

    complexFieldT tiltedPlane;
    complexFieldT focalPlane;
    complexFieldT sensorPlane;

    tiltedPlane.resize( m_wfSz, m_wfSz );
    focalPlane.resize( m_wfSz, m_wfSz );
    sensorPlane.resize( m_wfSz, m_wfSz );

    int nelem = m_wfSz * m_wfSz;

    complexT *tpm_Data = reinterpret_cast<complexT *>( tiltedPlane.data() );
    complexT *ppm_Data = reinterpret_cast<complexT *>( pupilPlaneCF->data() );
    complexT *opdm_Data = reinterpret_cast<complexT *>( m_opdMask.data() );
    complexT *fpm_Data = reinterpret_cast<complexT *>( focalPlane.data() );

    BREAD_CRUMB;

    //---------------------------------------------
    // Apply NO modulator tilt
    //---------------------------------------------
    memCopy<complexT, cudaGPU>( tpm_Data, ppm_Data, nelem );
    /*for( int ii = 0; ii < nelem; ++ii )
    {
        tpm_Data[ii] = ppm_Data[ii];
    }*/

    BREAD_CRUMB;

    //---------------------------------------------
    // Propagate to Pyramid tip
    //---------------------------------------------
    m_frProp.propagatePupilToFocal( focalPlane, tiltedPlane, true );

    BREAD_CRUMB;

    //---------------------------------------------
    // Extract the tip image.
    //---------------------------------------------
    wfp::extractIntensityImageAccum<realArrayT, complexArrayT, cudaGPU>( m_wfsTipImage.image,
                                                                         0,
                                                                         m_wfSz,
                                                                         0,
                                                                         m_wfSz,
                                                                         focalPlane,
                                                                         0,
                                                                         0 );

    BREAD_CRUMB;

    //---------------------------------------------
    // Now apply the pyramid OPD
    //---------------------------------------------
    elWiseProduct<complexT, cudaGPU>( fpm_Data, fpm_Data, opdm_Data, nelem );
    /*for( int ii = 0; ii < nelem; ++ii )
    {
        fpm_Data[ii] = fpm_Data[ii] * opdm_Data[ii];
    }*/

    BREAD_CRUMB;

    //---------------------------------------------
    // Propagate to sensor plane
    //---------------------------------------------
    m_frProp.propagateFocalToPupil( sensorPlane, focalPlane, true );

    BREAD_CRUMB;

    //---------------------------------------------
    // Extract the image.
    //---------------------------------------------
    wfp::extractIntensityImageAccum<realArrayT, complexArrayT, cudaGPU>( m_wfsImage.image,
                                                                         0,
                                                                         m_imageSz,
                                                                         0,
                                                                         m_imageSz,
                                                                         sensorPlane,
                                                                         0.5 * m_wfSz - m_imageSz / 2,
                                                                         0.5 * m_wfSz - m_imageSz / 2 );

    BREAD_CRUMB;

#endif
}

} // namespace sim
} // namespace AO
} // namespace mx

#endif // mx_AO_sim_pyramidSensor_hpp
