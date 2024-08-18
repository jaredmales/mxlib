
#include <gsl/gsl_integration.h>
#include <gsl/gsl_errno.h>

#include "constants.hpp"
#include "logInterpolator.hpp"

namespace mx
{
namespace math
{

template <typename interpT>
typename interpT::realT logRadProfIntegrationF( typename interpT::realT x, void *params )
{
    mx::math::logInterpolator<interpT> *lI = (mx::math::logInterpolator<interpT> *)params;

    return ( *lI )(x)*x;
}

/// Integrate a numerical radial profile using logarithmic interpolation
/** Useful for steep power-law like functions like power-spectra
 *
 * \returns the value of the integral over the limits
 */
template <typename interpT>
typename interpT::realT
logRadProfIntegrator( const std::vector<typename interpT::realT>
                          &x, ///< [in] the x values of the function.  Must be positive definite (can not contain 0).
                      const std::vector<typename interpT::realT>
                          &y, ///< [in] the y values of the function.  Must be positive definite (can not contain 0).
                      typename interpT::realT x0, ///< [in] [optional] the lower limit of integration
                      typename interpT::realT xf  ///< [in] [optional] the uper limit of integration
)
{
    typedef typename interpT::realT realT;

    mx::math::logInterpolator<interpT> logInterp;

    logInterp.setup( x, y );
    gsl_function func;
    func.function = &logRadProfIntegrationF<interpT>;
    func.params = &logInterp;

    realT result;
    realT abserr;
    size_t neval;

    gsl_set_error_handler_off();
    gsl_integration_workspace *w = gsl_integration_workspace_alloc( 1e6 );

    int rv = gsl_integration_qag( &func, x0, xf, 1e-3, 1e-3, 1e6, 6, w, &result, &abserr );

    result *= two_pi<realT>();

    gsl_integration_workspace_free( w );

    return result;
}

/// Integrate a numerical radial profile using logarithmic interpolation
/** Useful for steep power-law like functions like power-spectra
 *
 * \returns the value of the integral over the entire domain given by \p x
 */
template <typename interpT>
typename interpT::realT
logRadProfIntegrator( const std::vector<typename interpT::realT>
                          &x, ///< [in] the x values of the function.  Must be positive definite (can not contain 0).
                      const std::vector<typename interpT::realT>
                          &y ///< [in] the y values of the function.  Must be positive definite (can not contain 0).
)
{
    return logRadProfIntegrator<interpT>( x, y, x[0], x.back() );
}

} // namespace math
} // namespace mx
