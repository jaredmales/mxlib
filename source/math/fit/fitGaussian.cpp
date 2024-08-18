/** \file fitGaussian.cpp
 * \author Jared R. Males
 * \brief Tools for fitting Gaussians to data.
 * \ingroup fitting_files
 *
 */

//***********************************************************************//
// Copyright 2021 Jared R. Males (jaredmales@gmail.com)
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

#include "math/fit/fitGaussian.hpp"

namespace mx
{
namespace math
{
namespace fit
{

template class fitGaussian1D<float>;
template class fitGaussian1D<double>;

template struct gaussian1D_fitter<float>;
template struct gaussian1D_fitter<double>;

template class fitGaussian2D<mx::math::fit::gaussian2D_sym_fitter<float>>;
template class fitGaussian2D<mx::math::fit::gaussian2D_sym_fitter<double>>;
template class fitGaussian2D<mx::math::fit::gaussian2D_gen_fitter<float>>;
template class fitGaussian2D<mx::math::fit::gaussian2D_gen_fitter<double>>;

template int guessGauss2D_ang<float>( float &Ag,
                                      float &xg,
                                      float &yg,
                                      float &xFWHM,
                                      float &yFWHM,
                                      float &angG,
                                      mx::improc::eigenImage<float> &im,
                                      float maxWidth,
                                      float widthWidth,
                                      float nAngs,
                                      float xg0,
                                      float yg0 );

template int guessGauss2D_ang<double>( double &Ag,
                                       double &xg,
                                       double &yg,
                                       double &xFWHM,
                                       double &yFWHM,
                                       double &angG,
                                       mx::improc::eigenImage<double> &im,
                                       double maxWidth,
                                       double widthWidth,
                                       double nAngs,
                                       double xg0,
                                       double yg0 );

} // namespace fit
} // namespace math

} // namespace mx
