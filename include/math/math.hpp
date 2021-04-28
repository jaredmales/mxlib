/** \file math.hpp
  * \author Jared R. Males
  * \brief Library include for the math module.
  * \ingroup gen_math_files
  *
  */

//***********************************************************************//
// Copyright 2020 Jared R. Males (jaredmales@gmail.com)
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

#ifndef math_hpp
#define math_hpp

#include "cuda/cudaPtr.hpp"
#include "cuda/templateCublas.hpp"
#include "cuda/templateCuda.hpp"
#include "cuda/templateCufft.hpp"
#include "cuda/templateCurand.hpp"
#include "fft/fft.hpp"
#include "fft/fftwEnvironment.hpp"
#include "fft/fftwTemplates.hpp"
#include "fit/array2FitGaussian2D.hpp"
#include "fit/fitAiry.hpp"
#include "fit/fitGaussian.hpp"
#include "fit/fitMoffat.hpp"
#include "fit/levmarInterface.hpp"
#include "fit/templateLevmar.hpp"
#include "func/airyPattern.hpp"
#include "func/bessel.hpp"
#include "func/factorial.hpp"
#include "func/gamma.hpp"
#include "func/gaussian.hpp"
#include "func/jinc.hpp"
#include "func/logistic.hpp"
#include "func/moffat.hpp"
#include "func/precision.hpp"
#include "func/sign.hpp"
#include "func/weibull.hpp"
#include "plot/gnuPlot.hpp"
#include "constants.hpp"
#include "eigenLapack.hpp"
#include "geo.hpp"
#include "gslInterpolation.hpp"
#include "histogramUniform.hpp"
#include "randomSeed.hpp"
#include "randomT.hpp"
#include "roots.hpp"
#include "templateBLAS.hpp"
#include "templateLapack.hpp"
#include "vectorUtils.hpp"

#endif //math_hpp
