/** \file sigproc.hpp
 * \author Jared R. Males
 * \brief Library include for the sigproc module
 * \ingroup signal_processing_files
 *
 */

//***********************************************************************//
// Copyright 2023 Jared R. Males (jaredmales@gmail.com)
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

#ifndef sigproc_hpp
#define sigproc_hpp

#include "autocorrelation.hpp"
#include "averagePeriodogram.hpp"
#include "basisUtils2D.hpp"
#include "circularBuffer.hpp"
#include "fourierModes.hpp"
#include "gramSchmidt.hpp"
#include "levinsonRecursion.hpp"
#include "linearPredictor.hpp"
#include "psdFilter.hpp"
#include "psdUtils.hpp"
#include "psdVarMean.hpp"
#include "signalWindows.hpp"
#include "zernike.hpp"

#endif // sigproc_hpp
