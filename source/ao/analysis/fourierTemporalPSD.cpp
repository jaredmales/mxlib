/** \file fourierTemporalPSD.cpp
  * \author Jared R. Males (jaredmales@gmail.com)
  * \brief Implementation of calculation of the temporal PSD of Fourier modes.
  * \ingroup mxAOm_files
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

#include "ao/analysis/fourierTemporalPSD.hpp"
#include "ao/analysis/aoSystem.hpp"

namespace mx
{
namespace AO
{
namespace analysis
{


/*template
struct fourierTemporalPSD<float, aoSystem<float, vonKarmanSpectrum<float>, std::ostream>>;*/

template
struct fourierTemporalPSD<double, aoSystem<double, vonKarmanSpectrum<double>, std::ostream>>;

/*
template
struct fourierTemporalPSD<long double, aoSystem<long double, vonKarmanSpectrum<long double>, std::ostream>>;

#ifdef HASQUAD
template
struct fourierTemporalPSD<__float128, aoSystem<__float128, vonKarmanSpectrum<__float128>, std::ostream>>;
#endif
*/

} //namespace analysis
} //namespace AO
} //namespace mx

