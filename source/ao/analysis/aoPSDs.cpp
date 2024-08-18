/** \file aoPSDs.cpp
 * \author Jared R. Males (jaredmales@gmail.com)
 * \brief Implementation of spatial power spectra used in adaptive optics.
 * \ingroup mxAO_files
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

#include "ao/analysis/aoPSDs.hpp"

namespace mx
{
namespace AO
{
namespace analysis
{

namespace PSDComponent
{

std::string compName( int cc )
{
    if( cc == phase )
        return "phase";
    if( cc == amplitude )
        return "amplitude";
    if( cc == dispPhase )
        return "dispPhase";
    if( cc == dispAmplitude )
        return "dispAmplitude";

    return "unknown";
}

int compNum( const std::string &name )
{
    if( name == "phase" )
        return phase;
    else if( name == "amplitude" )
        return amplitude;
    else if( name == "dispPhase" )
        return dispPhase;
    else if( name == "dispAmplitude" )
        return dispAmplitude;

    return -1;
}

} // namespace PSDComponent

template struct vonKarmanSpectrum<float>;

template struct vonKarmanSpectrum<double>;

template struct vonKarmanSpectrum<long double>;

#ifdef HASQUAD
template struct vonKarmanSpectrum<__float128>;
#endif

} // namespace analysis
} // namespace AO
} // namespace mx
