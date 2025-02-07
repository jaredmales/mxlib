/** \file point2D.hpp
 * \author Jared R. Males
 * \brief A 2D point in space
 * \ingroup gen_math_files
 *
 */

//***********************************************************************//
// Copyright 2025 Jared R. Males (jaredmales@gmail.com)
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

#ifndef math_point2D_hpp
#define math_point2D_hpp

namespace mx
{
namespace math
{

/// A point in 2-dimensional space
/**
 *
 * \tparam _angleT specifies the angle units, either radiansT<realT> or degreesT<realT>.
 */
template <typename _angleT>
class point2D
{

  public:
    typedef _angleT angleT;
    typedef typename angleT::realT realT;

    realT x; ///< The x-coordinate
    realT y; ///< The y-coordinate

    /// Default c'tor
    point2D();

    /// Construct with values for x and y
    point2D( realT nx, /**< [in] the x-coordinate*/
             realT ny /**< [in] the y-coordinate*/ );

    /// Construct with values for x and y and rotate them by an angle
    point2D( realT nx, /**< [in] the x-coordinate*/
             realT ny, /**< [in] the y-coordinate*/
             realT ang /**< [in] the angle by which to rotate.  Degrees or radians specified by angleT*/ );

    /// Construct with values for x and y and rotate them by an angle specified by cosine and sine
    point2D( realT nx, /**< [in] the x-coordinate*/
             realT ny, /**< [in] the y-coordinate*/
             realT c,  /**< [in] cosine of the angle by which to rotate.*/
             realT s /**< [in] sine of the angle by which to rotate.*/ );

    /// Rotate the point about the origin by a given angle
    /**
     * \returns a point2D at the rotated coordinates
     */
    point2D rotate( realT ang /**< [in] the angle by which to rotate.  Degrees or radians specified by angleT*/ );

    /// Rotate the point about the origin by a given angle
    /** In-place: operates on the current point, changing its coordinates
     *
    */
    void rotateInPlace( realT ang /**< [in] the angle by which to rotate.  Degrees or radians specified by angleT*/ );

    /// Rotate the point about the origin by an angle specified by its cosine and sine
    /**
     * \returns a point2D at the rotated coordinates
     */
    point2D rotate( realT c, /**< [in] cosine of the angle by which to rotate.*/
                    realT s /**< [in] sine of the angle by which to rotate.*/ );

    /// Rotate the point about the origin by an angle specified by its cosine and sine
    /** In-place: operates on the current point, changing its coordinates
     *
    */
    void rotateInPlace( realT c, /**< [in] cosine of the angle by which to rotate.*/
                        realT s /**< [in] sine of the angle by which to rotate.*/ );
};

template <typename angleT>
point2D<angleT>::point2D() : x( 0 ), y( 0 )
{
}

template <typename angleT>
point2D<angleT>::point2D( realT nx, realT ny ) : x( nx ), y( ny )
{
}

template <typename angleT>
point2D<angleT>::point2D( realT nx, realT ny, realT ang ) : x( nx ), y( ny )
{
    rotateInPlace( ang );
}

template <typename angleT>
point2D<angleT>::point2D( realT nx, realT ny, realT c, realT s ) : x( nx ), y( ny )
{
    x = nx * c - ny * s;
    y = nx * s + ny * c;

}

template <typename angleT>
point2D<angleT> point2D<angleT>::rotate( realT ang )
{
    realT c = cos( ang * angleT::radians );
    realT s = sin( ang * angleT::radians );

    return rotate( c, s );
}

template <typename angleT>
void point2D<angleT>::rotateInPlace( realT ang )
{
    realT c = cos( ang * angleT::radians );
    realT s = sin( ang * angleT::radians );

    return rotateInPlace( c, s );
}

template <typename angleT>
point2D<angleT> point2D<angleT>::rotate( realT c, realT s )
{
    point2D<angleT> p;

    p.x = x * c - y * s;
    p.y = x * s + y * c;

    return p;
}

template <typename angleT>
void point2D<angleT>::rotateInPlace( realT c, realT s )
{
    realT tx = x * c - y * s;
    realT ty = x * s + y * c;

    x = tx;
    y = ty;
}

} // namespace math
} // namespace mx

#endif // math_point2D_hpp
