//
//
// MIT License
//
// Copyright (c) 2020 Stellacore Corporation.
//
// Permission is hereby granted, free of charge, to any person obtaining
// a copy of this software and associated documentation files (the
// "Software"), to deal in the Software without restriction, including
// without limitation the rights to use, copy, modify, merge, publish,
// distribute, sublicense, and/or sell copies of the Software, and to
// permit persons to whom the Software is furnished to do so, subject
// to the following conditions:
//
// The above copyright notice and this permission notice shall be
// included in all copies or substantial portions of the Software.
//
// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY
// KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE
// WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
// NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS
// BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN
// AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR
// IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
// THE SOFTWARE.
//
//


#ifndef peridetic_INCL_
#define peridetic_INCL_


#include <array>
#include <limits>


/*! \brief Enclosing namespace for Peridetic functions (and implementation).
 *
 * This peridetic.h header file provides transformation of Geographic
 * from/into Cartesian coordinates with algorithms suitable for use
 * in applications involving:
 * \arg Terrestrial - Anywhere on land
 * \arg Bathymetric - Throughout all ocean depths
 * \arg Atmospheric - To the edge of space (also beyond, but less optimally)
 *
 * The algorithms are designed and tested specifically for use at
 * altitudes within +/-100[km] of an Earth ellipsoid surface. Within this
 * domain the transformations:
 * \arg Provide sub[mm] accuracy
 *		- Self consistent (forward/backward transformations) within 7.5[nm]
 *		- Agree with NOAA CORS station values within the 1[mm] published
 *		report precision (for 10 select stations from around globe)
 * \arg Are easy use (e.g. copy two header files, compile, and go)
 * \arg Are quite computationally fast and efficient (ref note below)
 *
 * The principal transformations are:
 * \arg peri::lpaForXyz() - Geographic Lon/Par(lat)/Alt from Cartesian X/Y/Z
 * \arg peri::xyzForLpa() - Cartesian X/Y/Z from Geographic Lon/Par/Alt
 *
 * Principal data structures are standard C++ aggregates:
 * \arg peri::XYZ - std::array<double, 3u> interpreted
 *    as "xMeters", "yMeters", "zMeters"
 * \arg peri::LPA - std::array<double, 3u> interpreted
 *    as "lonRadians", "par(lat)Radians", "altMeters"
 *
 * To use the transformations in consuming code:
 * \arg Copy header files ("peridetic.h" and also "periDetail.h")
 *   to desired location
 * \arg Specify that location in compile path
 * \arg Include "peridetic.h" in your source code as with any other header
 * \arg Call the transformation functions as desired.
 *
 * \note This file (perdetic.h) directly includes the implementation detail
 * header file, periDetail.h, which also must be in the compiler's include
 * search path.
 *
 * A complete and fully functioning header-only installation comprises
 * the following:
 * \arg peridetic.h - Public interface specification
 *    (which includes periDetail.h)
 * \arg periDetail.h - Underlying implementation
 *    (which includes only std:: namespace header files)
 *
 * \note For \b performance compile *WITH OPTIMIZATION* enabled. The
 * implementation detail involves many, many inline functions and methods.
 * If these are not optimized away, running times can be (many) order(s) of
 * magnitude slower!!
 *
 */
namespace peri
{
} // [peri]


namespace peri
{
	/*! \brief Type alias for Cartesian coordinate triple (aka ECEF)
	 *
	 * Coordinate are relative to origin at center of reference ellipsoid.
	 * \arg Units for expressing all three coordinate values are \b meters.
	 * \arg I.e. Create along the lines of:
	 * \code XYZ const xyzLoc{ xMeters, yMeters, zMeters };
	 * \endcode
	 *
	 * For example, given equatorial radius, aRad, and polar radius bRad:
	 * \arg North pole at XYZ{ 0., 0.,  bRad }
	 * \arg South pole at XYZ{ 0., 0., -bRad }
	 * \arg Point in South Atlantic Ocean at XYZ{ aRad, 0., 0. }
	 * \arg Point in Northern Indian Ocean at XYZ{ 0., aRad, 0. }
	 */
	using XYZ = std::array<double, 3u>;

	/*! \brief Type alias for geodetic coordinates (Longitude, Parallel, Alt)
	 *
	 * NOTE: units are mixed!
	 * \arg Angular units for expressing Longitude and Parallel are \b radians
	 * \arg Linear units for expressing altitude are \b meters
	 * \arg I.e. Create along the lines of:
	 * \code LPA const lpaLoc{ lonRadians, parRadians, altMeters };
	 * \endcode
	 *
	 * Uses standard sign/direction convention:
	 * \arg Longitude is positive Eastward and negative to West.
	 * \arg Parallel(Latitude) is positive Northward and negative to South.
	 * \arg (Ellipsoidal)Altitude is positive outward (up) from ellipsoid
	 *
	 * For example
	 * \arg North pole at LPA{ anyAngle,  pi/2., 0. }
	 * \arg South pole at LPA{ anyAngle, -pi/2., 0. }
	 * \arg Point in Atlantic Ocean West of Africa at LPA{ 0., 0., 0. }
	 * \arg Point in Indian Ocean at LPA{ pi/4., 0., 0. }
	 */
	using LPA = std::array<double, 3u>;

} // [peri]


#include "periDetail.h" // Inlined implementation detail


// Published interface

namespace peri
{

//! Static instances of ellipsoidal Shapes commonly used in Geodesy
namespace shape
{
	/*! \brief Defining parameters for GRS80 ellipsoid.
	 *
	 * Note that the shape of the GRS80 ellipsoid is _defined_ in terms
	 * of the equatorial axis and the second order harmonic (J2). This
	 * means that the flattening factors are _derived_ quantities that
	 * need to be computed.
	 *
	 * Ref:
	 * \arg http://geoweb.mit.edu/~tah/12.221_2005/grs80_corr.pdf
	 * , 1/f = 298.257222101 // Significant figs to about 6um at pole
	 * \arg https://iag.dgfi.tum.de/media/archives/HB2000/part4/grs80_corr.htm
	 * , 1/f = 298.257222101 // Moritz is one of the accepted definitions
	 * \arg https://geodesy.noaa.gov/library/pdfs/NOAA_Manual_NOS_NGS_0005.pdf
	 * , 1/f = 298.25722210088 // Good to about 16-digits at pole
	 * \arg https://en.wikipedia.org/wiki/Geodetic_Reference_System_1980
	 * , 1/f = 298.257 222 100 882 711 243;
	 *
	 * Flattening factor needs about 11 decimal digits after the decimal
	 * point (so about 14 overall) to provide full 'double' type precision
	 * in computed polar radius.
	 *
	 * The values in (the above copy of) report by Moritz differs from the
	 * longer precision values by about 10[nm] at the pole.
	 *
	 */
	static Shape const sGRS80
		{ Shape::fromMajorInvFlat
			( 6378137.0 // set by definition
			, 298.257222100883 // this precision provides 16-digits at pole
			)
		};

	/*! \brief Defining parameters for WGS84 ellipsoid.
	 *
	 * The WGS84 ellipsoid shape uses _both_ the equatorial radius and
	 * the (first) flattening factor to _define_ its shape.
	 *
	 * Ref:
	 * \arg
	 * ftp://ftp.nga.mil/pub2/gandg/website/wgs84/NGA.STND.0036_1.0.0_WGS84.pdf
	 * (pg 3-4): a == 6378137.0 [m], 1/f == 298.257223563 [-].
	 * \arg https://earth-info.nga.mil/GandG/publications/tr8350.2/wgs84fin.pdf
	 * (pg 3-2) a == 6378137.0 [m], 1/f == 298.257223563 [-].
	 */
	 // * 		- GM == 3.986004418e+14 [m^3/s^2]
	 // * 		- omega == 7.292115 × 10−05 [rad/s]
	static Shape const sWGS84
		{ Shape::fromMajorInvFlat
			( 6378137.0 // set by definition
			, 298.257223563 // set by definition
			)
		};

} // [peri::shape]


/*! \brief Static instances of commonly used EarthModels
 *
 * Static instances that are available to consuming code. Others
 * can be created simply by constructing a peri::Shape instance with
 * the desired size and shape values and then using that shape to
 * construct a peri::EarthModel for use with transformations.
 *
 * E.g.:
 * \code
 *	static peri::EarthModel const myEarthModel
		( peri::Shape::fromMajorInvFlat
 *			( myEquatorialRadius
 *			, myInverseFlattening
 *			)
 *		);
 * \endcode
 */
namespace model
{
	//! \brief Earth model based on GRS80 ellipsoid
	static EarthModel const GRS80(shape::sGRS80);

	//! \brief Earth model based on WGS84 ellipsoid
	static EarthModel const WGS84(shape::sWGS84);

} // [peri::model]


	/*! \brief LPA geodetic coordinates associated with cartesian location.
	 *
	 * Conversion is performed using the ellipsoid geometry provided by
	 * the earth model argument (ref peri::model to use a predefined earth
	 * model and /or create one with a custom Shape instance and provide
	 * that to an EarthModel construction)
	 *
	 * Simple Example
	 * \code
	 * peri::XYZ const locXYZ{ -1834182.999, -1131997.250, 5982812.006 };
	 * peri::LPA const gotLPA{ peri::lpaForXyz(locXYZ) }; // default WGS84
	 * \endcode
	 *
	 * Return value is combination of angle values in \b radians and
	 * altitude value in \b meters. The angles are expressed in their
	 * respective principal domains:
	 *
	 * I.e. for return value "lpa":
	 * \arg the longitude angle nominally satisfies (-pi <= lpa[0] < pi).
	 *      However the exact end conditions (of hard or soft inequalities)
	 *      depend on the specific local implementation of std::atan2().
	 * \arg the parallel(latitude) angle satisfies (-pi/2. <= lpa[1] <= pi/2.)
	 * \arg design/test region for altitude is (-1e+5 <= lpa[2] <= +1e+5).
	 *
	 * Explicit Example
	 * \code
	 * double const myEquatorRadius{ 6379000.000 };
	 * double const myInvFlattening{ 298.000000000 };
	 * peri::EarthModel const myEarthModel
	 * 	 (peri::Shape::fromMajorInvFlat(myEquatorRadius, myInvFlattening));
	 * double const xMeters{ -1834182.999 };
	 * double const yMeters{ -1131997.250 };
	 * double const zMeters{  5982812.006 };
	 * peri::XYZ const locXYZ{ xMeters, yMeters, zMeters };
	 * peri::LPA const gotLPA{ peri::lpaForXyz(locXYZ, myEarthModel) };
	 * double const & lonRadians = gotLPA[0];
	 * double const & parRadians = gotLPA[1];
	 * double const & altMeters = gotLPA[2];
	 * \endcode
	 *
	 * Singular Value Handling:
	 * For locations along the polar axis, longitude values are not
	 * defined (this is a shortcoming in the concept of geographic coordinates
	 * not a limitation of the algorithm).
	 *
	 * If a location is exactly on the polar axis, then the
	 * returned longitude (lpa[0]) is set to zero. Otherwise, the
	 * result of std::atan2() function is returned verbatim.
	 *
	 * Algorithm is designed and tested for values that are
	 * "in vicinity of" Earth surface (within +/-100km of it) - i.e. this
	 * package is designed for "terrestrial, bathymetric and aerial
	 * applications"
	 *
	 * \note If you want a conversion valid for arbitrarily large negative
	 * altitudes (i.e. approaching to center of Earth) or optimized
	 * transformations for space-based applications, then probably this is
	 * not the package to use.  You may wish to investigate a more
	 * general and complete package such as:
	 * - GeographicLib at https://sourceforge.net/projects/geographiclib
	 * - NGS https://www.ngs.noaa.gov/TOOLS/XYZ/xyz.html
	 * However, for general land surface, ocean bathymetry and/or
	 * atmospheric applications, these transformations provide high precision
	 * and accuracy (and are quite fast).
	 *
	 */
	inline
	LPA
	lpaForXyz
		( XYZ const & xyzLoc
		, EarthModel const & earthModel = model::WGS84
		)
	{
		return earthModel.lpaForXyz(xyzLoc);
	}

	/*! \brief XYZ Cartesian (ECEF) coordinates for Geodetic location.
	 *
	 * Conversion is performed using the ellipsoid geometry provided by
	 * the earth model argument (ref peri::model to use a predefined earth
	 * model and /or create one with a custom Shape instance provide that
	 * to an EarthModel construction)
	 *
	 * Simple Example
	 * \code
	 * peri::LPA const locLPA{ 3.123456, 1.56789, 123.456 };
	 * peri::XYZ const gotXYZ{ peri::xyzForLpa(locLPA) }; // default WGS84
	 * \endcode
	 *
	 * Input longitude angle, lpaLoc[0], and parallel(latitude)
	 * angle, lpaLoc[1], must be expressed in \b radians while the
	 * altitude (lpaLoc[2]) must be expressed in (ISO standard) \b meters.
	 *
	 * Note that the LPA results are well defined for any numerically
	 * valid input values.
	 *
	 * The longitude (lpaLoc[0]) and parallel (aka latitude, lpa[1])
	 * angles are wrapped into principal angle intervals for interpretation.
	 * To be within the design/test domain, the ellipsoid altitude
	 * values (lpaLoc[2]) should be within about 100[km] of earth
	 * surface (e.g. -1e+5<lpaLoc[2]<+1e+5).
	 *
	 * Explicit Example
	 * \code
	 * double const myEquatorRadius{ 6379000.000 };
	 * double const myInvFlattening{ 298.000000000 };
	 * peri::EarthModel const myEarthModel
	 * 	 (peri::Shape::fromMajorInvFlat(myEquatorRadius, myInvFlattening));
	 * double const lonRadians{ 3.123456 };
	 * double const parRadians{ 1.56789 };
	 * double const altMeters{ 12345.678 };
	 * peri::LPA const locLPA{ lonRadians, parRadians, altMeters };
	 * peri::XYZ const gotXYZ{ peri::xyzForLpa(locLPA, myEarthModel) };
	 * double const & xECEF = gotXYZ[0];
	 * double const & yECEF = gotXYZ[1];
	 * double const & zECEF = gotXYZ[2];
	 * \endcode
	 *
	 */
	inline
	XYZ
	xyzForLpa
		( LPA const & lpaLoc
		, EarthModel const & earthModel = model::WGS84
		)
	{
		return earthModel.xyzForLpa(lpaLoc);
	}

} // [peri]


#endif // peridetic_INCL_

