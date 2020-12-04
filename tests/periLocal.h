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


#ifndef peri_Local_INCL_
#define peri_Local_INCL_


#include "peridetic.h"

#include <cmath>
#include <cstdint>
#include <iomanip>
#include <iostream>
#include <limits>
#include <sstream>

/* -- consider adding

	* TODO - principalValueForLpa() // put angles in principal domain
	* TODO - isInOptimalDomainLpa() // within +/-100[km] of ellipsoid surface
	* TODO - isInOptimalDomainXyz() // expensive, does conversion to check
	* TODO - isFeasibleLpa() // checks negative altitude vs minimum radius
	* TODO - isValid() // not a nan

*/

namespace peri
{
	/*! \brief Angular value small enough to be negligible near Earth surface.
	 *
	 * Nominally, this value is small enough so that, when multiplied by
	 * an Earth radius, the associated linear arc is small (comparable to
	 * #sSmallLinear value.
	 *
	 * The value { 1. /1024./1024./1024./1024./1024. } is:
	 * \arg (~= 9e-16 [rad] < 1[fRad]) - as pure angle
	 * \arg (!= 5.7e-9 [m] < 6[nm]) - linear displacement near Earth Surface
	 */
	constexpr double sSmallAngular{ 1. /1024./1024./1024./1024./1024. };

	/*! \brief Linear value small enough to be negligible for applications.
	 *
	 * The value { 1. /1024./1024./128. } is:
	 * \arg (~= 7.45e-9 < 7.5[nm])
	 */
	constexpr double sSmallLinear{ 1. /1024./1024./128. };

	//! Nominal ellipsoid altitude +/- limit for design/testing verification.
	constexpr double sDomainAltDelta{ 100000. }; // ~equator+100[km]

	//! Nominal upper domain limit for design/testing verification.
	// ~equator+100[km]
	constexpr double sDomainRadiusMax{ 6380000. + sDomainAltDelta };

	//! Nominal lower domain limit for design/testing verification.
	// ~pole-100[km]
	constexpr double sDomainRadiusMin{ 6350000. - sDomainAltDelta };

	//
	// data validity operations
	//

	//! True if data value is not a sNan
	inline
	bool
	isValid
		( double const & value
		)
	{
		return (! std::isnan(value));
	}

	//! True if all three components are isValid()
	inline
	bool
	isValid
		( std::array<double, 3u> const & array
		)
	{
		return (isValid(array[0]) && isValid(array[1]) && isValid(array[2]));
	}

	//
	// math operations
	//

	//! PI computed to full precision
	inline
	double
	pi
		()
	{
		return (4. * std::atan(1.));
	}

	//! Principal value of anyAngle such that (-pi <= principalValue() < pi)
	inline
	double
	principalAngle
		( double const & anyAngle
		)
	{
		double angle{ anyAngle };
		double const maxAngle{ pi() };
		double const minAngle{ - maxAngle };
		if ((anyAngle < minAngle) || (! (anyAngle < maxAngle)))
		{
			using namespace std;
			angle = atan2(sin(anyAngle), cos(anyAngle));
			// some libraries have atan2 sloppy about +/- pi
			if (maxAngle == angle)
			{
				angle = -minAngle;
			}
		}
		return angle;
	}

	//! Quadrature angle (pi/2) 
	inline
	double
	qtrTurn
		()
	{
		return (.5 * pi());
	}

	//! Radian equivalent of value expressed in degrees
	inline
	double
	radForDeg
		( double const & valueDeg
		)
	{
		return ((pi()/180.) * valueDeg);
	}

	//! Degree equivalent of radian value
	inline
	double
	degForRad
		( double const & valueRad
		)
	{
		return ((180./pi()) * valueRad);
	}

	//! True if value is valid and in half open interval (min <= value < max)
	inline
	bool
	isValidWithin
		( double const & value
		, double const & minIncluded
		, double const & maxExcluded
		)
	{
		return
			{  isValid(value)
			&& (! (minIncluded < value))
			&& (value < maxExcluded)
			};
	}

	//! Radians magnitude for DMS specification with all DMS fields positive
	inline
	double
	radMagFromDMS
		( double const & deg
		, double const & min
		, double const & sec
		)
	{
		double radMag{ sNan };
		double const degMag{ std::abs(deg) };
		// check for Babylonian consistency and get down with sexagesimalism
		if (  isValidWithin(degMag, 0., 180.)
		   && isValidWithin(min, 0.,  60.)
		   && isValidWithin(sec, 0.,  60.)
		   )
		{
			double const degMag{ deg + (min + sec/60.)/60. };
			radMag = radForDeg(degMag);
		}
		// return with same algebraic sign as on deg input
		return std::copysign(radMag, deg);
	}

	//
	// numeric operations
	//

	//! True if value[AB] difference is absolutely less than tolerance
	inline
	bool
	sameEnough
		( double const & valueA
		, double const & valueB
		, double const & tolAbs
		)
	{
		return (std::abs(valueA - valueB) < tolAbs);
	}

} // [peri]

//! String encoding useful for geodetic data values
namespace peri::string
{

	//! Value encoded as fixed format string suitable for Cartesian values
	inline
	std::string
	fixedLinear
		( double const & value
		, std::string const & title = {}
		)
	{
		constexpr std::size_t numDigits{ 3u }; //!< After decimal point
		std::ostringstream oss;
		if (! title.empty())
		{
			oss << std::setw(16u) << title << " ";
		}
		oss << std::fixed
			<< std::setw(9u+numDigits) << std::setprecision(numDigits)
			<< value;
		return oss.str();
	}

	//! Value encoded as fixed format string suitable for Geodetic angle values
	inline
	std::string
	fixedAngular
		( double const & value
		, std::string const & title = {}
		)
	{
		constexpr std::size_t numDigits{ 10u }; //!< After decimal point
		std::ostringstream oss;
		if (! title.empty())
		{
			oss << std::setw(16u) << title << " ";
		}
		oss << std::fixed
			<< std::setw(3u+numDigits) << std::setprecision(numDigits)
			<< value;
		return oss.str();
	}

	//! Value encoded as maximum (8-byte) precision e-notation
	inline
	std::string
	allDigits
		( double const & value
		, std::string const & title = {}
		)
	{
		std::ostringstream oss;
		if (! title.empty())
		{
			oss << std::setw(16u) << title << " ";
		}
		oss << std::scientific
			<< std::setw(23u) << std::setprecision(15u) << value;
		return oss.str();
	}

	//! Values in e-notation at full precision (8-byte IEEE-754)
	inline
	std::string
	allDigits
		( std::array<double, 3u> const & values
		, std::string const & title = {}
		, std::size_t const tPad = 16u //!< space to put title
		)
	{
		std::ostringstream oss;
		if (! title.empty())
		{
			oss << std::setw(tPad) << title << " ";
		}
		oss << string::allDigits(values[0])
			<< " "
			<< string::allDigits(values[1])
			<< " "
			<< string::allDigits(values[2])
			;
		return oss.str();
	}

} // [peri::string]

//! Special functions for Cartesian coordinates
namespace peri::xyz
{
	constexpr peri::XYZ sNull{ peri::sNull };

	//! True if all components of xyz[AB] agree closer than tolerance
	inline
	bool
	sameEnough
		( XYZ const & xyzA
		, XYZ const & xyzB
		, double const & tolLinear = sSmallLinear
		)
	{
		return
			(  peri::sameEnough(xyzA[0], xyzB[0], tolLinear)
			&& peri::sameEnough(xyzA[1], xyzB[1], tolLinear)
			&& peri::sameEnough(xyzA[2], xyzB[2], tolLinear)
			);
	}

	//! String encoding of Cartesian value with proper significant precision
	inline
	std::string
	infoString
		( XYZ const & xyz
		, std::string const & title = {}
		, std::size_t const tPad = 16u //!< space to put title
		)
	{
		std::ostringstream oss;
		if (! title.empty())
		{
			oss << std::setw(tPad) << title << " ";
		}
		oss << string::fixedLinear(xyz[0])
			<< " "
			<< string::fixedLinear(xyz[1])
			<< " "
			<< string::fixedLinear(xyz[2])
			;
		return oss.str();
	}

} // [peri::xyz]

//! Special functions for Geodetic coordinates
namespace peri::lpa
{
	constexpr peri::LPA sNull{ peri::sNull };

	//! True if all components agree closer than design tolerance
	inline
	bool
	sameEnough
		( LPA const & lpaA
		, LPA const & lpaB
		, double const & tolAngular = sSmallAngular
		, double const & tolLinear = sSmallLinear
		)
	{
		return
			(  peri::sameEnough(lpaA[0], lpaB[0], tolAngular)
			&& peri::sameEnough(lpaA[1], lpaB[1], tolAngular)
			&& peri::sameEnough(lpaA[2], lpaB[2], tolLinear)
			);
	}

	//! String encoding of geodetic value with proper significant precision
	inline
	std::string
	infoString
		( LPA const & lpa
		, std::string const & title = {}
		, std::size_t const tPad = 16u //!< space to put title
		)
	{
		std::ostringstream oss;
		if (! title.empty())
		{
			oss << std::setw(tPad) << title << " ";
		}
		oss << string::fixedAngular(lpa[0])
			<< " "
			<< string::fixedAngular(lpa[1])
			<< " "
			<< string::fixedLinear(lpa[2])
			;
		return oss.str();
	}

} // [peri::lpa]

#endif // peri_Local_INCL_
