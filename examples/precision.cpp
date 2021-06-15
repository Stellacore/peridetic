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

/*! \file Simple code spike to display data type precision at Earth radius.
 *
 * Display (to std::out) impact of floating point data type resolutions
 * in context of Earth radius.
 */

#include <cmath>
#include <iomanip>
#include <iostream>
#include <limits>
#include <sstream>


//! Overall scaling factors - assumes ISO units (e.g. [m]).
namespace scale
{
	//! Nominal earth radius (as float)
	template <typename Type>
	constexpr
	Type
	earthRadius
		()
	{
		constexpr std::uint32_t eRadMeters{ 6370000u };
		return static_cast<Type>(eRadMeters);
	}

} // [scale]

// useful formatting options
namespace info
{
	//! Size for title padding
	constexpr std::size_t sTitlePad{ 24u };

	//! Number formatting/display widths
	namespace wide
	{
		constexpr std::size_t sSign{ 1u };
		constexpr std::size_t sInt{ 9u };
		constexpr std::size_t sDot{ 1u };
		constexpr std::size_t sFrac{ 16u };
		constexpr std::size_t sSize{ sSign + sInt + sDot + sFrac };

	} // [wide]

	//! Full precision display of double value
	std::string
	asDouble
		( double const & dubValue
		, std::string const & title
		)
	{
		std::ostringstream oss;
		oss
			<< std::setw(sTitlePad) << title
			<< " " << std::setw(wide::sSize)
			<< std::fixed << std::setprecision(16u)
			<< dubValue
			;
		return oss.str();
	}

	//! Cast anyValue to double, multiply by sphere radius and display
	template <typename Type>
	std::string
	atSurface
		( Type const & anyValue
		, std::string const & title
		, double const & radValue = 6370000.
		)
	{
		double const dubValue{ static_cast<double>(anyValue) };
		double const surfValue{ dubValue * radValue };
		return asDouble(surfValue, title);
	}

	//! Cast anyValue to double and display
	template <typename Type>
	std::string
	general
		( Type const & anyValue
		, std::string const & title
		)
	{
		std::ostringstream oss;
		double const dubValue{ static_cast<double>(anyValue) };
		oss
			<< asDouble(dubValue, "relativeTo1:" + title)
			<< " "
			<< atSurface(dubValue, "atEarthRad:" + title)
			;
		return oss.str();
	}

} // [info]


int
main
	()
{
	// Evaluate for double precision (e.g. 53-bit mantissa for 64-bit type)
	constexpr double dEps{ std::numeric_limits<double>::epsilon() };
	constexpr double dAt{ scale::earthRadius<double>() * (1. + dEps) };

	// Evaluate for float precision (e.g. 24-bit mantissa for 32-bit type)
	constexpr double fEps{ std::numeric_limits<float>::epsilon() };
	constexpr double fAt{ scale::earthRadius<float>() * (1. + fEps) };

	std::cout << '\n';
	std::cout << "Example Values:" << '\n';
	std::cout << info::asDouble(dAt, "dAt") << '\n';
	std::cout << info::asDouble(fAt, "fAt") << '\n';
	std::cout << "Smallest Increments:" << '\n';
	std::cout << info::general(dEps, "dEps") << '\n';
	std::cout << info::general(fEps, "fEps") << '\n';
	std::cout << std::endl;

	return 0;
}

