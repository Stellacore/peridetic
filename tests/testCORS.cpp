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

#include "peridetic.h"

#include "corsDataParser.h"
#include "corsDataPairs.h"

#include "periLocal.h"

#include <iostream>


namespace
{
	//! Check decoding of CORS formatted data strings (cors::DataParser)
	int
	test0
		()
	{
		int errCount{ 0 };

		// same syntax used in cors::CoordPairs.h
		static std::string const sTextTest
		{
		" ITRF2014 POSITION (EPOCH 2010.0)   "
		" Faked in Nov 2020 using madeup data CORS reporting style.  "
		"     X =   1000000.001 m     latitude    =  01 23 45.00006 S  "
		"     Y =   2000000.002 m     longitude   = 076 54 32.00001 E  "
		"     Z =   3000000.003 m     ellipsoid height = 9876.543   m  "
		};

		// expected values
		constexpr double expX{ 1000000.001 };
		constexpr double expY{ 2000000.002 };
		constexpr double expZ{ 3000000.003 };
		using peri::radForDeg;
		double const expP{ -radForDeg( 1. + 23./60. + 45.00006/60./60.) };
		double const expL{  radForDeg(76. + 54./60. + 32.00001/60./60.) };
		constexpr double expA{ 9876.543 };
		//
		peri::XYZ const expXYZ{ expX, expY, expZ };
		peri::LPA const expLPA{ expL, expP, expA };

		peri::cors::DataParser const parser
			{ peri::cors::DataParser::from(sTextTest) };

		peri::XYZ const & gotXYZ = parser.theXYZ;
		peri::LPA const & gotLPA = parser.theLPA;

		// check Cartesian values
		if (! peri::xyz::sameEnough(gotXYZ, expXYZ))
		{
			std::cerr << "Failure of XYZ parse test" << '\n';
			std::cerr << peri::xyz::infoString(expXYZ, "expXYZ") << '\n';
			std::cerr << peri::xyz::infoString(gotXYZ, "gotXYZ") << '\n';
			++errCount;
		}

		// check geodetic values
		if (! peri::lpa::sameEnough(gotLPA, expLPA))
		{
			std::cerr << "Failure of LPA parse test" << '\n';
			std::cerr << peri::lpa::infoString(expLPA, "expLPA") << '\n';
			std::cerr << peri::lpa::infoString(gotLPA, "gotLPA") << '\n';
			++errCount;
		}

		return errCount;
	}

}


//! Check Peridetic transformations
int
main
	()
{
	int errCount{ 0 };
	errCount += test0();
	// errCount += test1();
	return errCount;
}

