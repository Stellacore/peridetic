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

#include "corsDataPairs.h"
#include "corsDataParser.h"

#include <iostream>
#include <vector>



namespace
{
	//! Check XYZ from LPA transformations.
	int
	test1a
		()
	{
		int errCount{ 0 };

		// CORS data based on the GRS80 ellipsoid
		peri::EarthModel const & earth = peri::model::GRS80;

		for (std::string const & staText : peri::cors::sStationTexts)
		{
			using peri::cors::DataParser;
			DataParser const parser{ DataParser::from(staText) };
			// access data elements
			peri::LPA const & locLPA = parser.theLPA;
			peri::XYZ const & expXYZ = parser.theXYZ;
			// evaluate transform
			peri::XYZ const gotXYZ{ peri::xyzForLpa(locLPA, earth) };
			// check results
		//	double const & tolLin = peri::sSmallLinear;
			double const tolLin{ 1./1024. }; // CORS files only good to [mm]
			if (! peri::xyz::sameEnough(gotXYZ, expXYZ, tolLin))
			{
				using peri::operator-;
				peri::XYZ const difXYZ{ gotXYZ - expXYZ };
				std::cerr << "Failure of CORS XYZ from LPA test" << '\n';
				std::cerr << peri::lpa::infoString(locLPA, "locLPA") << '\n';
				std::cerr << peri::xyz::infoString(expXYZ, "expXYZ") << '\n';
				std::cerr << peri::xyz::infoString(gotXYZ, "gotXYZ") << '\n';
				std::cerr << peri::string::allDigits(expXYZ, "expXYZ") << '\n';
				std::cerr << peri::string::allDigits(gotXYZ, "gotXYZ") << '\n';
				std::cerr << peri::string::allDigits(difXYZ, "difXYZ") << '\n';
				std::cerr << peri::string::allDigits(tolLin, "tolLin") << '\n';
				++errCount;
				break;
			}
		}
		return errCount;
	}

	//! Check LPA from XYZ transformations.
	int
	test1b
		()
	{
		int errCount{ 0 };

		// CORS data based on the GRS80 ellipsoid
		peri::EarthModel const & earth = peri::model::GRS80;

		for (std::string const & staText : peri::cors::sStationTexts)
		{
			using peri::cors::DataParser;
			DataParser const parser{ DataParser::from(staText) };
			// access data elements
			peri::XYZ const & locXYZ = parser.theXYZ;
			peri::LPA const & expLPA = parser.theLPA;

			// evaluate transform
			peri::LPA const gotLPA{ peri::lpaForXyz(locXYZ, earth) };
			// check results
			// NOTE: relax test conditions here to accommodate limited
			// precision in NOAA/CORS published values.
			constexpr double tolAng{ 172. / 1024./1024./1024./1024. };
			constexpr double tolLin{ 1./1024. }; // CORS files only good to [mm]
			if (! peri::lpa::sameEnough(gotLPA, expLPA, tolAng, tolLin))
			{
				using peri::operator-;
				peri::LPA const difLPA{ gotLPA - expLPA };
				std::cerr << "Failure of CORS LPA from XYZ test" << '\n';
				std::cerr << peri::xyz::infoString(locXYZ, "locXYZ") << '\n';
				std::cerr << peri::lpa::infoString(expLPA, "expLPA") << '\n';
				std::cerr << peri::lpa::infoString(gotLPA, "gotLPA") << '\n';
				std::cerr << peri::string::allDigits(expLPA, "expLPA") << '\n';
				std::cerr << peri::string::allDigits(gotLPA, "gotLPA") << '\n';
				std::cerr << peri::string::allDigits(difLPA, "difLPA") << '\n';
				std::cerr << peri::string::allDigits(tolAng, "tolAng") << '\n';
				std::cerr << peri::string::allDigits(tolLin, "tolLin") << '\n';
				++errCount;
			//	break;
			}
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
	errCount += test1a(); // Cartesian from Geographic w/ CORS examples
	errCount += test1b(); // Geographic from Cartesian w/ CORS examples
	return errCount;
}
