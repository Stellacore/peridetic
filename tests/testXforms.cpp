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

#include "periLocal.h"
#include "periSim.h"

#include <iostream>
#include <vector>


namespace
{
	//! True if *all* components of values are inValid
	bool
	isAllNan
		( std::array<double, 3u> const & values
		)
	{
		return
			{  (! peri::isValid(values[0]))
			&& (! peri::isValid(values[1]))
			&& (! peri::isValid(values[2]))
			};
	}

	//! Check null value handling
	int
	test0
		()
	{
		int errCount{ 0 };
		using namespace peri;
		double const & nan = sNan;
		std::array<LPA, 7u> const badLPAs
			{ LPA{ nan,  0.,  0. }
			, LPA{  0., nan,  0. }
			, LPA{  0.,  0., nan }
			, LPA{  0., nan, nan }
			, LPA{ nan,  0., nan }
			, LPA{ nan, nan,  0. }
			, LPA{ nan, nan, nan }
			};
		std::array<XYZ, 7u> const badXYZs
			{ XYZ{ nan,  0.,  0. }
			, XYZ{  0., nan,  0. }
			, XYZ{  0.,  0., nan }
			, XYZ{  0., nan, nan }
			, XYZ{ nan,  0., nan }
			, XYZ{ nan, nan,  0. }
			, XYZ{ nan, nan, nan }
			};
		for (std::size_t nn{0u} ; nn < 7u ; ++nn)
		{
			LPA const & badLPA = badLPAs[nn];
			XYZ const gotXYZ{ xyzForLpa(badLPA) };
			if (! isAllNan(gotXYZ))
			{
				std::cerr << "Failure of nullXYZ test" << '\n';
				std::cerr << lpa::infoString(badLPA, "badLPA") << '\n';
				std::cerr << xyz::infoString(gotXYZ, "gotXYZ") << '\n';
				++errCount;
			}
			LPA const & badXYZ = badXYZs[nn];
			LPA const gotLPA{ lpaForXyz(badXYZ) };
			if (! isAllNan(gotLPA))
			{
				std::cerr << "Failure of nullLPA test" << '\n';
				std::cerr << lpa::infoString(badXYZ, "badXYZ") << '\n';
				std::cerr << xyz::infoString(gotLPA, "gotLPA") << '\n';
				++errCount;
			}
		}
		return errCount;
	}

	//! Functor for conducting round-trip transformation consistency test
	struct RoundTripper
	{
		peri::EarthModel const earth{ peri::model::WGS84 };
		double const theTolAng{ peri::sSmallAngular };
		double const theTolLin{ peri::sSmallLinear };

		//! Conduct round-trip test at requested LPA location
		int
		errCountAt
			( peri::XYZ const & expLPA
			) const
		{
			int errCount{ 0 };
			peri::XYZ const gotXYZ{ peri::xyzForLpa(expLPA, earth) };
			peri::XYZ const gotLPA{ peri::lpaForXyz(gotXYZ, earth) };
			peri::XYZ const chkXYZ{ peri::xyzForLpa(gotLPA, earth) };
			using peri::string::allDigits;
			if (! peri::lpa::sameEnough(gotLPA, expLPA, theTolAng, theTolLin))
			{
				std::cerr << "Failure of round-trip LPA test" << std::endl;
				using peri::operator-;
				peri::LPA const difLPA{ gotLPA - expLPA };
				std::cerr << peri::lpa::infoString(expLPA, "expLPA") << '\n';
				std::cerr << peri::lpa::infoString(gotLPA, "gotLPA") << '\n';
				std::cerr << allDigits(expLPA, "expLPA") << '\n';
				std::cerr << allDigits(gotLPA, "gotLPA") << '\n';
				std::cerr << allDigits(difLPA, "difLPA") << '\n';
				std::cerr << allDigits(theTolAng, "tolAng") << '\n';
				std::cerr << allDigits(theTolLin, "tolLin") << '\n';
				++errCount;
			}
			if (! peri::xyz::sameEnough(gotXYZ, chkXYZ, theTolLin))
			{
				std::cerr << "Failure of round-trip XYZ test" << std::endl;
				using peri::operator-;
				peri::XYZ const difXYZ{ gotXYZ - chkXYZ };
				std::cerr << peri::lpa::infoString(expLPA, "expLPA") << '\n';
				std::cerr << peri::xyz::infoString(gotXYZ, "gotXYZ") << '\n';
				std::cerr << peri::xyz::infoString(chkXYZ, "chkXYZ") << '\n';
				std::cerr << allDigits(gotXYZ, "gotXYZ") << '\n';
				std::cerr << allDigits(chkXYZ, "chkXYZ") << '\n';
				std::cerr << allDigits(difXYZ, "difXYZ") << '\n';
				std::cerr << allDigits(theTolLin, "tolLin") << '\n';
				++errCount;
			}
			return errCount;
		}

	}; // RoundTripper

	//! Check round-trip consistency between forward/inverse transformations
	int
	test2
		()
	{
		int errCount{ 0 };

		// commonly used Earth model
		peri::EarthModel const & earth = peri::model::WGS84;

		// collection of values spanning design domain
		constexpr std::size_t numLon{  53u };
		constexpr std::size_t numPar{  67u };
		constexpr std::size_t numAlt{  73u };
		std::vector<peri::LPA> const expLPAs
			{ peri::sim::bulkSamplesLpa(numLon, numPar, numAlt) };

		// test at full precision
		RoundTripper const rt
			{ earth, peri::sSmallAngular, peri::sSmallLinear };
		for (peri::LPA const & expLPA : expLPAs)
		{
			errCount += rt.errCountAt(expLPA);
		}

		return errCount;
	}


	/*! \brief Check round-trip consistency at GNSS satellite altitudes.
	 *
	 * Experimental evaluation of round-trip precision at large altitudes
	 * comparable to GNSS satellite system orbits.
	 *
	 * \note Ref internal code tolerance values used for converging testing.
	 *
	 */
	int
	test3a
		()
	{
		int errCount{ 0 };

		// commonly used Earth model
		peri::EarthModel const & earth = peri::model::GRS80;

		constexpr std::size_t numLon{  19u };
		constexpr std::size_t numPar{  37u };
		std::vector<double> const lonSamps{ peri::sim::bulkSamplesLon(numLon) };
		std::vector<double> const parSamps{ peri::sim::bulkSamplesPar(numPar) };
		std::vector<double> const altSamps
			{ -100.e+3 // lower design bound
			,  100.e+3 // upper design bound
			,   11.e+6 // approx upper limit on design precision
			,   25.e+6 // nominal GNSS altitudes // tol ~= 1.e-8
			};
		std::vector<peri::LPA> const expLPAs
			{ peri::sim::comboSamplesLpa(lonSamps, parSamps, altSamps) };

		// test with reduced precision threshold for the large distances
		constexpr double tolLinGNSS{ 2.e-8 };
		constexpr double tolAngGNSS{ peri::sSmallAngular };
		RoundTripper const rtGNSS{ earth, tolAngGNSS, tolLinGNSS };
		for (peri::LPA const & expLPA : expLPAs)
		{
			errCount += rtGNSS.errCountAt(expLPA);
		}

		return errCount;
	}

	/*! \brief Check round-trip consistency at Lunar altitudes.
	 *
	 * Experimental evaluation of round-trip precision at very large
	 * altitudes comparable to Lunar distances.
	 *
	 * \note Ref internal code tolerance values used for converging testing.
	 *
	 */
	int
	test3b
		()
	{
		int errCount{ 0 };

		// commonly used Earth model
		peri::EarthModel const & earth = peri::model::GRS80;

		constexpr std::size_t numLon{ 17u };
		constexpr std::size_t numPar{ 29u };
		std::vector<double> const lonSamps{ peri::sim::bulkSamplesLon(numLon) };
		std::vector<double> const parSamps{ peri::sim::bulkSamplesPar(numPar) };
		std::vector<double> const altSamps
			{ -100.e+3 // lower design bound
			,  100.e+3 // upper design bound
			,  410.e+6 // slightly beyond lunar apogee
			};
		std::vector<peri::LPA> const expLPAs
			{ peri::sim::comboSamplesLpa(lonSamps, parSamps, altSamps) };

		// test with reduced precision threshold for the large distances
		constexpr double tolLinGNSS{ 2.e-7 };
		constexpr double tolAngGNSS{ peri::sSmallAngular };
		RoundTripper const rtGNSS{ earth, tolAngGNSS, tolLinGNSS };
		for (peri::LPA const & expLPA : expLPAs)
		{
			errCount += rtGNSS.errCountAt(expLPA);
		}

		return errCount;
	}

	/*! \brief Check round-trip consistency at Earth interior altitudes.
	 *
	 * Experimental evaluation of round-trip precision inside Earth.
	 *
	 * \note Ref internal code tolerance values used for converging testing.
	 *
	 */
	int
	test3c
		()
	{
		int errCount{ 0 };

		// commonly used Earth model
		peri::EarthModel const & earth = peri::model::GRS80;

		constexpr std::size_t numLon{ 59u };
		constexpr std::size_t numPar{ 29u };
		std::vector<double> const lonSamps{ peri::sim::bulkSamplesLon(numLon) };
		std::vector<double> const parSamps{ peri::sim::bulkSamplesPar(numPar) };
		std::vector<double> const altSamps
			{ -100.e+3 // lower design bound
			,  100.e+3 // upper design bound
			, -3200.e+3 // about half way to center
			, -5800.e+3 // about 9% of an Earth radius
		//	, -6300.e+3 // approx limit on convergence
			};
		std::vector<peri::LPA> const expLPAs
			{ peri::sim::comboSamplesLpa(lonSamps, parSamps, altSamps) };

		// test with reduced precision threshold for the large distances
		constexpr double tolLinGNSS{ peri::sSmallLinear };
		constexpr double tolAngGNSS{ peri::sSmallAngular };
		RoundTripper const rtGNSS{ earth, tolAngGNSS, tolLinGNSS };
		for (peri::LPA const & expLPA : expLPAs)
		{
			errCount += rtGNSS.errCountAt(expLPA);
		}

		return errCount;
	}

	/*! \brief Evaluation playground for experimentation.
	 *
	 */
	void
	testLimits
		()
	{
		int errCount{ 0u };

		// commonly used Earth model
		peri::EarthModel const & earth = peri::model::GRS80;

		// trials along a simple profile
		constexpr std::size_t numLon{ 1u };
		constexpr std::size_t numPar{ 19u };
		std::vector<double> const lonSamps{ peri::sim::bulkSamplesLon(numLon) };
		std::vector<double> const parSamps{ peri::sim::bulkSamplesPar(numPar) };

		constexpr double lyDist{ 26700. };
		constexpr double mPerSec{ 3.e8 };
		constexpr double secPerYr{ 3600.*24.*365.25 };
		constexpr double mPerLy{ mPerSec * secPerYr };
		constexpr double mDist{ mPerLy * lyDist };

		std::vector<double> const altSamps
			{ -100.e+3 // lower design bound
			,  100.e+3 // upper design bound
			,  100.e+6 // upper design bound
			,  mDist // dist to center of galaxy
			};
		constexpr double tolLinGNSS{ 100.e+3 }; // most affected
		constexpr double tolAngGNSS{ peri::sSmallAngular };

		std::vector<peri::LPA> const expLPAs
			{ peri::sim::comboSamplesLpa(lonSamps, parSamps, altSamps) };

		// test with reduced precision threshold for the large distances
		RoundTripper const rtGNSS{ earth, tolAngGNSS, tolLinGNSS };
		for (peri::LPA const & expLPA : expLPAs)
		{
			errCount += rtGNSS.errCountAt(expLPA);
		}

		if (0 < errCount)
		{
			std::cout << "testLimits: errCount: " << errCount << std::endl;
		}
	}
}


//! Check Peridetic transformations
int
main
	()
{
	int errCount{ 0 };
	errCount += test0(); // Null values
	errCount += test2(); // RoundTrip consistency testing
	errCount += test3a(); // RoundTrip evaluation in near outer space
	errCount += test3b(); // RoundTrip evaluation in far outer space
	errCount += test3c(); // RoundTrip evaluation interior to Earth
	testLimits(); // For developer experimentation
	return errCount;
}
