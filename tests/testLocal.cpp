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


namespace
{
	//! Check basic XYZ coordinate operations
	int
	test0
		()
	{
		int errCount{ 0 };

		// construct values that are slightly same and different
		constexpr double mag{ peri::sDomainRadiusMax };
		constexpr double del{ peri::sSmallLinear };
		constexpr double eps{ .5 * peri::sSmallLinear };
		peri::XYZ const xyz{ mag+1000., -(mag+2000.), mag+3000. };
		peri::XYZ const xyzDiff{ xyz[0]+del, xyz[1]+del, xyz[2]+del };
		peri::XYZ const xyzSame{ xyz[0]+eps, xyz[1]+eps, xyz[2]+eps };

		// check values that should test as same
		if (! peri::xyz::sameEnough(xyzSame, xyz)) // sameEnough should be true
		{
			using peri::operator-;
			peri::XYZ const delta{ xyzSame - xyz };
			std::cerr << "Failure of XYZ Same test" << std::endl;
			std::cerr << peri::xyz::infoString(xyz, "xyz") << std::endl;
			std::cerr << peri::xyz::infoString(xyzSame, "xyzSame") << std::endl;
			std::cerr << peri::xyz::infoString(delta, "delta") << std::endl;
			std::cerr << peri::string::allDigits(delta, "delta") << std::endl;
			++errCount;
		}

		// check values that should test as different
		if (  peri::xyz::sameEnough(xyzDiff, xyz)) // sameEnough should be false
		{
			using peri::operator-;
			peri::XYZ const delta{ xyzDiff - xyz };
			std::cerr << "Failure of XYZ Diff test" << std::endl;
			std::cerr << peri::xyz::infoString(xyz, "xyz") << std::endl;
			std::cerr << peri::xyz::infoString(xyzDiff, "xyzDiff") << std::endl;
			std::cerr << peri::xyz::infoString(delta, "delta") << std::endl;
			std::cerr << peri::string::allDigits(delta, "delta") << std::endl;
			++errCount;
		}

		return errCount;
	}

	//! Check basic LPA coordinate operations
	int
	test1
		()
	{
		int errCount{ 0 };

		// construct values that are slightly same and different
		double const pi{ peri::pi() };
		constexpr double delLP{ peri::sSmallAngular };
		constexpr double delA{ peri::sSmallLinear };
		constexpr double epsLP{ .5 * peri::sSmallAngular };
		constexpr double epsA{ .5 * peri::sSmallLinear };
		peri::LPA const lpa{ pi, -pi, 123000. };
		peri::LPA const lpaDiff{ lpa[0]+delLP, lpa[1]+delLP, lpa[2]+delA };
		peri::LPA const lpaSame{ lpa[0]+epsLP, lpa[1]+epsLP, lpa[2]+epsA };

		// check values that should test as same
		if (! peri::lpa::sameEnough(lpaSame, lpa)) // sameEnough should be true
		{
			using peri::operator-;
			peri::LPA const delta{ lpaSame - lpa };
			std::cerr << "Failure of LPA Same test" << std::endl;
			std::cerr << peri::lpa::infoString(lpa, "lpa") << std::endl;
			std::cerr << peri::lpa::infoString(lpaSame, "lpaSame") << std::endl;
			std::cerr << peri::lpa::infoString(delta, "delta") << std::endl;
			std::cerr << peri::string::allDigits(delta, "delta") << std::endl;
			++errCount;
		}

		// check values that should test as different
		if (  peri::lpa::sameEnough(lpaDiff, lpa)) // sameEnough should be false
		{
			using peri::operator-;
			peri::LPA const delta{ lpaDiff - lpa };
			std::cerr << "Failure of LPA Diff test" << std::endl;
			std::cerr << peri::lpa::infoString(lpa, "lpa") << std::endl;
			std::cerr << peri::lpa::infoString(lpaDiff, "lpaDiff") << std::endl;
			std::cerr << peri::lpa::infoString(delta, "delta") << std::endl;
			std::cerr << peri::string::allDigits(delta, "delta") << std::endl;
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
	errCount += test1();
	return errCount;
}

