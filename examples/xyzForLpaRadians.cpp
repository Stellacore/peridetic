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

#include <cmath>
#include <cstdlib>
#include <iomanip>
#include <iostream>
#include <string>


//! xyzForLpa program: Compute/Report XYZ values given LPA arguments.
int
main
	( int argc
	, char ** argv
	)
{
	if (! (3 < argc))
	{
		std::cerr
			<< "Usage: <progname> lonInRadians parInRadians altInMeters"
			<< "\nlonInRadians: longitude angle in radian units"
			<< "\nparInRadians: parallel (latitude) angle in radian units"
			<< "\naltInMeters: altitude (ellipsoidal) in units of meters"
			<< '\n';
	}
	else
	{
		// Example: Important part

		// note: atof values are 0 if decoding error, but not the point here
		peri::LPA const locLPA
			{ std::atof(argv[1])
			, std::atof(argv[2])
			, std::atof(argv[3])
			};
		peri::XYZ const locXYZ{ peri::xyzForLpa(locLPA) };

		// Echo input values and report converted ones

		// For useful formatting functions, ref. .../tests/periLocal.h
		std::cout << std::fixed ;
		using std::setw;
		std::cout << "input: locLPA: "
			<< std::setprecision(9)
			<< " " << setw(12) << locLPA[0]
			<< " " << setw(12) << locLPA[1]
			<< std::setprecision(3)
			<< " " << setw(12) << locLPA[2]
			<< " " << "[rad,rad,m]"
			<< std::endl;
		std::cout << "equiv: locXYZ: "
			<< std::setprecision(3)
			<< " " << setw(12) << locXYZ[0]
			<< " " << setw(12) << locXYZ[1]
			<< " " << setw(12) << locXYZ[2]
			<< " " << "[m,m,m]"
			<< std::endl;
	}
	return 0;
}

