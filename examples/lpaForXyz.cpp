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


//! lpaForXyz program: Compute/Report LPA values given XYZ arguments.
int
main
	( int argc
	, char ** argv
	)
{
	if (! (3 < argc))
	{
		std::cerr
			<< "Usage: <progname> xInMeters yInMeters zInMeters"
			<< "\n[xyz]InMeters: Cartesian coordinates in units of meters"
			<< '\n';
	}
	else
	{
		// Example: Important part

		// note: atof values are 0 if decoding error, but not the point here
		peri::XYZ const locXYZ
			{ std::atof(argv[1])
			, std::atof(argv[2])
			, std::atof(argv[3])
			};
		peri::LPA const locLPA{ peri::lpaForXyz(locXYZ) };

		// Echo input values and report converted ones

		// For useful formatting functions, ref. .../tests/periLocal.h
		std::cout << std::fixed ;
		using std::setw;
		std::cout << "input: locXYZ: "
			<< std::setprecision(3)
			<< " " << setw(12) << locXYZ[0]
			<< " " << setw(12) << locXYZ[1]
			<< " " << setw(12) << locXYZ[2]
			<< " " << "[m,m,m]"
			<< std::endl;
		std::cout << "equiv: locLPA: "
			<< std::setprecision(9)
			<< " " << setw(12) << locLPA[0]
			<< " " << setw(12) << locLPA[1]
			<< std::setprecision(3)
			<< " " << setw(12) << locLPA[2]
			<< " " << "[rad,rad,m]"
			<< std::endl;
		double const degPerRad{ 45. / std::atan(1.) };
		std::cout << "equiv: locLPA: "
			<< std::setprecision(7)
			<< " " << setw(12) << degPerRad * locLPA[0]
			<< " " << setw(12) << degPerRad * locLPA[1]
			<< std::setprecision(3)
			<< " " << setw(12) << locLPA[2]
			<< " " << "[deg,deg,m]"
			<< std::endl;
	}
	return 0;
}

