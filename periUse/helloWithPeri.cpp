//
//
// MIT License
//
// Copyright (c) 2021 Stellacore Corporation.
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


#include <peridetic.h>  // include for Peridetic usage

// common data store and I/O code for both hello{Sans,With}Peri.cpp
#include "helloCommon.h"


int
main
	()
{
	// allocate same data space for both "with" and "sans" cases
	hello::WorkSpace workSpace;

	// perform data transformation operations - *WITH* Peridetic use
	workSpace.theDataTmp = peri::xyzForLpa(workSpace.theDataSrc);
	workSpace.theDataOut = peri::lpaForXyz(workSpace.theDataTmp);

	// output data values to confirm that
	// above code operations are actually present in executable
	using hello::operator<<;
	std::cout << workSpace << std::endl;
	return 0;
}

