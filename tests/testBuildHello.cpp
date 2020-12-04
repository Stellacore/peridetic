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


#include "peridetic.h" // check project include

#include <iostream>
#include <sstream>
#include <string>

//! Example unit test structure
int
main
	()
{
	// use error count as ctest return status indicator
	int errCount{ 0 };
	// conduct test
	std::ostringstream oss;
	std::string const expMsg("hello");
	oss << expMsg;
	std::string const gotMsg{ oss.str() };
	if (! (gotMsg == expMsg))
	{
		// use "ctest --verbose" to see these type messages
		std::cerr << "Failure of simple 'hello' test" << std::endl;
		std::cerr << "expMsg: '" << expMsg << "'" << std::endl;
		std::cerr << "gotMsg: '" << gotMsg << "'" << std::endl;
		// increment error count for all errors
		++errCount;
	}
	return errCount;
}

