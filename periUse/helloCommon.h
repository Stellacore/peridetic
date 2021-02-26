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


#include <array>  // standard data type used by Peridetic
#include <iostream> // for example data formatting


namespace hello
{
	//! Data space allocation for consuming application
	struct WorkSpace
	{
		std::array<double, 3u> const theDataSrc{ 1.0, 1.1, 1.2 };
		std::array<double, 3u> theDataTmp{};
		std::array<double, 3u> theDataOut{};

	}; // WorkSpace

	//! Put formatted data elements to stream
	std::ostream &
	operator<<
		( std::ostream & ostrm
		, std::array<double, 3u> const & data
		)
	{
		ostrm
			<< std::fixed
			<< " " << data[0]
			<< " " << data[1]
			<< " " << data[2]
			;
		return ostrm;
	}

	//! Put formatted workspace data
	std::ostream &
	operator<<
		( std::ostream & ostrm
		, WorkSpace const & workSpace
		)
	{
		ostrm
			<< std::fixed
			<< " " << "dataSrc: " << workSpace.theDataSrc
			<< '\n'
			<< " " << "dataTmp: " << workSpace.theDataTmp
			<< '\n'
			<< " " << "dataOut: " << workSpace.theDataOut
			;
		return ostrm;
	}
	
}


