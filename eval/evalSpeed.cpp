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


//#include "peridetic.h"

// #include "periLocal.h"
#include "periSim.h"

#include <algorithm>
#include <chrono>
#include <iostream>
#include <vector>


namespace eval
{
	using Array = std::array<double, 3u>;

	static peri::EarthModel const sEarth{ peri::model::WGS84 };

	//! Convert collection of LPAs into XYZs
	inline
	std::vector<peri::XYZ>
	xyzsFor
		( std::vector<peri::LPA> const & lpas
		)
	{
		std::vector<peri::XYZ> xyzs;
		xyzs.reserve(lpas.size());
		for (peri::LPA const & lpa : lpas)
		{
			xyzs.emplace_back(sEarth.xyzForLpa(lpa));
		}
		return xyzs;
	}

	//! A complete data set of pre-computed values
	struct DataSet
	{
		std::vector<peri::LPA> const theLpas;
		std::vector<peri::XYZ> const theXyzs;

		inline
		explicit
		DataSet
			( std::size_t const & numLon = 128u
			, std::size_t const & numPar = 128u
			, std::size_t const & numAlt = 128u
			)
			: theLpas{ peri::sim::bulkSamplesLpa(numLon, numPar, numAlt) }
			, theXyzs{ xyzsFor(theLpas) }
		{
		}

		inline
		std::size_t
		size
			() const
		{
			return theLpas.size();
		}

	}; // DataSet

	//! A pre-allocated (upon construction) workspace
	struct WorkSpace
	{
		std::vector<Array> theSpace{};

		inline
		explicit
		WorkSpace
			( std::size_t const & size
			)
			: theSpace{}
		{ 
			// ensure workspace is preallocated
			theSpace.reserve(size);
		}

		std::insert_iterator<std::vector<Array> >
		insert_iterator
			()
		{
			return std::insert_iterator(theSpace, theSpace.begin());
		}


	}; // WorkSpace

	//! A group of transformations useful for speed evaluation
	struct Transformer
	{
		DataSet theDataSet;
		WorkSpace theWorkSpace;

		inline
		explicit
		Transformer
			()
			: theDataSet{}
			, theWorkSpace(theDataSet.size())
		{
		}

		//! abs() of components - simple, but non-trivial to estimate overhead
		inline
		static
		Array
		absOf // arbitrary simple computation as reference for timing
			( Array const & array
			)
		{
			return
				{ std::sqrt(std::abs(array[0]))
				, std::sqrt(std::abs(array[1]))
				, std::sqrt(std::abs(array[2]))
				};
		}


		//! Sqrt(abs()) of each component - for complexity reference
		inline
		static
		Array
		sqrtAbsOf // arbitrary simple computation as reference for timing
			( Array const & array
			)
		{
			return
				{ std::sqrt(std::abs(array[0]))
				, std::sqrt(std::abs(array[1]))
				, std::sqrt(std::abs(array[2]))
				};
		}

		//! Perform forward computation (xyzForLpa())
		inline
		static
		Array
		toXyz
			( Array const & lpa
			)
		{
			return peri::xyzForLpa(lpa, sEarth);
		}

		//! Perform inverse computation (lpaForXyz())
		inline
		static
		Array
		toLpa
			( Array const & xyz
			)
		{
			return peri::lpaForXyz(xyz, sEarth);
		}

		//! Perform reference computations
		template <typename Func>
		inline
		void
		compute
			( Func const & func
			)
		{
			std::transform
				( theDataSet.theXyzs.cbegin()
				, theDataSet.theXyzs.cend()
				, theWorkSpace.insert_iterator()
				, func
				);
		}

		//! Perform reference computations
		inline
		void
		computeAbs
			()
		{
			compute(absOf);
		}

		//! Perform reference computations
		inline
		void
		computeSqrt
			()
		{
			compute(sqrtAbsOf);
		}

		//! Perform (easy) forward computations
		inline
		void
		computeToXyz
			()
		{
			compute(toXyz);
		}

		//! Perform (complex) inverse computations
		inline
		void
		computeToLpa
			()
		{
			compute(toLpa);
		}


	}; // Transformer

	namespace timer
	{
		using namespace std::chrono;

		//! Utilize std::chrono to get simple timing information
		struct ChronoAboration
		{
			high_resolution_clock theClock{};
			time_point<high_resolution_clock, high_resolution_clock::duration>
				theT0{};
			time_point<high_resolution_clock, high_resolution_clock::duration>
				theT1{};

			void
			start
				()
			{
				theT0 = theClock.now();
			}

			void
			stop
				()
			{
				theT1 = theClock.now();
			}

			//! Time between t0 and t1 in seconds
			double
			elapsed
				() const
			{
				duration<double> const delta{ theT1 - theT0 };
				return delta.count();
			}

		}; // ChronoAboration

	} // [timer]

	
	//! Return time (in sec) to run func()
	template <typename Func>
	inline
	double
	timeToRun
		( Func const & func
		)
	{
		timer::ChronoAboration timer{};
		timer.start();
		func();
		timer.stop();
		return timer.elapsed();
	}

	std::string
	timeString
		( double const & timeValue
		, std::string const & description = {}
		)
	{
		std::ostringstream oss;
		constexpr std::size_t numTailDigits{ 9u };
		oss
			<< std::fixed
			<< std::setw(5u+1u+numTailDigits)
			<< std::setprecision(numTailDigits)
			<< timeValue
			;
		if (! description.empty())
		{
			oss << "  : ";
			oss
				<< std::setw(40u)
				<< std::left
				<< description
				;
		}
		return oss.str();
	}

} // [annon]


int
main
	()
{
	constexpr std::size_t num1D{ 256u };
	constexpr std::size_t numLon{ num1D };
	constexpr std::size_t numPar{ num1D };
	constexpr std::size_t numAlt{ num1D };

	std::cout << "--- setup: " << std::endl;

	eval::DataSet const data(numLon, numPar, numAlt);
	eval::Transformer xformer{};

	// configure tests
	std::function<void(void)> const funcAbs
		{ std::bind(&eval::Transformer::computeAbs, xformer) };
	std::function<void(void)> const funcSqt
		{ std::bind(&eval::Transformer::computeSqrt, xformer) };
	std::function<void(void)> const funcXyz
		{ std::bind(&eval::Transformer::computeToXyz, xformer) };
	std::function<void(void)> const funcLpa
		{ std::bind(&eval::Transformer::computeToLpa, xformer) };

	std::string const nameAbs{ "Reference evaluations - abs(): " };
	std::string const nameSqt{ "Reference evaluations - sqrt(abs()): " };
	std::string const nameXyz{ "Cartesian from Geodetic - xyzForLpa(): " };
	std::string const nameLpa{ "Geodetic from Cartesian - lpaForXyz(): " };

	std::cout << "--- evaluating: " << std::endl;

	double const timeAbs{ eval::timeToRun(funcAbs) };
	double const timeSqt{ eval::timeToRun(funcSqt) };
	double const timeXyz{ eval::timeToRun(funcXyz) };
	double const timeLpa{ eval::timeToRun(funcLpa) };

	std::cout << "--- reporting: " << std::endl;

	double const dCount{ static_cast<double>(data.size()) };

	double const perEachAbs{ timeAbs / dCount };
	double const perEachSqt{ timeSqt / dCount };
	double const perEachXyz{ timeXyz / dCount };
	double const perEachLpa{ timeLpa / dCount };

	double const & perEachRef = perEachAbs;

	double const relativeAbs{ perEachAbs / perEachRef };
	double const relativeSqt{ perEachSqt / perEachRef };
	double const relativeXyz{ perEachXyz / perEachRef };
	double const relativeLpa{ perEachLpa / perEachRef };


	std::cout << std::endl;
	std::cout << "# Number samples tested: " << data.size() << std::endl;

	std::cout << std::endl;
	std::cout << "# absolute times ('wall-clock' time in sec)" << std::endl;
	std::cout << eval::timeString(timeAbs, nameAbs) << std::endl;
	std::cout << eval::timeString(timeSqt, nameSqt) << std::endl;
	std::cout << eval::timeString(timeXyz, nameXyz) << std::endl;
	std::cout << eval::timeString(timeLpa, nameLpa) << std::endl;

	std::cout << std::endl;
	std::cout << "# per each times (in sec)" << std::endl;
	std::cout << eval::timeString(perEachAbs, nameAbs) << std::endl;
	std::cout << eval::timeString(perEachSqt, nameSqt) << std::endl;
	std::cout << eval::timeString(perEachXyz, nameXyz) << std::endl;
	std::cout << eval::timeString(perEachLpa, nameLpa) << std::endl;

	std::cout << std::endl;
	std::cout << "# relative times (ratio to 'abs')" << std::endl;
	std::cout << eval::timeString(relativeAbs, nameAbs) << std::endl;
	std::cout << eval::timeString(relativeSqt, nameSqt) << std::endl;
	std::cout << eval::timeString(relativeXyz, nameXyz) << std::endl;
	std::cout << eval::timeString(relativeLpa, nameLpa) << std::endl;

	std::cout << std::endl;

	return 0;
}

