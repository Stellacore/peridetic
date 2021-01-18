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
			( std::size_t const & numLon
			, std::size_t const & numPar
			, std::size_t const & numAlt
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
			( DataSet const & dataSet
			)
			: theDataSet{ dataSet }
			, theWorkSpace(theDataSet.size())
		{
		}

		// should optimize away, e.g. for timing only data loop/move overhead
		//! Copy component values unchanged
		inline
		static
		Array
		funcCpy
			( Array const & array
			)
		{
			return array;
		}

		//! abs() of components - simple, but non-trivial to estimate overhead
		inline
		static
		Array
		funcMul
			( Array const & array
			)
		{
			constexpr double co0{ .1 };
			constexpr double co1{ .3 };
			constexpr double co2{ .7 };
			return
				{ co0 * array[0]
				, co1 * array[1]
				, co2 * array[2]
				};
		}

		//! Sqrt(abs()) of each component - for complexity reference
		inline
		static
		Array
		funcSqt
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
		funcXyz
			( Array const & lpa
			)
		{
			return peri::xyzForLpa(lpa, sEarth);
		}

		//! Perform inverse computation (lpaForXyz())
		inline
		static
		Array
		funcLpa
			( Array const & xyz
			)
		{
			return peri::lpaForXyz(xyz, sEarth);
		}

		//! Perform reference computations
		template <typename Func>
		inline
		void
		run
			( std::vector<Array> const & locData
			, Func const & func
			)
		{
			std::transform
				( locData.cbegin()
				, locData.cend()
				, theWorkSpace.insert_iterator()
				, func
				);
		}

		//! Run simple copy operation ('compute' should be optimized away)
		inline
		void
		runCpy
			()
		{
			run(theDataSet.theXyzs, funcCpy);
		}

		//! Perform simple multiplication on every element
		inline
		void
		runMul
			()
		{
			run(theDataSet.theXyzs, funcMul);
		}

		//! Perform square root (of abs()) on every element
		inline
		void
		runSqt
			()
		{
			run(theDataSet.theXyzs, funcSqt);
		}

		//! Perform (easy) forward computations
		inline
		void
		runXyz
			()
		{
			run(theDataSet.theLpas, funcXyz);
		}

		//! Perform (complex) inverse computations
		inline
		void
		runLpa
			()
		{
			run(theDataSet.theXyzs, funcLpa);
		}

	}; // Transformer

	//! basic support for simple 'wall-clock' style timing
	namespace timer
	{
		using namespace std::chrono;

		static high_resolution_clock sClock{};

		//! Utilize std::chrono to get simple timing information
		struct ChronoAboration
		{
			time_point<high_resolution_clock, high_resolution_clock::duration>
				theT0{};
			time_point<high_resolution_clock, high_resolution_clock::duration>
				theT1{};

			void
			start
				()
			{
				theT0 = sClock.now();
			}

			void
			stop
				()
			{
				theT1 = sClock.now();
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

	
} // [eval]


namespace report
{

	//! Return time (in sec) to run func()
	template <typename Func>
	inline
	double
	runTimeFor
		( Func const & func
		)
	{
		eval::timer::ChronoAboration timer{};
		timer.start();
		func();
		timer.stop();
		return timer.elapsed();
	}

	//! Consistently formatted time representation
	std::string
	timeString
		( double const & timeValue
		, std::string const & description = {}
		, std::size_t const & numTailDigits = 9u
		)
	{
		std::ostringstream oss;
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

	//! Pairing of time (sec) and test name
	using TimeName = std::pair<double, std::string>;

	//! relative run times (full values normalized by number of samples)
	std::vector<TimeName>
	perEachRunTimes
		( std::size_t const & numSamps
		, std::vector<TimeName> const & allTimeNames
		)
	{
		std::vector<TimeName> perTimeNames;
		perTimeNames.reserve(allTimeNames.size());
		double const perSamp{ 1. / static_cast<double>(numSamps) };
		for (TimeName const & allTimeName : allTimeNames)
		{
			double const perTime{ perSamp * allTimeName.first };
			std::string const & name = allTimeName.second;
			perTimeNames.emplace_back(std::make_pair(perTime, name));
		}
		return perTimeNames;
	}

	//! Report perTest timing info
	std::string
	absTimingInfo
		( std::vector<report::TimeName> const & allTimeNames
		, std::size_t const & numSamps
		)
	{
		std::ostringstream rpt;
		// report absolute timing (for full tests and average per operation)
		rpt << std::endl;
		rpt << "# Absolute times per test" << '\n';
		rpt << "# -- time values are 'wall-clock' elapsed [in sec]" << '\n';
		rpt << "# -- absolute total and 'per-each' times" << '\n';
		rpt << std::endl;
		std::size_t const numTests{ allTimeNames.size() };
		double const perSamp{ 1. / static_cast<double>(numSamps) };
		for (std::size_t nn{0u} ; nn < numTests ; ++nn)
		{
			double const & allTime = allTimeNames[nn].first;
			std::string const & name = allTimeNames[nn].second;
			double const perTime{ perSamp * allTime };
			rpt
				<< report::timeString(allTime)
				<< " "
				<< report::timeString(perTime, name)
				<< '\n';
		}
		return rpt.str();
	}

	//! Report ratios of test times with respect to each other
	std::string
	relTimeInfo
		( std::vector<report::TimeName> const & allTimeNames
		)
	{
		std::ostringstream rpt;
		rpt << std::endl;
		rpt << "# Relative test times" << '\n';
		rpt << "# -- times tests with respect to each other [ratio]" << '\n';
		rpt << "# -- column order matches row order" << '\n';
		rpt << std::endl;
		std::size_t const & numTests = allTimeNames.size();
		for (std::size_t nCurr{0u} ; nCurr < numTests ; ++nCurr)
		{
			double const & currTime = allTimeNames[nCurr].first;
			std::string const & currName = allTimeNames[nCurr].second;
			for (std::size_t nBase{0u} ; nBase < numTests ; ++nBase)
			{
				double const & baseTime = allTimeNames[nBase].first;
				double const relTime{ currTime / baseTime };
				constexpr std::size_t nDig{ 2u };
				rpt << " "
					<< std::fixed << std::setprecision(nDig)
					<< std::setw(3u+1u+nDig)
					<< relTime;
			}
			rpt << "  : " << currName << std::endl;
		}
		return rpt.str();
	}



} // [report]

int
main
	()
{
	// 32 leads to about same number points as sq-deg in a sphere (41253)
	// constexpr std::size_t num1D{ 32u };
	constexpr std::size_t num1D{ 256u };
	constexpr std::size_t numLon{ num1D };
	constexpr std::size_t numPar{ num1D };
	constexpr std::size_t numAlt{ num1D };

	std::cout << "--- setup: " << std::endl;

	// allocate data with pre-set values in both domains
	eval::DataSet const data(numLon, numPar, numAlt);
	std::size_t const numSamps{ data.size() };
	// allocate fixed workspace for available transformations
	eval::Transformer xformer(data);

	using Func_t = std::function<void(void)>;
	Func_t const funcCpy{ std::bind(&eval::Transformer::runCpy, xformer) };
	Func_t const funcMul{ std::bind(&eval::Transformer::runMul, xformer) };
	Func_t const funcSqt{ std::bind(&eval::Transformer::runSqt, xformer) };
	Func_t const funcXyz{ std::bind(&eval::Transformer::runXyz, xformer) };
	Func_t const funcLpa{ std::bind(&eval::Transformer::runLpa, xformer) };

	std::string const nameCpy{ "Reference evaluation - copy: " };
	std::string const nameMul{ "Reference evaluation - multiply: " };
	std::string const nameSqt{ "Reference evaluation - sqrt(abs()): " };
	std::string const nameXyz{ "Cartesian from Geodetic - xyzForLpa(): " };
	std::string const nameLpa{ "Geodetic from Cartesian - lpaForXyz(): " };

	// run each computation test and note (wall) time it takes
	double const timeCpy{ report::runTimeFor(funcCpy) };
	double const timeMul{ report::runTimeFor(funcMul) };
	double const timeSqt{ report::runTimeFor(funcSqt) };
	double const timeXyz{ report::runTimeFor(funcXyz) };
	double const timeLpa{ report::runTimeFor(funcLpa) };

	// gather results for use in reporting
	std::vector<report::TimeName> const allTimeNames
		{ std::make_pair(timeCpy, nameCpy)
		, std::make_pair(timeMul, nameMul)
		, std::make_pair(timeSqt, nameSqt)
		, std::make_pair(timeXyz, nameXyz)
		, std::make_pair(timeLpa, nameLpa)
		};

	std::cout << "--- reporting: " << std::endl;

	// report test stats
	std::ostringstream rpt;

	// report test stats
	rpt << std::endl;
	rpt << "# Number samples tested: " << numSamps << std::endl;

	// report absolute timing (for full tests and average per operation)
	rpt << report::absTimingInfo(allTimeNames, numSamps) << std::endl;

	// generate table of times relative to each other
	rpt << report::relTimeInfo(allTimeNames);

	// display results
	std::cout << rpt.str() << std::endl;

	return 0;
}

