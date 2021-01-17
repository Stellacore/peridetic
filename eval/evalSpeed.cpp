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
		computeReference
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
		, std::string const & title = {}
		)
	{
		std::ostringstream oss;
		if (! title.empty())
		{
			oss << std::setw(16u) << title << " ";
		}
		constexpr std::size_t numTailDigits{ 9u };
		oss
			<< std::fixed
			<< std::setw(5u+1u+numTailDigits)
			<< std::setprecision(numTailDigits)
			<< timeValue
			;
		return oss.str();
	}

} // [annon]


int
main
	()
{
	std::cout << "hi" << std::endl;

	eval::DataSet const data{};
	eval::Transformer xformer{};

	std::function<void(void)> const funcRef
		{ std::bind(&eval::Transformer::computeReference, xformer) };
	std::function<void(void)> const funcXyz
		{ std::bind(&eval::Transformer::computeToXyz, xformer) };
	std::function<void(void)> const funcLpa
		{ std::bind(&eval::Transformer::computeToLpa, xformer) };

	double const dCount{ static_cast<double>(data.size()) };
	double const perEachRef{ eval::timeToRun(funcRef) / dCount };
	double const perEachXyz{ eval::timeToRun(funcXyz) / dCount };
	double const perEachLpa{ eval::timeToRun(funcLpa) / dCount };

	double const relativeRef{ perEachRef / perEachRef };
	double const relativeXyz{ perEachXyz / perEachRef };
	double const relativeLpa{ perEachLpa / perEachRef };

	std::cout << eval::timeString(perEachRef, "perEachRef") << std::endl;
	std::cout << eval::timeString(perEachXyz, "perEachXyz") << std::endl;
	std::cout << eval::timeString(perEachLpa, "perEachLpa") << std::endl;

	std::cout << eval::timeString(relativeRef, "relativeRef") << std::endl;
	std::cout << eval::timeString(relativeXyz, "relativeXyz") << std::endl;
	std::cout << eval::timeString(relativeLpa, "relativeLpa") << std::endl;

	return 0;
}

