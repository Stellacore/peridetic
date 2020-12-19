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

namespace peri
{
	//! Evaluate scalar function: psi = sum(q_k^2/mu_k)-1.
	inline
	double
	psiValueAt
		( XYZ const & xVec
		, peri::Shape const & shape
		)
	{
		return
			{ sq(xVec[0]) / shape.theMuSqs[0]
			+ sq(xVec[1]) / shape.theMuSqs[1]
			+ sq(xVec[2]) / shape.theMuSqs[2]
			- 1.
			};
	}

	/*! \brief Distance from origin to surface of ellipsoid toward xyzAny.
	 *
	 * NOTE: xyzAny can be any \b non-zero arbitrary location.
	 */
	double
	ellipsoidRadiusToward
		( XYZ const & xyzAny
		, peri::Shape const & shape
		)
	{
		// ellipsoid shape parameters
		double const & alpha = shape.theRadA;
		double const & beta = shape.theRadB;
		double const boa{ beta / alpha };
		double const eSq{ 1. - boa*boa };

		// shorthand aliases
		double const & x1 = xyzAny[0];
		double const & x2 = xyzAny[1];
		double const & x3 = xyzAny[2];

		double const x3Sq{ x3*x3 }; // square of polar constituent
		double const hSq{ x1*x1 + x2*x2 }; // square of equatorial constituent
		double const xSq{ hSq + x3Sq }; // squared magnitude of xyzAny
		double const cosPar{ hSq / xSq }; // cosine of elevation from equator
		double const den{ 1. - eSq * cosPar };
		double const frac{ std::sqrt(1. / den) }; // fraction of polar radius
		double const rho{ beta * frac }; // radial magnitude
		return rho;
	}

	//! \brief Centric vector to surface of ellipsoid in direction of xyzAny.
	inline
	XYZ
	ellipsoidVectorToward
		( XYZ const & xyzAny
		, peri::Shape const & shape
		)
	{
		return ellipsoidRadiusToward(xyzAny, shape) * unit(xyzAny);
	}

	/*! \brief Compute (wastefully) various mathematical quantities at an LPA.
	 *
	 */
	struct Info
	{
		Shape const theShape;
		Ellipsoid const theEllip;
		EarthModel const theEarth;
		std::array<double, 3u> const theMus;
		std::array<double, 3u> const theMuSqs;
		LPA const theLpa; // Location under consideration
		XYZ const theVecX; // vector location for theLpa
		XYZ const theVecP; // perpendicular footer point on ellipsoid
		XYZ const theVecR; // radially aligned point on ellipsoid
		XYZ const theVecG; // Gradient vector (perpendicular to ellipsoid)
		double const theEta;
		double const thePsi; // scalar function (misclosure) value

		//! Compute values for various quantities and relationships at lpa.
		explicit
		Info
			( LPA const & lpa
			, Shape const & shape = shape::sWGS84 // .normalizedShape()
			)
			: theShape{ shape }
			, theEllip(theShape)
			, theEarth(theShape)
			, theMus
				{ theShape.theRadA
				, theShape.theRadA
				, theShape.theRadB
				}
			, theMuSqs
				{ sq(theMus[0])
				, sq(theMus[1])
				, sq(theMus[2])
				}
			, theLpa{ lpa }
			, theVecX{ xyzForLpa(theLpa, theEarth) }
			, theVecP{ xyzForLpa(LPA{ theLpa[0], theLpa[1], 0. }, theEarth) }
			, theVecR{ ellipsoidVectorToward(theVecX, theShape) }
			, theVecG{ theShape.gradientAt(theVecP) }
			//
			, theEta{ theLpa[2] }
			, thePsi{ psiValueAt(theVecX, theShape) }
		{
		}

		inline
		double
		magG
			() const
		{
			double sum{ 0. };
			for (std::size_t kk{0u} ; kk < 3u ; ++kk)
			{
				double const mu4{ sq(theMuSqs[kk]) };
				double const pSq{ sq(theVecP[kk]) };
				sum += pSq/mu4;
			}
			return std::sqrt(4. * sum);
		}

		inline
		XYZ
		upComputed
			() const
		{
			XYZ const grad{ theShape.gradientAt(theVecP) };
			return { (1./magG()) * grad };
		}

		inline
		double
		sigmaValue
			() const
		{
			return { 2.*theEta / magG() };
		}

		//! Radial altitude
		inline
		double
		rAlt
			() const
		{
			return magnitude(theVecX - theVecR);
		}

		//! Perpendicular (proper) altitude
		inline
		double
		pAlt
			() const
		{
			return magnitude(theVecX - theVecP);
		}

		//! Difference in altitude representations (rAlt-pAlt)
		inline
		double
		deltaAlt
			() const
		{
			return (rAlt() - pAlt());
		}

		inline
		std::array<double, 3u>
		sigmaMuSqs
			() const
		{
			double const sigma{ sigmaValue() };
			return
				{ sigma / theMuSqs[0]
				, sigma / theMuSqs[1]
				, sigma / theMuSqs[2]
				};
		}

		//! Descriptive information about this instance
		std::string
		infoString
			( std::string const & title = {}
			) const
		{
			std::ostringstream oss;
			if (! title.empty())
			{
				oss << title << std::endl;
			}
			oss << xyz::infoString(theMus, "theMus");
			oss << std::endl;
			oss << lpa::infoString(theLpa, "theLpa");
			oss << std::endl;
			oss << xyz::infoString(theVecX, "theVecX")
				<< " " << string::fixedLinear(magnitude(theVecX), "magX");
			oss << std::endl;
			oss << xyz::infoString(theVecR, "theVecR")
				<< " " << string::fixedLinear(magnitude(theVecR), "magR");
			oss << std::endl;
			oss << xyz::infoString(theVecP, "theVecP")
				<< " " << string::fixedLinear(magnitude(theVecP), "magP");
			oss << std::endl;
			oss << string::allDigits(theVecG, "theVecG")
				<< " " << string::allDigits(magnitude(theVecG), "magG");
			oss << std::endl;
			oss << xyz::infoString(upComputed(), "upComputed()")
				<< " " << string::allDigits(magnitude(upComputed()), "upMag");
			oss << std::endl;
			oss << xyz::infoString(unit(theVecG), "gHat");
			oss << std::endl;
			oss << xyz::infoString(upComputed(), "up()");
			oss << std::endl;
			oss << string::fixedLinear(rAlt(), "rAlt()")
				<< "  " << string::allDigits(rAlt(), "rAlt()");
			oss << std::endl;
			oss << string::fixedLinear(pAlt(), "pAlt()")
				<< "  " << string::allDigits(pAlt(), "pAlt()");
			oss << std::endl;
			oss << string::fixedLinear(deltaAlt(), "deltaAlt()")
				<< "  " << string::allDigits(deltaAlt(), "deltaAlt()");
			oss << std::endl;
			oss << string::fixedLinear(theEta, "theEta");
			oss << std::endl;
			oss << string::allDigits(theMuSqs, "theMuSqs");
			oss << std::endl;
			oss << string::allDigits(sigmaValue(), "sigmaValue()");
			oss << std::endl;
			oss << string::allDigits(sigmaMuSqs(), "sigmaMuSqs()");

			return oss.str();
		}

	}; // Info

} // [annon]

namespace
{
	//! Check formula for restoring physical altitude
	int
	checkEta
		( peri::Info const & info
		)
	{
		int errCount{ 0 };
		double const & expEta = info.theEta;

		double const sigma{ info.sigmaValue() };

		double sum{ 0. };
		for (std::size_t kk{0u} ; kk < 3u ; ++kk)
		{
			sum += peri::sq(info.theVecX[kk] / (info.theMuSqs[kk] + sigma));
		}
		double const gotEta{ sigma * std::sqrt(sum) };

		constexpr double tol{ 1.e-14 };
		if (! peri::sameEnough(gotEta, expEta, tol))
		{
			std::cerr << "FAILURE checkEta() error" << std::endl;
			std::cerr << peri::string::allDigits(expEta, "expEta") << '\n';
			std::cerr << peri::string::allDigits(gotEta, "gotEta") << '\n';
			std::cerr << peri::string::allDigits(tol, "tol") << '\n';
			++errCount;
		}
		else
		{
			std::cout << "Success checkEta()" << std::endl;
			std::cerr << peri::string::allDigits(expEta, "expEta") << '\n';
			std::cerr << peri::string::allDigits(gotEta, "gotEta") << '\n';
		}
		return errCount;
	}

	//! Check formula for restoring point of interest
	int
	checkX
		( peri::Info const & info
		)
	{
		int errCount{ 0 };

		peri::XYZ const & expX = info.theVecX;

		using peri::operator+;
		using peri::operator*;

		/* - okay
		double const & eta = info.theEta;
		peri::XYZ const up{ info.upComputed() };
		peri::XYZ const gotX{ info.theVecP + eta * up };
		*/

		/* - okay
		*/
		double const sigma = info.sigmaValue();
		peri::XYZ const grad{ info.theShape.gradientAt(info.theVecP) };
		peri::XYZ const gotX{ info.theVecP + .5 * sigma * grad };

		constexpr double tol{ 1.e-14 };
		if (! peri::xyz::sameEnough(gotX, expX, tol))
		{
			std::cerr << "FAILURE checkX() error" << std::endl;
			std::cerr << peri::string::allDigits(expX, "expX") << '\n';
			std::cerr << peri::string::allDigits(gotX, "gotX") << '\n';
			std::cerr << peri::string::allDigits(tol, "tol") << '\n';
			++errCount;
		}
		else
		{
			std::cout << "Success checkX()" << std::endl;
			std::cerr << peri::string::allDigits(expX, "expX") << '\n';
			std::cerr << peri::string::allDigits(gotX, "gotX") << '\n';
		}
		return errCount;
	}

	/*
	//! Check TODO
	int
	test1
		()
	{
		int errCount{ 0 };
		peri::LPA const lpa{ .5*peri::pi(), .25*peri::pi(), .25 };
		peri::Info const info(lpa);

		// check if footer point matches 0-altitude point

		std::cout << info.infoString("info") << std::endl;
		std::cout << "====" << '\n';

		std::cout << "checkEta:\n";
		errCount += checkEta(info);

		std::cout << "checkX:\n";
		errCount += checkX(info);

		return errCount;
	}
	*/
}


//! Check various ellipsoid traits and relationships
int
main
	()
{
	int errCount{ 0 };
	peri::LPA const lpa{ .5*peri::pi(), .25*peri::pi(), .25 };
	peri::Info const info(lpa);
	errCount += checkEta(info);
	errCount += checkX(info);
	return errCount;
}

