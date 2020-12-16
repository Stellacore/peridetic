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
			, Shape const & shape = shape::sWGS84
			)
			: theShape{ shape }
			, theEllip(theShape)
			, theEarth(theShape)
			, theMus
				{ theShape.theRadA
				, theShape.theRadA
				, theShape.theRadB
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
				double const mu4{ sq(sq(theMus[kk])) };
				double const pSq{ sq(theVecP[kk]) };
				sum += pSq/mu4;
			}
			return std::sqrt(sum);
		}

		inline
		XYZ
		upComputed
			() const
		{
			double const scale{ .5 / magG() };
			XYZ const grad{ theShape.gradientAt(theVecP) };
			return { scale * grad };
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
				{ sigma / sq(theMus[0])
				, sigma / sq(theMus[1])
				, sigma / sq(theMus[2])
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
			oss << xyz::infoString(theVecX, "theVecX");
			oss << std::endl;
			oss << string::fixedLinear(magnitude(theVecX), "magX");
			oss << std::endl;
			oss << xyz::infoString(theVecR, "theVecR");
			oss << std::endl;
			oss << xyz::infoString(theVecP, "theVecP");
			oss << std::endl;
			oss << string::allDigits(theVecG, "theVecG");
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
			oss << string::allDigits(sigmaValue(), "sigmaValue()");
			oss << std::endl;
			oss << string::allDigits(sigmaMuSqs(), "sigmaMuSqs()");

			return oss.str();
		}

	}; // Info

} // [annon]

namespace
{
	//! Check TODO
	int
	test1
		()
	{
		int errCount{ 0 };
		peri::LPA const lpa{ .5*peri::pi(), .25*peri::pi(), 1000. };
		peri::Info const info(lpa);

		// check if footer point matches 0-altitude point

		std::cout << info.infoString("info") << std::endl;
++errCount;

		/*
		// normalize values
		peri::XYZ const xyzNear
			{ info.theEarth.nearEllipsoidPointFor(info.theVecX) };
		std::cout << peri::xyz::infoString(xyzNear, "xyzNear") << std::endl;
		*/

		return errCount;
	}
}


//! Check various ellipsoid traits and relationships
int
main
	()
{
	int errCount{ 0 };
	errCount += test1();
	return errCount;
}

