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
			, Shape const & shape
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
			, theVecR{ ellip::vectorToward(theVecX, theShape) }
			, theVecG{ theShape.gradientAt(theVecP) }
			//
			, theEta{ theLpa[2] }
			, thePsi{ psiValueAt(theVecX, theShape) }
		{
		}

		inline
		double
		nuAve
			() const
		{
			return std::sqrt(theShape.theRadA * theShape.theRadB);
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
			double sum{ 0. };
			for (std::size_t kk{0u} ; kk < 3u ; ++kk)
			{
				double const mu4{ sq(theMuSqs[kk]) };
				double const pSq{ sq(theVecP[kk]) };
				sum += pSq/mu4;
			}
			return { theEta / std::sqrt(sum) };
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
			oss << string::allDigits(sigmaValue(), "sigma");
			oss << std::endl;
			oss << string::allDigits(nuAve()*theEta, "nu*eta");
			oss << std::endl;
			oss << string::allDigits(sigmaMuSqs(), "sigmaMuSqs()");

			return oss.str();
		}

	}; // Info

} // [annon]

namespace
{
	/*
	//! Check TODO
	int
	checkQ
		( peri::Info const & info
		)
	{
		int errCount{ 0 };

		double const expFunc{ 0. };
		double const gotFunc{ peri::sNan };

		constexpr double tol{ 1.e-14 };
		if (! peri::sameEnough(gotFunc, expFunc, tol))
		{
			std::cerr << "FAILURE checkFunc() error" << std::endl;
			std::cerr << peri::string::allDigits(expFunc, "expFunc") << '\n';
			std::cerr << peri::string::allDigits(gotFunc, "gotFunc") << '\n';
			std::cerr << peri::string::allDigits(tol, "tol") << '\n';
			++errCount;
		}
		else
		{
			std::cout << "Success checkFunc()" << std::endl;
			std::cerr << peri::string::fixedLinear(expFunc, "expFunc") << '\n';
			std::cerr << peri::string::fixedLinear(gotFunc, "gotFunc") << '\n';
		}

		return errCount;
	}
	*/

	//! Check that point, p, is on ellipsoid
	int
	checkOnEllip
		( peri::Info const & info
		)
	{
		int errCount{ 0 };

		double const expFunc{ 0. };

		// eqn: DefEllipsoid
		double sum{ 0. };
		for (std::size_t kk{0u} ; kk < 3u ; ++kk)
		{
			sum += peri::sq(info.theVecP[kk]) / info.theMuSqs[kk];
		}
		double const gotFunc{ sum - 1. };

		constexpr double tol{ 1.e-14 };
		if (! peri::sameEnough(gotFunc, expFunc, tol))
		{
			std::cerr << "FAILURE checkFunc() error" << std::endl;
			std::cerr << peri::string::allDigits(expFunc, "expFunc") << '\n';
			std::cerr << peri::string::allDigits(gotFunc, "gotFunc") << '\n';
			std::cerr << peri::string::allDigits(tol, "tol") << '\n';
			++errCount;
		}
		else
		{
			std::cout << "Success checkFunc()" << std::endl;
			std::cerr << peri::string::fixedLinear(expFunc, "expFunc") << '\n';
			std::cerr << peri::string::fixedLinear(gotFunc, "gotFunc") << '\n';
		}

		return errCount;
	}

	//! Check gradient computations
	int
	checkGrad
		( peri::Info const & info
		)
	{
		int errCount{ 0 };

		peri::XYZ const & expGrad = info.theVecG;
		// eqn: DefGradient
		peri::XYZ const gotGrad
			{ 2. * info.theVecP[0] / info.theMuSqs[0]
			, 2. * info.theVecP[1] / info.theMuSqs[1]
			, 2. * info.theVecP[2] / info.theMuSqs[2]
			};

		constexpr double tol{ 1.e-14 };
		if (! peri::xyz::sameEnough(gotGrad, expGrad, tol))
		{
			std::cerr << "FAILURE checkGrad() error" << std::endl;
			std::cerr << peri::string::allDigits(expGrad, "expGrad") << '\n';
			std::cerr << peri::string::allDigits(gotGrad, "gotGrad") << '\n';
			std::cerr << peri::string::allDigits(tol, "tol") << '\n';
			++errCount;
		}
		else
		{
			std::cout << "Success checkGrad()" << std::endl;
			std::cerr << peri::xyz::infoString(expGrad, "expGrad") << '\n';
			std::cerr << peri::xyz::infoString(gotGrad, "gotGrad") << '\n';
		}

		return errCount;
	}

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

		double const tol{ 1.e-14 * expEta };
		if (! peri::sameEnough(gotEta, expEta, tol))
		{
			double const dif{ gotEta - expEta };
			std::cerr << "FAILURE checkEta() error" << std::endl;
			std::cerr << peri::string::allDigits(expEta, "expEta") << '\n';
			std::cerr << peri::string::allDigits(gotEta, "gotEta") << '\n';
			std::cerr << peri::string::allDigits(dif, "dif") << '\n';
			std::cerr << peri::string::allDigits(tol, "tol") << '\n';
			++errCount;
		}
		else
		{
			std::cout << "Success checkEta()" << std::endl;
			std::cerr << peri::string::fixedLinear(expEta, "expEta") << '\n';
			std::cerr << peri::string::fixedLinear(gotEta, "gotEta") << '\n';
		}

		return errCount;
	}

	//! Check formula sigma evaluation
	int
	checkSigma
		( peri::Info const & info
		)
	{
		int errCount{ 0 };

		double const expSig{ info.sigmaValue() };

		double sum{ 0. };
		for (std::size_t jj{0u} ; jj < 3u ; ++jj)
		{
			sum += peri::sq(info.theVecP[jj] / info.theMuSqs[jj]);
		}
		double const gotSig{ info.theEta / std::sqrt(sum) };

		constexpr double tol{ 1.e-14 };
		if (! peri::sameEnough(gotSig, expSig, tol))
		{
			std::cerr << "FAILURE checkSig() error" << std::endl;
			std::cerr << peri::string::allDigits(expSig, "expSig") << '\n';
			std::cerr << peri::string::allDigits(gotSig, "gotSig") << '\n';
			std::cerr << peri::string::allDigits(tol, "tol") << '\n';
			++errCount;
		}
		else
		{
			std::cout << "Success checkSig()" << std::endl;
			std::cerr << peri::string::fixedLinear(expSig, "expSig") << '\n';
			std::cerr << peri::string::fixedLinear(gotSig, "gotSig") << '\n';
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
		double const tol{ 1.e-14 * peri::magnitude(expX) };

		// - okay
		// eqn: x = p + eta*up
		double const & eta = info.theEta;
		peri::XYZ const up{ info.upComputed() };
		peri::XYZ const gotXEta{ info.theVecP + eta * up };

		if (! peri::xyz::sameEnough(gotXEta, expX, tol))
		{
			std::cerr << "FAILURE checkX() error (eta*up)" << std::endl;
			std::cerr << peri::string::allDigits(expX, "expX") << '\n';
			std::cerr << peri::string::allDigits(gotXEta, "gotXEta") << '\n';
			std::cerr << peri::string::allDigits(tol, "tol") << '\n';
			++errCount;
		}
		else
		{
			std::cout << "Success checkX()" << std::endl;
			std::cerr << peri::xyz::infoString(expX, "expX") << '\n';
			std::cerr << peri::xyz::infoString(gotXEta, "gotXEta") << '\n';
		}

		/* - okay
		*/
		// eqn: xInTwoSteps_sg : x = p +.5*sigma*grad : checks sigma,grad
		double const sigma{ info.sigmaValue() };
		peri::XYZ const grad{ info.theShape.gradientAt(info.theVecP) };
		peri::XYZ const gotXSig{ info.theVecP + .5 * sigma * grad };

		if (! peri::xyz::sameEnough(gotXSig, expX, tol))
		{
			std::cerr << "FAILURE checkX() error (sigma*grad)" << std::endl;
			std::cerr << peri::string::allDigits(expX, "expX") << '\n';
			std::cerr << peri::string::allDigits(gotXSig, "gotXSig") << '\n';
			std::cerr << peri::string::allDigits(tol, "tol") << '\n';
			++errCount;
		}
		else
		{
			std::cout << "Success checkX()" << std::endl;
			std::cerr << peri::xyz::infoString(expX, "expX") << '\n';
			std::cerr << peri::xyz::infoString(gotXSig, "gotXSig") << '\n';
		}

		return errCount;
	}

	//! Check point, p, components computation
	int
	checkP
		( peri::Info const & info
		)
	{
		int errCount{ 0 };

		double const sigma{ info.sigmaValue() };

		peri::XYZ const & expP = info.theVecP;
		// eqn: pkFuncOfSigma
		peri::XYZ const gotP
			{ (info.theVecX[0] * info.theMuSqs[0]) / (info.theMuSqs[0] + sigma)
			, (info.theVecX[1] * info.theMuSqs[1]) / (info.theMuSqs[1] + sigma)
			, (info.theVecX[2] * info.theMuSqs[2]) / (info.theMuSqs[2] + sigma)
			};

		constexpr double tol{ 1.e-14 };
		if (! peri::xyz::sameEnough(gotP, expP, tol))
		{
			std::cerr << "FAILURE checkP() error" << std::endl;
			std::cerr << peri::string::allDigits(expP, "expP") << '\n';
			std::cerr << peri::string::allDigits(gotP, "gotP") << '\n';
			std::cerr << peri::string::allDigits(tol, "tol") << '\n';
			++errCount;
		}
		else
		{
			std::cout << "Success checkP()" << std::endl;
			std::cerr << peri::xyz::infoString(expP, "expP") << '\n';
			std::cerr << peri::xyz::infoString(gotP, "gotP") << '\n';
		}

		return errCount;
	}

	//! Check inverse reconstruction recipie
	int
	checkLpa
		( peri::Info const & info
		)
	{
		int errCount{ 0 };

		// assuming this root solution is available
		double const sigma{ info.sigmaValue() };
		double const & x1 = info.theVecX[0];
		double const & x2 = info.theVecX[1];
		double const & x3 = info.theVecX[2];
		double const & mu1 = info.theMus[0];
		double const & mu2 = info.theMus[1];
		double const & mu3 = info.theMus[2];
		double const frac1{ x1 / (peri::sq(mu1) + sigma) };
		double const frac2{ x2 / (peri::sq(mu2) + sigma) };
		double const frac3{ x3 / (peri::sq(mu3) + sigma) };

		double const fSq1{ peri::sq(frac1) };
		double const fSq2{ peri::sq(frac2) };
		double const fSq3{ peri::sq(frac3) };
		double const fhSq{ fSq1 + fSq2 };
		double const fxSq{ fhSq + fSq3 };

		double const lam{ std::atan2(frac2, frac1) };
		double const phi{ std::atan2(frac3, std::sqrt(fhSq)) };
		double const eta{ sigma * std::sqrt(fxSq) };

		peri::LPA const & expLpa = info.theLpa;
		peri::LPA const gotLpa{ lam, phi, eta };

		constexpr double tolAng{ 1.e-14 };
		double const tolLin{ 1.e-14 * expLpa[2] };
		if (! peri::lpa::sameEnough(gotLpa, expLpa, tolAng, tolLin))
		{
			std::cerr << "FAILURE checkLpa() error" << std::endl;
			std::cerr << peri::string::allDigits(expLpa, "expLpa") << '\n';
			std::cerr << peri::string::allDigits(gotLpa, "gotLpa") << '\n';
			std::cerr << peri::string::allDigits(tolAng, "tolAng") << '\n';
			std::cerr << peri::string::allDigits(tolLin, "tolLin") << '\n';
			++errCount;
		}
		else
		{
			std::cout << "Success checkLpa()" << std::endl;
			std::cerr << peri::lpa::infoString(expLpa, "expLpa") << '\n';
			std::cerr << peri::lpa::infoString(gotLpa, "gotLpa") << '\n';
		}

		return errCount;
	}

	//! Check approximations to ellipsoid radius computation
	int
	checkEllipRad
		( peri::Info const & info
		)
	{
		int errCount{ 0 };

		// expected value from exact formula
		peri::XYZ const & xyz = info.theVecX;

		// common terms
		peri::ellip::SquaredParms const sqParms(xyz);
		double const & beta = info.theShape.theRadB;
		double const & hoxSq = sqParms.the_hoxSq;
		double const eSq{ peri::ellip::squaredEccentricity(info.theShape) };
		double const ehoxSq{ eSq * hoxSq };

		// use normalized value for comparison
		double const expRad{ peri::ellip::radiusToward(xyz, info.theShape) };
		double const expApx{ (1./beta) * expRad };

		// check second order approximation
		double const got2nd{ 1 + .5*ehoxSq };
		constexpr double tol2nd{ 5.e-6 };

		// check fourth order approximation
		double const got4th{ 1 + .5*ehoxSq * (1 + (3./4.)*ehoxSq) };
		constexpr double tol4th{ 5.e-8 };

		// check sixth order approximation
		double const got6th
			{ 1 + .5*ehoxSq*(1. + (3./4.)*ehoxSq*(1 + (5./6.)*ehoxSq)) };
		constexpr double tol6th{ 5.e-11};

		std::array<double, 3u> const gots{ got2nd, got4th, got6th };
		constexpr std::array<double, 3u> tols{ tol2nd, tol4th, tol6th };
		for (std::size_t nn{0u} ; nn < gots.size() ; ++nn)
		{
			std::size_t const order{ 2u * nn };
			double const & gotApx = gots[nn];
			double const & tolApx = tols[nn];
			double const difApx{ gotApx - expApx };
			if (! peri::sameEnough(gotApx, expApx, tolApx))
			{
				std::cerr << "FAILURE checkApx() error: order ="
					<< " " << order << '\n';
				std::cerr << peri::string::allDigits(expApx, "expApx") << '\n';
				std::cerr << peri::string::allDigits(gotApx, "gotApx") << '\n';
				std::cerr << peri::string::allDigits(difApx, "difApx") << '\n';
				std::cerr << peri::string::allDigits(tolApx, "tolApx") << '\n';
				++errCount;
			}
			else
			{
				std::cout << "Success checkApx() order = " << order << '\n';
				std::cout << peri::string::allDigits(difApx, "difApx") << '\n';
			}
		}

		return errCount;
	}

}


//! Check various ellipsoid traits and relationships
int
main
	()
{
	int errCount{ 0 };
	// Non-normalized data values to stress formulae
	peri::Shape const shape{ peri::shape::sWGS84/*.normalizedShape()*/ };
	peri::LPA const lpa{ .5*peri::pi(), .25*peri::pi(), 4.e8 };
	peri::Info const info(lpa, shape);
	std::cout << "------" << '\n';
	std::cout << info.infoString("info") << '\n';
	std::cout << "------" << '\n';

	errCount += checkOnEllip(info);
	errCount += checkGrad(info);
	errCount += checkEta(info);
	errCount += checkSigma(info);
	errCount += checkX(info);
	errCount += checkP(info);
	errCount += checkLpa(info);
	errCount += checkEllipRad(info);

	return errCount;
}

