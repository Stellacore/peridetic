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
#include "periSim.h"

#include <algorithm>
#include <fstream>
#include <iostream>
#include <utility>
#include <vector>


namespace
{
	//! Compute constants associated with 'zeta' perturbation polynomials
	struct ZetaFSN
	{
		double const theFgkInv{ peri::sNan };
		double const theS1k{ peri::sNan };
		double const theN1k{ peri::sNan };
		double const theN1SqPerMuSq{ peri::sNan };

		ZetaFSN
			() = default;

		inline
		explicit
		ZetaFSN
			( double const & xk
			, double const & grMag
			, double const & muSqk
			, double const & eta0
			)
			: theFgkInv{ .5 * grMag * muSqk }
			, theS1k{ 1. / (theFgkInv + eta0) }
			, theN1k{ theS1k * theFgkInv * xk }
			, theN1SqPerMuSq{ peri::sq(theN1k) / muSqk }
		{
		}

		//! Increment associated with coefficient 'Ck'
		inline
		double
		kN1MuS0
			() const
		{
			return { theN1SqPerMuSq };
		}

		//! Increment associated with coefficient 'Bk'
		inline
		double
		kN1MuS1
			() const
		{
			return { theN1SqPerMuSq * theS1k };
		}

		//! Increment associated with coefficient 'Ak'
		inline
		double
		kN1MuS2
			() const
		{
			return { theN1SqPerMuSq * theS1k * theS1k };
		}

	}; // ZetaPolyQuad::ZetaFSN

	//! Perturbation 'zeta'-polynomial (for testing: is NOT performant)
	struct ZetaPoly
	{
		peri::XYZ const theVecX{};
		double const theEta0{};
		double const theGrMag{};
		std::array<double, 3u> const theMuSqs{};
		std::array<ZetaFSN, 3u> const theFSNs{};

		explicit
		ZetaPoly
			( peri::XYZ const & xVec
			, double const & eta0
			, double const & grMag
			, std::array<double, 3u> const & muSqs
			)
			: theVecX{ xVec }
			, theEta0{ eta0 }
			, theGrMag{ grMag }
			, theMuSqs{ muSqs }
			, theFSNs
				{ ZetaFSN(xVec[0], theGrMag, theMuSqs[0], theEta0)
				, ZetaFSN(xVec[1], theGrMag, theMuSqs[1], theEta0)
				, ZetaFSN(xVec[2], theGrMag, theMuSqs[2], theEta0)
				}
		{ }

		// Constant Term
		//! Sum of (n1k/muk)^2
		inline
		double
		sumN1MuS0
			() const
		{
			return
				( theFSNs[0].kN1MuS0()
				+ theFSNs[1].kN1MuS0()
				+ theFSNs[2].kN1MuS0()
				);
		}

		// Linear Term
		//! Sum of (n1k/muk)^2*s1k
		inline
		double
		sumN1MuS1
			() const
		{
			return
				( theFSNs[0].kN1MuS1()
				+ theFSNs[1].kN1MuS1()
				+ theFSNs[2].kN1MuS1()
				);
		}

		// Quadratic Term
		//! Sum of (n1k/muk)^2*s1k*s1k
		inline
		double
		sumN1MuS2
			() const
		{
			return
				( theFSNs[0].kN1MuS2()
				+ theFSNs[1].kN1MuS2()
				+ theFSNs[2].kN1MuS2()
				);
		}

	}; // ZetaPoly

	//! Perturbation 'zeta'-polynomial expanded to only 1st order (~10^-8)
	struct ZetaPolyLinear : public ZetaPoly
	{
		//! Value construction (for underlying class)
		explicit
		ZetaPolyLinear
			( peri::XYZ const & xVec
			, double const & eta0
			, double const & grMag
			, std::array<double, 3u> const & muSqs
			)
			: ZetaPoly(xVec, eta0, grMag, muSqs)
		{ }

		//! Coefficients 'alpha,betat,gamma' for quadratic zeta polynomial
		inline
		std::array<double, 2u>
		coCBs
			() const
		{
			return std::array<double, 2u>
				{ (sumN1MuS0() - 1.) // constant
				, sumN1MuS1() // linear
				};
		}

		//! Root computed exactly (linear transform)
		inline
		double
		rootExact
			( std::array<double, 2u> const & coCBs
			) const
		{
			// shorthand notation
			double const & coC = coCBs[0];
			double const & coB = coCBs[1];
			// xArg is a relatively small quantity (used in series expansion)
			double const fracCoB{ coC / coB };
			double const zetaExact
				{ .5 * fracCoB };
			return zetaExact;
		}

	}; // ZetaPolyLinear

	//! Perturbation 'zeta'-polynomial function to quadratic order (~10^-16)
	struct ZetaPolyQuad : public ZetaPoly
	{
		//! Value construction (for underlying class)
		explicit
		ZetaPolyQuad
			( peri::XYZ const & xVec
			, double const & eta0
			, double const & grMag
			, std::array<double, 3u> const & muSqs
			)
			: ZetaPoly(xVec, eta0, grMag, muSqs)
		{ }

		//! Coefficients 'alpha,betat,gamma' for quadratic zeta polynomial
		inline
		std::array<double, 3u>
		coCBAs
			() const
		{
			return std::array<double, 3u>
				{ (sumN1MuS0() - 1.) // constant
				, sumN1MuS1() // linear
				, 3.*sumN1MuS2() // quadratic
				};
		}

		//! Root computed exactly (with square root function)
		inline
		double
		rootExact
			( std::array<double, 3u> const & coCBAs
			) const
		{
			// shorthand notation
			double const & coC = coCBAs[0];
			double const & coB = coCBAs[1];
			double const & coA = coCBAs[2];
			// xArg is a relatively small quantity (used in series expansion)
			double const fracCoB{ coC / coB };
			double const xArg{ fracCoB * (coA / coB) };
			double const radical{ std::sqrt(1. - xArg) };
			double const zetaExact
				{ (coB / coA) * (1. - radical) };
			return zetaExact;
		}

		//! Root computed using linear expansion of sqrt()
		inline
		double
		rootApprox1st // Note: same as ZetaPolyLinear::rootExact()
			( std::array<double, 3u> const & coCBAs
			) const
		{
			// shorthand notation
			double const & coC = coCBAs[0];
			double const & coB = coCBAs[1];
			double const fracCoB{ coC / coB };
			// first order approximation
			// precise to about 1.e-7 [m] (between +/- 1000[km])
			double const zetaApx1
				{ (.5 * fracCoB) };
			return zetaApx1;
		}

		//! Root computed using quadratic expansion of sqrt()
		inline
		double
		rootApprox2nd
			( std::array<double, 3u> const & coCBAs
			) const
		{
			// shorthand notation
			double const & coC = coCBAs[0];
			double const & coB = coCBAs[1];
			double const & coA = coCBAs[2];
			// xArg is a relatively small quantity (used in series expansion)
			double const fracCoB{ coC / coB };
			double const xArg{ fracCoB * (coA / coB) };
			// second order approximation
			// precise to about 1.e-8 [m] (between +/- 1000[km])
			double const zetaApx2
				{ (.5 * fracCoB) * (1 + (1./4.)*xArg) };
			return zetaApx2;
		}

		//! Root computed using cubic expansion of sqrt()
		inline
		double
		rootApprox3rd // NOTE: results insignificantly different from 2nd order
			( std::array<double, 3u> const & coCBAs
			) const
		{
			// shorthand notation
			double const & coC = coCBAs[0];
			double const & coB = coCBAs[1];
			double const & coA = coCBAs[2];
			// xArg is a relatively small quantity (used in series expansion)
			double const fracCoB{ coC / coB };
			double const xArg{ fracCoB * (coA / coB) };
			// third order approximation
			// precision UNmeasureable
			// second order approx is already in computation noise
			double const zetaApx3
				{ (.5 * fracCoB) * (1. + ((1./4.) + (1./8.)*xArg)*xArg) };
			return zetaApx3;
		}

	}; // ZetaPolyQuad


	int
	test1
		( peri::sim::SampleSpec const & radSpec
		, peri::sim::SampleSpec const & parSpec
		, peri::EarthModel const & earth
		, std::ostream & ofsExcess
		)
	{
		int errCount{ 0u };

		peri::Shape const & shape = earth.theEllip.theShapeOrig;

		using namespace  peri::sim;
		std::vector<peri::XYZ> const xyzs
			{ meridianPlaneSamples(radSpec, parSpec) };

		std::vector<double> extras{};
		extras.reserve(xyzs.size());
		for (peri::XYZ const & xyz : xyzs)
		{
			peri::XYZ const & xVec = xyz;
			peri::LPA const xLpa{ peri::lpaForXyz(xVec, earth) };
			double const & eta = xLpa[2];

			peri::XYZ const xDir{ peri::unit(xVec) };

			peri::LPA const pLpa{ xLpa[0], xLpa[1], 0. };
			peri::XYZ const pVec{ peri::xyzForLpa(pLpa, earth) };

			double const rMag{ peri::ellip::radiusToward(xDir, shape) };
			using peri::operator*;
			peri::XYZ const rVec{ rMag * xDir };

			using peri::operator-;
			peri::XYZ const xrVec{ xVec - rVec };
			peri::XYZ const xpVec{ xVec - pVec };

			double const xrMag{ peri::magnitude(xrVec) };
			double const xpMag{ peri::magnitude(xpVec) };

			double const extra{ xrMag - xpMag };
			extras.emplace_back(extra);

			// gradients
			peri::XYZ const gpVec{ shape.gradientAt(pVec) };
			double const gpMag{ peri::magnitude(gpVec) };
			peri::XYZ const grVec{ shape.gradientAt(rVec) };
			double const grMag{ peri::magnitude(grVec) };
			double const gRatio{ grMag / gpMag };
			double const rEps{ gRatio - 1. };
		//	peri::XYZ const gDifVec{ grVec - gpVec };
		//	double const gDifMag{ peri::magnitude(gDifVec) };

			// compute deltaEta, deltaSigma
			double const delEta{ xrMag - xpMag };
			double const dEtaPerR{ delEta / rMag };

			// info for analysis
			ofsExcess
			//	<< " xVec: " << peri::xyz::infoString(xVec)
			//	<< " xLpa: " << peri::lpa::infoString(xLpa)
			//	<< " Lon: " << peri::string::fixedAngular(xLpa[0])
				<< " Par: " << peri::string::fixedAngular(xLpa[1])
				<< " Alt: " << peri::string::allDigits(eta)
			//	<< " " << peri::string::fixedLinear(xrMag, "xrMag")
			//	<< " " << peri::string::fixedLinear(xpMag, "xpMag")
				<< " " << peri::string::fixedLinear(extra, "extra")
				<< " " << peri::string::allDigits(rEps, "rEps")
				<< " " << peri::string::allDigits(delEta, "delEta")
				<< " " << peri::string::allDigits(dEtaPerR, "dEtaPerR")
				<< std::endl;
		}

		if (! extras.empty())
		{
			std::sort(extras.begin(), extras.end());
			using peri::string::fixedLinear;
			ofsExcess << "# minExcess: " << fixedLinear(extras.front()) << '\n';
			ofsExcess << "# maxExcess: " << fixedLinear(extras.back()) << '\n';
		}

		// needs work to define what is a reasonable test
		// ++errCount;

		return errCount;
	}

	//! Formula to evaluate vector, p, from 'zeta' root
	peri::XYZ
	pVecFor
		( peri::XYZ const & xVec
		, double const & zeta
		, double const & eta0
		, double const & grMag
		, peri::Shape const & shape
		)
	{
		std::array<double, 3u> const & muSqs = shape.theMuSqs;
		// compute point on ellipse
		double const correction{ 2. * (zeta + eta0) / grMag };
		peri::XYZ const pVec
			{ xVec[0] / (1. + (correction/muSqs[0]))
			, xVec[1] / (1. + (correction/muSqs[1]))
			, xVec[2] / (1. + (correction/muSqs[2]))
			};
		return pVec;
	}

	//! Compute point on ellipsoid, p, using perturbation expansion
	peri::XYZ
	pVecViaExcess
		( peri::XYZ const & xVec
		, peri::EarthModel const & earth
		)
	{
		using namespace peri;
		Shape const & shape = earth.theEllip.theShapeOrig;

		// radial point on ellipsoid
		double const rho{ peri::ellip::radiusToward(xVec, shape) };
		XYZ const rVec{ rho * unit(xVec) };

		// gradient at radial point
		XYZ const grVec{ shape.gradientAt(rVec) };
		double const grMag{ magnitude(grVec) };

		// radial pseudo-altitude
		double const xMag{ magnitude(xVec) };
		double const rMag{ magnitude(rVec) };
		double const eta0{ xMag - rMag };

		// form and solve 'zeta' linear
		ZetaPolyLinear const zpoly1{ xVec, eta0, grMag, shape.theMuSqs };
		std::array<double, 2u> const coCBs{ zpoly1.coCBs() };
		double const zetaLineExact{ zpoly1.rootExact(coCBs) };

		// form and solve 'zeta' quadratic
		ZetaPolyQuad const zpoly2{ xVec, eta0, grMag, shape.theMuSqs };
		std::array<double, 3u> const coCBAs{ zpoly2.coCBAs() };
		// double const zetaQuadExact{ zpoly2.rootExact(coCBAs) };
		double const zetaQuadApx1{ zpoly2.rootApprox1st(coCBAs) };
		double const zetaQuadApx2{ zpoly2.rootApprox2nd(coCBAs) };
		// double const zetaQuadApx3{ zpoly2.rootApprox3rd(coCBAs) };
		double const zeta{ zetaQuadApx2 };

		// 1st order approx to quadratic poly == exact soln to linear poly
		double const difOrder1{ zetaLineExact - zetaQuadApx1 };
		if (! (0. == difOrder1))
		{
			std::cerr << "Implementation error (LineExact != QuadApprox1st)"
				<< std::endl;
			std::cerr << peri::string::allDigits(difOrder1, "difOrder1")
				<< std::endl;
			exit(8);
		}

		// compute point on ellipse
		peri::XYZ const pVec{ pVecFor(xVec, zeta, eta0, grMag, shape) };

		return pVec;
	}

	//! Evaluate math equations for ellipsoidal excess at this point
	int
	checkXYZ
		( peri::XYZ const & xVecExp
		, peri::EarthModel const & earth
		, std::ofstream & ostrm
		)
	{
		int errCount{ 0u };
		using namespace peri;

		// compute pLocation based on perturbation expansion
		XYZ const pVecGot{ pVecViaExcess(xVecExp, earth) };

		// quantities for checking values
		LPA const xLpaExp{ lpaForXyz(xVecExp, earth) };
		XYZ const pLpaExp{ xLpaExp[0], xLpaExp[1], 0. };
		XYZ const pVecExp{ xyzForLpa(pLpaExp, earth) };

		// error amount
		XYZ const pVecDif{ pVecGot - pVecExp };

		double const pMagDif{ magnitude(pVecDif) };

		// check that precision of computation
		double const pMagExp{ magnitude(pVecExp) };
		double const pMagTol{ pMagExp * 1.e-15 };
		if (! (pMagDif < pMagTol))
		{
			std::cerr << "FAILURE of precision test on pMagDif" << '\n';
			std::cerr << xyz::infoString(pVecExp, "pVecExp") << '\n';
			std::cerr << xyz::infoString(pVecGot, "pVecGot") << '\n';
			std::cerr << string::allDigits(pVecDif, "pVecDif") << '\n';
			std::cerr << string::allDigits(pMagDif, "pMagDif") << '\n';
			std::cerr << string::allDigits(pMagTol, "pMagTol") << '\n';
			++errCount;
		}

		// save evaluation data
		ostrm
			<< lpa::infoString(xLpaExp, "xLpaExp")
			<< " "
			<< xyz::infoString(pVecDif, "pVecDif")
			<< " "
			<< string::allDigits(pMagDif, "pMagDif")
			<< '\n';

		return errCount;
	}

	//! Check equations on sampling of points
	int
	test2
		( peri::sim::SampleSpec const & radSpec
		, peri::sim::SampleSpec const & parSpec
		, peri::EarthModel const & earth
		, std::ofstream & ofsDifPVec
		)
	{
		int errCount{ 0u };

		using namespace  peri::sim;
		std::vector<peri::XYZ> const xyzs
			{ meridianPlaneSamples(radSpec, parSpec) };
		for (peri::XYZ const & xyz : xyzs)
		{
			errCount += checkXYZ(xyz, earth, ofsDifPVec);
		}
		return errCount;

	}

	peri::XYZ
	rVecFor
		( peri::XYZ const & xVec
		, peri::Shape const & shape
		)
	{
		// radial point on ellipsoid
		double const rho{ peri::ellip::radiusToward(xVec, shape) };
		using peri::operator*;
		peri::XYZ const rVec{ rho * peri::unit(xVec) };
		return rVec;
	}

	double
	grMagFor
		( peri::XYZ const & xVec
		, peri::Shape const & shape
		)
	{
		// gradient at radial point
		peri::XYZ const grVec{ shape.gradientAt(rVecFor(xVec, shape)) };
		double const grMag{ peri::magnitude(grVec) };
		return grMag;
	}

	double
	eta0For
		( peri::XYZ const & xVec
		, peri::Shape const & shape
		)
	{
		// radial point on ellipsoid
		peri::XYZ const rVec{ rVecFor(xVec, shape) };
		// radial pseudo-altitude
		double const xMag{ peri::magnitude(xVec) };
		double const rMag{ peri::magnitude(rVec) };
		double const eta0{ xMag - rMag };
		return eta0;
	}

	void
	saveZetaDiffs
		( peri::sim::SampleSpec const & radSpec
		, peri::sim::SampleSpec const & parSpec
		, peri::EarthModel const & earth
		, std::ostream & ofsZetaDiff
		)
	{
		peri::Shape const & shape = earth.theEllip.theShapeOrig;

		std::vector<peri::XYZ> const xyzs
			{ peri::sim::meridianPlaneSamples(radSpec, parSpec) };
		for (peri::XYZ const & xyz : xyzs)
		{
			peri::XYZ const & xVec = xyz;

			// gradient at radial point
			double const grMag{ grMagFor(xVec, shape) };

			// radial pseudo-altitude
			double const eta0{ eta0For(xVec, shape) };

			// form and solve 'zeta' quadratic
			ZetaPolyQuad const zpoly2{ xVec, eta0, grMag, shape.theMuSqs };
			std::array<double, 3u> const coCBAs{ zpoly2.coCBAs() };
			double const zetaQuadApx1{ zpoly2.rootApprox1st(coCBAs) };
			double const zetaQuadApx2{ zpoly2.rootApprox2nd(coCBAs) };

			// expected POE
			peri::XYZ const pVecExp{ earth.nearEllipsoidPointFor(xVec) };

			// computed POE with different approximation orders
			peri::XYZ const pVecApx1
				{ pVecFor(xVec, zetaQuadApx1, eta0, grMag, shape) };
			peri::XYZ const pVecApx2
				{ pVecFor(xVec, zetaQuadApx2, eta0, grMag, shape) };

			using peri::operator-;
			double const difMagApx1{ peri::magnitude(pVecApx1 - pVecExp) };
			double const difMagApx2{ peri::magnitude(pVecApx2 - pVecExp) };

			double const xMag{ peri::magnitude(xVec) };
			double const & zVal = xVec[2];
			double const hSq{ peri::sq(xVec[0]) + peri::sq(xVec[1]) };
			double const hVal{ std::sqrt(hSq) };
			double const theta{ std::atan2(zVal, hVal) };
			ofsZetaDiff
				<< peri::string::allDigits(theta, "theta")
				<< " "
				<< peri::string::allDigits(xMag, "xMag")
				<< " "
				<< peri::string::allDigits(difMagApx1, "difMagApx1")
				<< " "
				<< peri::string::allDigits(difMagApx2, "difMagApx2")
			//	<< " "
			//	<< peri::string::allDigits(xVec, "xVec")
				<< '\n';
		}

	}

} // [annon]


//! TODO
int
main
	()
{
	int errCount{ 0 };

	// optionally save data to files for analysis
	std::ofstream ofsExcess
		(
		"/dev/null"
	//	"/dev/stdout"
		);
	std::ofstream ofsDifPVec
		(
		"/dev/null"
	//	"pvecDiff.dat"
		);
	std::ofstream ofsZeta
		(
	//	"/dev/null"
		"zetaDiff.dat"
		);

	constexpr std::size_t numRad{  32u + 1u };
	constexpr std::size_t numPar{ 256u + 1u }; // odd number hits at 45-deg Lat
#define UseNorm
#if defined(UseNorm)
	peri::Shape const shape(peri::shape::sWGS84.normalizedShape());
	constexpr double altLo{ -(100./6370.) };
	constexpr double altHi{  (100./6370.) };
#else
	peri::Shape const shape(peri::shape::sWGS84);
	constexpr double altLo{ -100. * 1.e+3 };
	constexpr double altHi{  100. * 1.e+3 };
#endif

	peri::EarthModel const earth(shape);
	peri::Ellipsoid const & ellip = earth.theEllip;

	double const radEarth{ ellip.lambda() };
	double const radMin{ radEarth + altLo };
	double const radMax{ radEarth + altHi };
	using Range = std::pair<double, double>;
	double const halfPi{ .5 * peri::pi() };
	peri::sim::SampleSpec const radSpec{ numRad, Range{ radMin, radMax } };
	peri::sim::SampleSpec const parSpec{ numPar, Range{ -halfPi, halfPi } };

	errCount += test1(radSpec, parSpec, earth, ofsExcess);
	errCount += test2(radSpec, parSpec, earth, ofsDifPVec);

	saveZetaDiffs(radSpec, parSpec, earth, ofsZeta);

	return errCount;
}

