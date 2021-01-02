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

	//! Perturbation 'zeta'-polynomial function
	struct ZetaPolyQuad
	{
		double const theEta0{};
		double const theGrMag{};
		std::array<double, 3u> const theMuSqs{};

		ZetaPolyQuad
			() = default;

		//! Compute constants associated with 'zeta' perturbation polynomials
		struct FSN
		{
			double const theFgkInv{ peri::sNan };
			double const theS1k{ peri::sNan };
			double const theN1k{ peri::sNan };
			double const theN1SqPerMuSq{ peri::sNan };

			inline
			explicit
			FSN
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

			//! Increment associated with coefficient 'Ak'
			inline
			double
			dA
				() const
			{
				return { theN1SqPerMuSq * theS1k * theS1k };
			}

			//! Increment associated with coefficient 'Bk'
			inline
			double
			dB
				() const
			{
				return { theN1SqPerMuSq * theS1k };
			}

			//! Increment associated with coefficient 'Ck'
			inline
			double
			dC
				() const
			{
				return { theN1SqPerMuSq };
			}

		}; // ZetaPolyQuad::FSN

		//! Coefficients 'alpha,betat,gamma' for quadratic zeta polynomial
		inline
		static
		std::array<double, 3u>
		coABCs
			( peri::XYZ const & xVec
			, double const & eta0
			, double const & grMag
			, std::array<double, 3u> const & muSqs
			)
		{
			FSN const fsn0(xVec[0], grMag, muSqs[0], eta0);
			FSN const fsn1(xVec[1], grMag, muSqs[1], eta0);
			FSN const fsn2(xVec[2], grMag, muSqs[2], eta0);
			return std::array<double, 3u>
				{ 3. * (fsn0.dA() + fsn1.dA() + fsn2.dA())
				, fsn0.dB() + fsn1.dB() + fsn2.dB()
				, (fsn0.dC() + fsn1.dC() + fsn2.dC()) - 1.
				};
		}

		//! Coefficients 'alpha,betat,gamma' for quadratic zeta polynomial
		inline
		std::array<double, 3u>
		coABCs
			( peri::XYZ const & xVec
			) const
		{
			return coABCs(xVec, theEta0, theGrMag, theMuSqs);
		}

		//! Root computed exactly (with square root function)
		inline
		double
		rootExact
			( std::array<double, 3u> const & coABCs
			) const
		{
			// shorthand notation
			double const & coA = coABCs[0];
			double const & coB = coABCs[1];
			double const & coC = coABCs[2];
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
		rootApprox1st
			( std::array<double, 3u> const & coABCs
			) const
		{
			// shorthand notation
			double const & coB = coABCs[1];
			double const & coC = coABCs[2];
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
			( std::array<double, 3u> const & coABCs
			) const
		{
			// shorthand notation
			double const & coA = coABCs[0];
			double const & coB = coABCs[1];
			double const & coC = coABCs[2];
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
			( std::array<double, 3u> const & coABCs
			) const
		{
			// shorthand notation
			double const & coA = coABCs[0];
			double const & coB = coABCs[1];
			double const & coC = coABCs[2];
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

		// form and solve 'zeta' quadratic
		ZetaPolyQuad const zpoly{ eta0, grMag, shape.theMuSqs };
		std::array<double, 3u> const coABCs{ zpoly.coABCs(xVec) };
		// double const zeta{ zpoly.rootExact(coABCs) };
		// double const zeta{ zpoly.rootApprox1st(coABCs) };
		double const zeta{ zpoly.rootApprox2nd(coABCs) };
		// double const zeta{ zpoly.rootApprox3rd(coABCs) };

		// compute correction factor for adjusting location coordinates
		double const correction{ 2. * (zeta + eta0) / grMag };

		// compute point on ellipse
		std::array<double, 3u> const & muSqs = shape.theMuSqs;
		XYZ const pVec
			{ xVec[0] / (1. + (correction/muSqs[0]))
			, xVec[1] / (1. + (correction/muSqs[1]))
			, xVec[2] / (1. + (correction/muSqs[2]))
			};

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

	constexpr std::size_t numRad{  32u + 1u };
	constexpr std::size_t numPar{ 256u + 1u }; // odd number hits at 45-deg Lat
//#define UseNorm
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
	peri::sim::SampleSpec const radSpec{ numRad, Range{ radMin, radMax } };
	peri::sim::SampleSpec const parSpec{ numPar, Range{ 0.,  .5*peri::pi()} };

	errCount += test1(radSpec, parSpec, earth, ofsExcess);
	errCount += test2(radSpec, parSpec, earth, ofsDifPVec);

	return errCount;
}

