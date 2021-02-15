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


#include <iostream>


//! Supporting evaluations
namespace eqns
{

	//! Precomputed constants associated with a particular ellipsoid shape
	struct EllipConsts
	{
		std::array<double, 3u> const theMus{};
		double const theRho{};
		std::array<double, 3u> const theMuSqs{};
		std::array<double, 3u> const theUpsilons{};

		EllipConsts
			() = default;

		explicit
		EllipConsts
			( peri::Shape const & shape
			)
			: theMus{ shape.theRadA, shape.theRadA, shape.theRadB }
			, theRho{ std::sqrt(theMus[0] * theMus[2]) }
			, theMuSqs{ shape.theMuSqs }
			, theUpsilons
				{ peri::sq(theRho / theMuSqs[0])
				, peri::sq(theRho / theMuSqs[1])
				, peri::sq(theRho / theMuSqs[2])
				}
		{ }

	}; // EllipConsts

	//! Array of element-by-element squared values
	std::array<double, 3u>
	squaredValues
		( peri::XYZ const & xvec
		)
	{
		return
			{ peri::sq(xvec[0])
			, peri::sq(xvec[1])
			, peri::sq(xvec[2])
			};
	}

	//! Initial value for pseudo altitude parameter, "sigma"
	double
	initSigma0
		( std::array<double, 3u> const & xSqs
		, EllipConsts const & ecs
		)
	{
		double sumX2{ 0. };
		double sumUpX2{ 0. };
		for (std::size_t nn{0u} ; nn < 3u ; ++nn)
		{
			sumX2 += xSqs[nn];
			sumUpX2 += ecs.theUpsilons[nn] * xSqs[nn];
		}
		double const num{ sumX2 - ecs.theRho * std::sqrt(sumX2) };
		double const den{ std::sqrt(sumUpX2) };
		double const sig0{ num / den };
		return sig0;
	}


} // [eqns]


//! Q&D program to spot-check formulae presented in doc/perideticSummary.pdf
int
main
	()
{
	peri::EarthModel const earth{ peri::model::WGS84 };

	// use (working) peridetic code as reference for document content check
	// (i.e. summary document written after code implementation verification)
	peri::LPA const periLPA{ 2., .75, 12345. };
	peri::XYZ const periXYZ{ peri::xyzForLpa(periLPA, earth) };

	// report expected values (assuming verified implementation"
	std::cout << '\n';
	std::cout << peri::lpa::infoString(periLPA, "periLPA") << '\n';
	std::cout << peri::xyz::infoString(periXYZ, "periXYZ") << '\n';

	// precomputed constants and short-hand notations
	eqns::EllipConsts const ecs(earth.theEllip.theShapeNorm);
	std::array<double, 3u> const & mu1s = ecs.theMus; // ellip radii
	std::array<double, 3u> const & mu2s = ecs.theMuSqs; // squared ellip radii

	// formulae apply to normalized coordinate data
	peri::XYZ const xvec{ earth.theEllip.xyzNormFrom(periXYZ) };

	// squared coordinate values for a couple spots
	std::array<double, 3u> const xSqs{ eqns::squaredValues(xvec) };

	// iterative solution
	double currSigma{ eqns::initSigma0(xSqs, ecs) };
	double prevSigma{ peri::sNan };
	double bestSigma{ peri::sNan };
	for (std::size_t nn{0u} ; nn < 8u ; ++nn)
	{
		// intermediate values
		std::array<double, 3u> const xis
			{ mu1s[0] / (mu2s[0] + currSigma)
			, mu1s[1] / (mu2s[1] + currSigma)
			, mu1s[2] / (mu2s[2] + currSigma)
			};
		std::array<double, 3u> const psis
			{ peri::sq(xis[0]) * xSqs[0]
			, peri::sq(xis[1]) * xSqs[1]
			, peri::sq(xis[2]) * xSqs[2]
			};

		// apply linearized iteration formula
		double const num{ psis[0] + psis[1] + psis[2] - 1. };
		double const den{ xis[0]*psis[0] + xis[1]*psis[1] + xis[2]*psis[2] };
		prevSigma = currSigma;
		currSigma = prevSigma + (.5 * num / den);

		// parms for evaluating/reporting convergence
		constexpr double tol{ 1.e-15 };
		double const prevEval{ 1. + prevSigma };
		double const currEval{ 1. + currSigma };
		double const diffEval{ currEval - prevEval};

		// watch convergence properties
		std::cout << '\n';
		std::cout
			<< " " << peri::string::allDigits(currSigma, "currSigma")
			<< " " << peri::string::allDigits(prevSigma, "prevSigma")
			<< '\n';
		std::cout
			<< " " << peri::string::allDigits(currEval, "currEval")
			<< " " << peri::string::allDigits(prevEval, "prevEval")
			<< " " << peri::string::allDigits(diffEval, "diffEval")
			<< '\n';

		// test for convergence
		if (std::abs(currEval - prevEval) < tol)
		{
			bestSigma = currSigma;
			break;
		}
	}

	// Compute fraction values shared between various geodetic expressions
	std::array<double, 3u> const fracs
		{ xvec[0] / (mu2s[0] + bestSigma)
		, xvec[1] / (mu2s[1] + bestSigma)
		, xvec[2] / (mu2s[2] + bestSigma)
		};

	// compute the geodetic representation parameters
	using peri::sq;
	double const altNorm
		{ bestSigma * std::sqrt(sq(fracs[0]) + sq(fracs[1]) + sq(fracs[2])) };
	// restore linear units
	double const altOrig{ altNorm * earth.theEllip.lambdaOrig() };

	double const lonOrig
		{ std::atan2(fracs[1], fracs[0]) };
	double const parOrig
		{ std::atan2(fracs[2], std::sqrt(sq(fracs[1]) + sq(fracs[0]))) };

	/*
	std::cout << peri::string::fixedLinear(altOrig, "altOrig") << '\n';
	std::cout << peri::string::fixedAngular(lonOrig, "lonOrig") << '\n';
	std::cout << peri::string::fixedAngular(parOrig, "parOrig") << '\n';
	*/

	peri::LPA const gotLPA{ lonOrig, parOrig, altOrig };

	std::cout << '\n';
	std::cout << peri::lpa::infoString(gotLPA, "gotLPA") << '\n';
	std::cout << std::endl;

}

