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


#ifndef cors_DataParser_INCL_
#define cors_DataParser_INCL_


#include "peridetic.h"

#include "periLocal.h"

#include <sstream>
#include <string>
#include <vector>


namespace peri
{
namespace cors
{

	//! Parse cut/paste strings into usable data
	struct DataParser
	{
		//! Cartesian XYZ ECEF values
		XYZ const theXYZ{ xyz::sNull };

		//! Geodetic LPA values
		LPA const theLPA{ lpa::sNull };

		//! Value of (+/-1.) depending on direction (W:-1,E:+1,S:-1,N:+1)
		static
		double
		pmOneForDir
			( std::string::value_type const & dirWENS
			)
		{
			double pmOne{ sNan };
			switch (dirWENS)
			{
				case 'W':
				case 'S':
					pmOne = -1.;
					break;
				case 'E':
				case 'N':
					pmOne = +1.;
					break;
			}
			return pmOne;
		}

		//! Extract values from text - assuming precisely above formatting
		static
		DataParser
		from
			( std::string const & text
			)
		{
			XYZ xyz{ xyz::sNull };
			LPA lpa{ lpa::sNull };

			// Format cut/paste from NGS reports (sans vertical bars)
			/*
			{
			" ITRF2014 POSITION (EPOCH 2010.0)   "
			" Faked in Nov 2020 using madeup data CORS reporting style.  "
			"     X =   1000000.001 m     latitude    =  01 22 33.00004 N  "
			"     Y =   2000000.002 m     longitude   = 023 44 55.00006 E  "
			"     Z =   3000000.003 m     ellipsoid height = 4000.004   m  "
			};
			*/

			// find start of data portion of text
			std::string::size_type const pos{ text.find("X =") };
			if (! (std::string::npos == pos))
			{
				std::istringstream iss(text.substr(pos));
				std::string na{};
				double locX{}, locY{}, locZ{};
				double degP{}, minP{}, secP{};
				std::string::value_type dirSN{};
				double degL{}, minL{}, secL{};
				std::string::value_type dirWE{};
				double alt{};
				iss
					>> na >> na >> locX >> na
					>> na >> na >> degP >> minP >> secP >> dirSN
					>> na >> na >> locY >> na
					>> na >> na >> degL >> minL >> secL >> dirWE
					>> na >> na >> locZ >> na
					>> na >> na >> na >> alt >> na
					;
				if (! iss.fail())
				{
					double const signL{ pmOneForDir(dirWE) };
					double const signP{ pmOneForDir(dirSN) };
					XYZ const tmpXYZ
						{ locX, locY, locZ };
					LPA const tmpLPA
						{ signL*radMagFromDMS(degL, minL, secL)
						, signP*radMagFromDMS(degP, minP, secP)
						, alt
						};
					if (peri::isValid(tmpXYZ) && peri::isValid(tmpLPA))
					{
						xyz = tmpXYZ;
						lpa = tmpLPA;
					}
				}

				/*
				std::cout << std::endl;
				std::cout << std::fixed << std::setprecision(6);
				std::cout << "degP: " << degP << std::endl;
				std::cout << "minP: " << minP << std::endl;
				std::cout << "secP: " << secP << std::endl;
				std::cout << "degL: " << degL << std::endl;
				std::cout << "minL: " << minL << std::endl;
				std::cout << "secL: " << secL << std::endl;
				std::cout << " alt: " << alt << std::endl;
				std::cout << xyz::infoString(xyz, "xyz") << std::endl;
				std::cout << lpa::infoString(lpa, "lpa") << std::endl;
				std::cout << std::endl;
				*/

			}
			return DataParser{ xyz, lpa };
		}

		//! True if all members have valid content
		bool
		isValid
			() const
		{
			return
				(  peri::isValid(theXYZ)
				&& peri::isValid(theLPA)
				);
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
				oss << title << '\n';
			}
			if (isValid())
			{
				oss << xyz::infoString(theXYZ, "theXYZ");
				oss << '\n';
				oss << lpa::infoString(theLPA, "theLPA");
			}
			else
			{
				oss << " <null>";
			}
			return oss.str();
		}

	}; // DataParser


} // [cors]
} // [peri]

// Definitions for functions
// #include "corsDataParser.inl"


#endif // cors_DataParser_INCL_

