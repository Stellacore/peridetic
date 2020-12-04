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


/*! \brief XYZ/LPA coordinate pairs from US NOAA/NGS CORS stations
 *
 * Coordinates cut/pasted from site information obtained via
 * \arg https://geodesy.noaa.gov/CORS_Map/
 *
 */


#include <string>
#include <vector>

#include <sstream>


namespace peri::cors
{
/* -- Note from bottom of coordinate descriptions
 * Latitude, longitude and ellipsoid height are computed from their
   corresponding cartesian coordinates using dimensions for the
   GRS 80 ellipsoid: semi-major axis = 6,378,137.0 meters
                          flattening = 1/298.257222101...
*/

// PARAKOU (BJPA),  BORGOU
static std::string const sTextBJPA
{
" ITRF2014 POSITION (EPOCH 2010.0)                                            "
" Computed in Jun 2019 using data through gpswk 1933.                         "
"     X =   6287630.499 m     latitude    =  09 21 27.03254 N                 "
"     Y =    288340.869 m     longitude   = 002 37 32.35800 E                 "
"     Z =   1030266.214 m     ellipsoid height =  423.917   m                 "
};

// WESTEAST__AK2008 (AC60),  ALASKA
static std::string const sTextAC60
{
" ITRF2014 POSITION (EPOCH 2010.0)                                            "
" Computed in Jul 2020 using 14 days of data.                                 "
"     X =  -3851330.396 m     latitude    =  52 42 52.63061 N                 "
"     Y =    399608.571 m     longitude   = 174 04 34.56698 E                 "
"     Z =   5051382.453 m     ellipsoid height =   18.309   m                 "
};

// EAST DOCK (EDOC),  ALASKA
static std::string const sTextEDOC
{
" ITRF2014 POSITION (EPOCH 2010.0)                                            "
" Computed in Jun 2019 using data through gpswk 1933.                         "
"     X =  -1834182.999 m     latitude    =  70 18 36.97329 N                 "
"     Y =  -1131997.250 m     longitude   = 148 19 06.54187 W                 "
"     Z =   5982812.006 m     ellipsoid height =   22.414   m                 "
};

// AMERICAN SAMOA (ASPA),  AMERICAN SAMOA
static std::string const sTextASPA
{
" ITRF2014 POSITION (EPOCH 2010.0)                                            "
" Computed in Jun 2019 using data through gpswk 1933.                         "
"     X =  -6100260.080 m     latitude    =  14 19 33.93679 S                 "
"     Y =   -996503.222 m     longitude   = 170 43 20.77040 W                 "
"     Z =  -1567977.593 m     ellipsoid height =   53.552   m                 "
};

// FORTALEZA 2005 (BRFT),  UNIDENTIFIED STATE OF BRAZIL
static std::string const sTextBRFT
{
" ITRF2014 POSITION (EPOCH 2010.0)                                            "
" Published by the IGS in May 2020.                                           "
"     X =   4985393.521 m     latitude    =  03 52 38.80678 S                 "
"     Y =  -3954993.444 m     longitude   = 038 25 31.93497 W                 "
"     Z =   -428426.657 m     ellipsoid height =   21.679   m                 "
};

// IRAQ SURVY BASRAH (ISBS),  IRAQ
static std::string const sTextISBS
{
" ITRF2014 POSITION (EPOCH 2010.0)                                            "
" Computed in Jun 2019 using data through gpswk 1933.                         "
"     X =   3695492.629 m     latitude    =  30 29 13.24384 N                 "
"     Y =   4074925.387 m     longitude   = 047 47 44.00284 E                 "
"     Z =   3217012.649 m     ellipsoid height =   -2.376   m                 "
};

// ZANDERIJ (SRZN),  UNIDENTIFIED DISTRICT OF SURINAM
static std::string const sTextSRZN
{
" ITRF2014 POSITION (EPOCH 2010.0)                                            "
" Computed in Jun 2019 using data through gpswk 1933.                         "
"     X =   3623419.977 m     latitude    =  05 27 20.31676 N                 "
"     Y =  -5214015.457 m     longitude   = 055 12 11.07483 W                 "
"     Z =    602359.251 m     ellipsoid height =  -17.251   m                 "
};

	// Altitudes

// BOULDER (DSRC),  COLORADO
static std::string const sTextDSRC
{
" ITRF2014 POSITION (EPOCH 2010.0)                                            "
" Computed in Jun 2019 using data through gpswk 1933.                         "
"     X =  -1288338.794 m     latitude    =  39 59 29.15054 N                 "
"     Y =  -4721988.544 m     longitude   = 105 15 39.72063 W                 "
"     Z =   4078321.096 m     ellipsoid height = 1656.287   m                 "
};

// STEAMBOAT SPRINGS (STBT),  COLORADO
static std::string const sTextSTBT
{
" ITRF2014 POSITION (EPOCH 2010.0)                                            "
" Computed in Aug 2019 using 256 days of data.                                "
"     X =  -1409244.102 m     latitude    =  40 30 45.79500 N                 "
"     Y =  -4648589.132 m     longitude   = 106 51 53.83536 W                 "
"     Z =   4122789.882 m     ellipsoid height = 2087.327   m                 "
};

// MEXICO CITY WAAS (MMX1),  DISTRITO FEDERAL
static std::string const sTextMMX1
{
" ITRF2014 POSITION (EPOCH 2010.0)                                            "
" Computed in Oct 2020 using 14 days of data.                                 "
"     X =   -948701.102 m     latitude    =  19 25 53.95219 N                 "
"     Y =  -5943935.692 m     longitude   = 099 04 06.20239 W                 "
"     Z =   2109212.719 m     ellipsoid height = 2235.680   m                 "
};


	//! Collection of text strings to use for coordinate testing
	static std::vector<std::string> const sStationTexts
		{ sTextBJPA
		, sTextAC60
		, sTextEDOC
		, sTextASPA
		, sTextBRFT
		, sTextISBS
		, sTextSRZN
		, sTextDSRC
		, sTextSTBT
		, sTextMMX1
		};

} // [peri::cors]


