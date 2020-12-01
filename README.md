[]( #  )
[]( #  )
[]( #  MIT License )
[]( #  )
[]( #  Copyright (c) 2020 Stellacore Corporation. )
[]( #  )
[]( #  Permission is hereby granted, free of charge, to any person obtaining )
[]( #  a copy of this software and associated documentation files (the )
[]( #  "Software"), to deal in the Software without restriction, including )
[]( #  without limitation the rights to use, copy, modify, merge, publish, )
[]( #  distribute, sublicense, and/or sell copies of the Software, and to )
[]( #  permit persons to whom the Software is furnished to do so, subject )
[]( #  to the following conditions: )
[]( #  )
[]( #  The above copyright notice and this permission notice shall be )
[]( #  included in all copies or substantial portions of the Software. )
[]( #  )
[]( #  THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY )
[]( #  KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE )
[]( #  WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND )
[]( #  NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS )
[]( #  BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN )
[]( #  AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR )
[]( #  IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN )
[]( #  THE SOFTWARE. )
[]( #  )
[]( #  )

# Peridetic - Near-Surface Geodetic Transformations

Peridetic is an extremely lightweight, easy-to-use, simple C++ header
implementation providing precise, accurate and fast transformations
between Geographic (lon/lat/alt) and Cartesian (x/y/z) coordinate
expressions.

## Peridetic - Key Points

_Easy_: Simplicity, simplicity, simplicity

* Install (or simply clone or copy) the (two) header files, compile and go.

_Freedom_:

* Peridetic is a liberally licensed (MIT/X11) for any uses including
commercial ones.

_Useful_:

Provides the two fundamental and ubiquitous geodetic/rectangular
coordinate transformations:

* Longitude/Parallel(latitude)/(ellipsoidal)Altitude from Cartesian (XYZ)

* Cartesian XYZ coordinates from Longitude/Parallel/Altitude

_Applications_:

* Land/Terrestrial - E.g. infrastructure, vehicles, pedestrians,
mountain climbers, etc.

* Bathymetric - E.g. ships, submarines, undersea cables, buoys, etc.

* Aviation - E.g. aircraft, UAV/drones, balloons, instrumented birds, etc.

_Testing/Verification_:

* Optimal Domain is within +/-100[km] altitude from Earth's (ellipsoid)
surface.

* Precision: Better than ~7.6[nm] within optimal domain (this corresponds
	with computational relative precision of approximately 1.e-15).

* Accuracy: Agrees with select globally distributed NOAA/NGS CORS stations
within their published coordinate precision of 1[mm].

_Limitations and Cautions_:

* Not optimized for use in outer space (although
	internal precision remains under approximately 0.2[um] at lunar
	distances).

* Not optimized for use deep within Earth interior (although
	transformations still meet design precision limits until
	approximately below -5800[km] depths. 

	* Locations involving altitudes below -6300[km], still transform
	with sub[mm] precision (roundtrip consistency better than
	100[um]). However, if this region is important to you, then you
	are working in a extremely specialized and novel domain that is
	not within scope here.
	
* For detail on precision inside and outside the optimal domain
refer to the section on [Transformation Precision](#xform-precision)


## Peridetic - Project

### Project Motivation

GNSS (Global Navigation Satellite System) receivers are everywhere
(e.g. in phones, cars, etc), and existing geographic information is
ubiquitous (e.g. maps, infrastructure databases, etc). Therefore,
software operates with either or both types of data often needs to
convert coordinate locations expressed in one system into equivalent
coordinates expressed in the other.

Geographic transformations provide the means to convert the GNSS Cartesian
coordinate data values (aka "ECEF", "XYZ") into and/or from coordinates
expressed in conventional terrestrial Geographic terms (aka "lon/lat" or
"longitude, latitude, altitude").

Software development projects, even simple ones, often need to perform
these transformations. However, existing transformation options are
typically only available as a small part of a much larger code project.
This can be inconvenient in requiring a full installation processes,
climbing the package learning curve, and, most critically, can create
the need to carry the (potentially large) dependencies forward in the
(often small) project work at hand.

Peridetic addresses this need by offering the two essential geographic
transformations without the overhead of installing or learning a new
environment and without creating any unnecessary dependencies for software
development or runtime environments.

If "All you want to do is convert XYZ into/from LPA values and vice versa",
and you're using C++, then this Peridetic package has been created just
for you.

### Project Concept

Peridetic provides (MIT license) C++ (header) code that performs
precise and accurate transformation from/into geographic coordinates
(longitude, latitude, and ellipsoidal altitude) and Cartesian "XYZ"
(Earth Centered Earth Fixed - ECEF). The transformations are optimized
for locations "near" (within +/- 100[km]) of Earth's ellipsoidal surface .

### Project Name

The project name, Peridetic, is a catenation of "PERI" and "geoDETIC"
and is reflective of geodetic operations designed to be valid within
the important practical domain of operation that is "near" (within
approximately +/- 100[km]) of Earth's surface (ref: [domain of
validity](#domain-of-validity)).

### Project Contributions

Your consideration to use Peridetic is appreciated and respected and your
comments, suggests, ideas and criticism  are all welcome with open arms
feedback is most welcome at the email link below.

The simple "vision" behind Peridetic is "simplicity".

Feedback pertaining to further simplification of the code, its description
and/or project structure are especially welcome (esp to improve documentation).

If you graciously wish to offer specific changes, please fork the
repository and issue a pull request (from a rebased branch).


## Peridetic - Getting Started

To get started refer to:

* Documentation

	* Here

	* TODO - Doxygen (github pages?)

* Installation/Use

	* Quick Hack: "git clone" this repo into a directory that
		is in the compiler include path.

	* A bit more civilized: Copy the two header files
		(peridetic.h and periDetail.h) from here into
		whatever directory you want (e.g. your own project
		or local or system include directories as you like).

	* To build example programs (command line batch coordinate conversions),
		please refer to the related github project TODO - perideticExamples.

	* TODO - packages

		* TODO - header-only install package

		* TODO - examples package

		* TODO - documentation package

		* TODO - development and test environment package

* Usage

	* Simple Examples - ref further below -- TODO

	* Detail Examples - ref Doxygen documents - TODO

	* Full Program Example -- TODO

* Questions and Feedback

	* TODO - FAQ

	* (Email)[mailto://peridetic@stellacore.com]

	* TODO - Feature requests (ref email)



## Peridetic - General Use

Peridetic transformations are very easy to use. Include the source
header file and call one of the transformation functions that are in the
"peri" namespace. Functions use standard C++ structures for argument and
return types. All data values are interpreted consistently in standard
units: _radians_ for angles; _meters_ for distances.

Illustrative code snippet:

	// include header file
    #include "peridetic.h" // indirectly includes implementation periDetail.h

	// standard C++ data types
	//	(generally, just use peri::LPA and peri::XYZ as defined in header)
	using LPA = std::array<double, 3u>;
	using XYZ = std::array<double, 3u>;

	// Start using the transformations (angle *Radians*, linear *meters*)
	// Note: these numeric values only approximately match each other
	//       (ref: Precision further below for full detail on precision)
	LPA const gotLPA{peri::lpaFromXyz(XYZ{-1266326., -4725993., 4877986.})};
	XYZ const gotXYZ{peri::xyzFromLpa(LLA{.6981317, -1.8325957, 0.})};

A complete demonstration examples may be found here

	* TODO - [here](http://example.com).

	* TODO - link to doxygen generated outputs



## Peridetic - Technical Detail

### Terminology

A comment on terminology and notation:

* "XYZ" is used herein to denote classic Cartesian coordinates. In
technical math speak - the three components of a vector as it is expressed
with respect to a orthonormal dextral (right-handed) basis. For code, all
data values (input arguments and return values) are interpreted as _meters_.

* "LPA" is used herein to denote geographic coordinates. The letters stand
for "Longitude", "Parallel" (of latitude), and (ellipsoidal) "Altitude".
The use of the "P" (instead of a second 'L') provides an easy way to
distinguish the two geographic angles in notation.

Note that "Altitude" is used herein to mean "ellipsoidal height" - I.e.
the distance between point of interest and the (outer facing) ellipsoidal
surface.

### Basic Geodesy

This package provides two fundamental transformations:

* LPA coordinates from XYZ coordinates - Geographic from Cartesian

* XYZ coordinates from LPA coordinates - Cartesian from Geographic

The concept of geographic location (the LP parts) is only meaningful if
it is associated with a very specific "Figure of Earth".

For standard geographic coordinate conversions, the Figure of Earth, is
accepted to be an ellipsoid of revolution (specifically an oblate one).
Peridetic header files include definitions for two shapes,
the WGS84 and GRS80 ellipsoids, which are commonly encountered with GNSS
data and modern geodetic activities. Other and custom ellipsoids are
easily created and can be supplied to these transformations.

The transformations offer the WGS84 as a default ellipsoid so that they
are out-of-box compatible with most modern GNSS data.
Note that, in practical terms, the WGS84 ellipsoid is the essentially
same shape as the GRS80 ellipsoid (within 0.1[mm] at the North pole).

### Peridetic - Software Considerations

This Peridetic project (in it's simplest form, is a pair of public/private
C++ header files) provides the two most useful geodetic transformations
that provide high accuracy and precision results and fast runtime 
performance within this operational domain.

#### Use and Integration

This project can be used in your own code in two ways:

* Copy the "peridetic.h" public header file, and the "periDetail.h" 
implementation detail header file wherever you wish, then include
the "peridetic.h" file in your own source code, compile and go.

* -or- TODO - Incorporate this project into another development effort
using the ["cmake" paradigms](https://cmake.org/documentation/)

#### Licensing

This code may be used for any purpose (including commercial) provided
the copy right notice is retained as described in the
[MIT (aka X11) License)](https://directory.fsf.org/wiki/License:X11)

#### Software Environment

Software development points include:

Language

* Standard C++17 or later

	* TODO - check if C++11 works (should)

Compilers

* GCC - primary

* Clang - TODO

Processors

* Precision relies on 8-byte (64-bit) or larger "double" type.

* Tested on

	* x86 - (primary: AMD Ryzen 2700)

Error Handling:

* No exceptions are involved or utilized within the code
	(none used, none thrown)

* Internal code supports quite_NaN propagation
	(e.g. if any input data are NaN, then an all NaN result value is returned)

Thread safety:

* Functions are both re-entrant and thread safe

	* TODO - should be, but needs verification

#### Data Preconditions

No argument value testing is performed since any function arguments
generally have a valid interpretation if they contain non-degenerate
numeric values (i.e.  provided values are not NaN or infinity values -
in which case NaN values are returned).

##### Crazy Inputs

Internal checking of values is limited to conditions that may occur
in the course of processing meaningful data. I.e. there is no checking
for unusual conditions such as providing denormalized values as input
and so on.

There intentionally is no effort (nor overhead) to detect infinity or
subnormal data values. Instead, such values are propagated through
computations based on the properties and characteristics of local
compiler, libraries and hardware.

If this is important, consider wrapping the transforms functions inside
a test along the lines of:

	double const wildValueComponent{ ... };
	if ((0. == wildValueComponent) || std::isnormal(wildValueComponent))
	{
		// Here, wildValueComponent is a valid LPA or XYZ component
	}

NaN values are used internally to represent "null" (i.e. undefined or
uninitialized values). In general, these are not visible to consuming
code. However, if the calling code provides any NaN values as
input data then NaN populated data instances will be returned.

Consuming code can verify that return data values are valid
(which they should always been unless input values are out of control) by
checking for quite_NaN. For invalid return values, if any component is
bad, then *all* components are set to NaN values, so that only the first
array element need be tested. E.g.

	peri::XYZ const ecefLoc{ peri::xyzForLpa(...) };
	if (! std::isnan(ecefLoc[0]))
	{
		// Life is good
	}

##### Data Interpretation and Integrity

Large magnitude input angle values, will be wrapped into principal angle
domains, and interpreted accordingly.

Therefore, beware of the easy potential mistake of providing an angle in
units of (wrong)degrees instead of correct units of radians. For example
(e.g. a numeric value of 90.0 is an angle of ~1.7575743 not the potentially
(incorrectly)expected quarter turn of ~1.5707963 radians).

Overall, consuming code is expected to be nominally responsible
for itself in terms of providing values that are meaningful in a geodetic
sense.

For example, negative altitude values with magnitude greater than
the ellipsoid's polar radius are not valid geographic coordinates.
To avoid the overhead in testing for such a silly condition, Peridetic
assume the consuming code is sufficiently responsible to avoid this.

Overall, as long as the calling code is reasonably responsible with its
own data values, then everything should be fine. If you are uncomfortable
with this level of responsibility, you might consider utilizing a few
utility functions from the project test environment: TODO-

	* TODO - principalValueFor()
	* TODO - isInDomain()
	* TODO - isValid()



## Peridetic - Transformation Details

### Transformations for Location Representations

* Geodetic Surface Location (LPA) from Earth Centered Cartesian (XYZ)

* Earth Centered Cartesian (XYZ) from Geodetic Surface Location (LPA)

The standard LPA and XYZ coordinate representations are both specified
as a triple of numeric values (std::array) but with dramatically
different interpretations. The specific [XYZ](#XYZ-Coordinates) and
[LPA](#LPA-Coordinates) coordinate systems are described in detail
further below.

#### Coordinate Frame Definitions

* Reference Ellipsoids - are used to define the origin and alignment of
the underlying coordinate frames for both XYZ and LPA data values.

* Specific reference ellipsoid may be provided to each transformation. If
none is provided then a current (as of 2020) WGS84 ellipsoid is used.

#### XYZ <a id=#XYZ-Coordinates></a>

The abbreviation "XYZ" is used to denote coordinates in a Cartesian
Coordinate system that are also known commonly as the "ECEF" (Earth
Centered, Earth Fixed) coordinates.

Ref: ["coordinate, Cartesian"](https://www.ngs.noaa.gov/PUBS_LIB/GEOID/Droman_papers/glossary/xml/C.xml)

##### Remarks on XYZ Coordinates

The XYZ (ECEF) Cartesian coordinate frame is often useful when working
with GNSS (Global Navigation Satellite System) observations and/or other
computations that are global in scope.

* Origin -  Is associated with the geometric center of the reference
ellipsoid. For most ellipsoids, this is approximately at the centroid of
Earth's mass distribution.

* Z - is the axis orthogonal to the equator (The plane of rotation
associated to the ellipsoid under consideration). On Earth, points toward
the North pole.

* X - is axis orthogonal to "Z" (in the equatorial plane) and directed
toward the prime meridian. On Earth points approximately toward Gulf
of Guinea.

* Y - is axis mutually orthogonal to Z and X in "right-hand" sense. On
Earth, points approximately south of Bay of Bengal.

The XYZ is a classic orthonormal rectangular coordinate system. It is
associated with three basis vectors (aka "axes"). Each basis vector
has unit magnitude (is "normalized"). The three axes are mutually
perpendicular (are "orthogonal") and have a dextral (aka "right-handed")
chirality interpretation (in the order: "X", "Y", "Z").

Using XYZ coordinates, distances and angles can be computed directly
from the coordinate component values (e.g. via Pythagorean theorem,
the law of cosines, etc).

#### LPA <a id=#LPA-Coordinates></a>

Geodetic Surface Location is expressed by three values denoted as
"LPA", where

* Origin.LP - is on the surface of ellipsoid at a point on its equator
chosen (arbitrarily) to be identified with the prime meridian.

* L(ongitude) - is an azimuthal angle (in radians), positive Eastward
from the prime meridian.

* P(arallel) - (of latitude) is an elevation angle (in radians),
positive Northward from the equator.

* Origin.A - is dynamically defined to lie on the surface of the ellipsoid
at a location specified by the LP coordinates.

* A(ltitude) - is the distance of a point from the Origin.A location.
The "A" value is associated with a distance value that is "stationary"
with respect to small changes in LP coordinates.  In practice, dual
solutions are resolved by selecting the "A" to be the shortest distance
from the ellipsoid surface to the point of interest.

Ref: ["coordinate, geocentric"](https://www.ngs.noaa.gov/PUBS_LIB/GEOID/Droman_papers/glossary/xml/C.xml)

Note that Peridetic uses the term "altitude" to be clear that the
interpretation is applied to an ellipsoid (not a specialized Geoid surface).

##### Remarks on LPA Coordinates

The LPA ellipsoidal coordinate frame is often useful when dealing with
local topocentric applications and often used for navigation on and near
the surface of Earth.

The LPA system is an orthogonal curvilinear coordinate system. Three
unique basis vectors (aka axes) can be defined for each individual point
in space.  However, for these basis directions to be defined uniquely,
requires specification (or assumption) of a specific underlying
ellipsoid. I.e. the LPA coordinate values are a determined _BOTH_ by
the point location _AND_ by the ellipsoid in acting in combination with
each other.

Note that, even for a given fixed ellipsoid, the LPA values are
*_NOT_ unique*. As just one illustrative example, a point with
LPA of (0,0,0) can also be expressed as (pi,0,-2b) where 'b' is an
equatorial radius of the ellipsoid. Every point in space has this dual
LPA representation. Also, the origin (center of Earth) has an infinite
number of LPA representations.

Peridetic transformations are concerned only with the locations on and
"near" the surface of Earth (ref: [domain of validity](#domain-of-validity)).
In this case, the singularity at the center of Earth is mostly irrelevant
(ref: [special cases](#special-cases)).  Of the two dual LPA representations,
only one is in the domain of validity (the one with altitude that has the
smallest absolute value).

### Domain of Validity <a id=domain-of-validity></a>

The quality of the results produced by the code in this project is
associated with a particular domain of validity. Results are fairly useful
outside of this domain, although the quality may be less than when operating
inside the specified domain.

Practical applications are typically concerned with locations and motions
within the atmosphere, at sea, or on the surface of Earth. A useful
(conservative) definition of this region may be arbitrarily specified as
the volume of space within +/- 100[km] of an ellipsoidal shape that best-fits
the Earth surface.

There is nothing magic about this arbitrary +/- 100[km] threshold. However
it is a useful limit at which to evaluate and express transformation
precision and accuracy since it encompasses the overwhelming majority of
practical use-cases. Higher altitudes are conventionally associated with
"outer space" while lower altitudes are physically inaccessible (with
current technology).

There are many subtleties in exactly what a "best-fit" ellipsoid means
when it pertains to the
["Figure of the Earth"](https://en.wikipedia.org/wiki/Figure_of_the_Earth).
However most commonly used ellipsoid definitions have dimensions that
are within a few 100[m] of each other and all of these common ellipsoids
differ from a pure spherical model by less than about +/- 11[km]

## Peridetic - Technical Deep Dive

### Transformation Algorithms

In general, geodetic transformations are non-linear operations for
which computation can be somewhat complex and involved. However, most
practical applications of Geodesy involve spatial locations that are
fairly "near" to Earth's surface (within +/-100[km]). Within this domain,
the transformations can exploit techniques that remain highly accurate
and precise while being simpler to implement and faster to run.

There are a number of papers that espouse the use of "closed form"
solutions for the LPA from XYZ conversion. Certainly this is one approach
that can be taken for implementation. However, these solutions require
evaluation of n-th root computations (i.e. cube roots and square roots)

The root evaluations are not algebraic operations but are themselves
evaluated by iterative algorithms (generally within math library source
code, processor firmware, or in transistor hardware structures).

Peridetic algorithms utilize a simple, fast and precise direct
iterative approach followed by the use of two trig function (std::atan2())
evaluations used to express the solution in angular units to
return the longitude/parallel coordinates as angle values.

The XYZ from LPA algorithm utilizes a single std::sqrt() call.

### Transformation Precision <a id=xform-precision></a>

For purposes here, "precision" is defined loosely as "self-consistency",
"repeatability", "computational significance". It is an intrinsic
metric that describes the quality of the computations more so than
the quality associated with data values.

Peridetic transformations are self-consistent with a precision on
the order of 7.6[nm] for locations within within the operational optimal
design domain (i.e. within +/-100[km] from Earth surface). This level
of precision extends considerably farther although runtime performance
efficiency may drop slightly if operating outside the optimal domain.

The precision (worst case for any coordinate component) changes as
function of point location altitudes (Alt values):

 * Alt < -6300[km]: -- Out of range! Don't do - expect garbage. If you
 	expect to be using data in this range, then put guard code before
	or after calling peri::xyzForLpa() to ensure it is only called or
	results are only used when the altitude value is in range.

 * -6300[km] <= Alt < -5800[km]: -- Reduced precision, < 100[um]

 * -5800[km] <= Alt < +11[Mm]: -- Meets design precision, < 7.6[nm]

 * +11[Mm] <= Alt < +405[Mm]: -- Reduced precision, < 0.2[um]

 * Beyond lunar distances, precision will continue to drop as
	altitude increases.  For example of the extreme: at the altitude
	of the black hole, Sagittarius A-star at the center of the
	Milkyway, transform precision reduces to < 100[km]. On the order
	of 1/4 of the way to Andromeda gallaxy, a 64-bit double completely
	loses all precision for expressing distances relative to the size
	of Earth.

For precision testing, the transformations are evaluated for
self-consistency.  Note that self-consistency does _not_ constitute a
proof of correctness, but only provides a measure of numeric/computational
noise involved.  Therefore, even if transformation computations
are precise, there is the completely independent question if they are
accurate (correct).

The question of accuracy is addressed in the section
[Transform Accuracy](#xform-accuracy) below.

### Transformation Accuracy <a id=xform-accuracy></a>

To determine transformation accuracy it is necessary to compare results
with "absolute" or "known" values.

This current implementation of Peridetic is evaluated using the
published data sources described in the following sections.

#### NOAA NGS CORS network comparisions:

A sampling of published CORS station locations is used to evaluate the
accuracy of the geographic from/into transformations.

The National Geodetic Survey ([NGS](https://www.ngs.noaa.gov/)) is
a branch of the United States National Oceanic and Atmospheric Administration
[NOAA](https://www.noaa.gov/) responsible for:

	"NOAA’s National Geodetic Survey (NGS) provides the framework
	for all positioning activities in the Nation."

The NGS manages the Continuously Operating Reference Stations
([CORS](https://www.ngs.noaa.gov/CORS/)) network of survey locations.
Although concentrated within the U.S. the network includes a number
of stations scattered around the globe.

The following stations were selected with an emphasis on variation 
of geographic location (longitude and latitude) and elevation
(altitude).

	* // PARAKOU (BJPA),  BORGOU
	* // WESTEAST__AK2008 (AC60),  ALASKA
	* // EAST DOCK (EDOC),  ALASKA
	* // AMERICAN SAMOA (ASPA),  AMERICAN SAMOA
	* // FORTALEZA 2005 (BRFT),  UNIDENTIFIED STATE OF BRAZIL
	* // IRAQ SURVY BASRAH (ISBS),  IRAQ
	* // ZANDERIJ (SRZN),  UNIDENTIFIED DISTRICT OF SURINAM
	* // BOULDER (DSRC),  COLORADO
	* // STEAMBOAT SPRINGS (STBT),  COLORADO
	* // MEXICO CITY WAAS (MMX1),  DISTRITO FEDERAL

For each station, NGS publishes both geographic and Cartesian coordinates.
The following is an example of the published data values obtainable
from the interactive [NGS CORS Map](https://geodesy.noaa.gov/CORS_Map/):

	> ITRF2014 POSITION (EPOCH 2010.0)
	> Computed in Jul 2020 using 14 days of data.
	>     X =  -3851330.396 m     latitude    =  52 42 52.63061 N
	>     Y =    399608.571 m     longitude   = 174 04 34.56698 E
	>     Z =   5051382.453 m     ellipsoid height =   18.309   m

These values are used to evaluate Peridetic transforms via a process
that includes:

* Cut-n-paste data from the interactive map site into a header
	file (in the peridetic development and testing codebase).

* Unit test reads data from the header file and populates variables
	with these "expected" data values (after decoding d-m-s values
	into Radians).

* The expected values are transformed in each direction to obtain "got"
	values. E.g.

	* gotXYZ=xyzForLpa(expLPA)

	* gotLPA=lpaForXYZ(expXYZ).

* The "got" values are compared with the corresponding "expected" values
	and differences from zero are compared against tolerance values that
	reflect the provided data precision. I.e.

	* Angular tolerance of { 172. / 1024./1024./1024./1024. };
		This is <.16[nRad] (published values have resolution of
		approximately .049[nRad], but the +/-1[mm] surface distance
		for comparison corresponds with the larger angle ~.16[nRad]
		which is therefore used for testing.

	* Linear tolerance of { 1./1024. }; // CORS files only good to [mm]


* Meet domain design/testing specifications (ref above).

Overall, these transforms are thought to be correct and it is hoped
that the open source transparency of the implementation code and
opportunity for peer review should help provides assurance on its
agreement with external conventions and expectations.

However, please _note_ the _"AS IS"_ clause of the associated license.

### Transformation Quality and Testing

TODO

* Unit tests in ./test module (TODO integrate with project)

* Verification in ./verify module (TODO integrate with project)


