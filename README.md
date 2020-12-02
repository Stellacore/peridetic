
	MIT License
	Copyright (c) 2020 Stellacore Corporation.


# Peridetic - Easy and Effective Geodetic Transformations

Peridetic is an extremely lightweight, easy-to-use, simple C++ header
implementation providing precise, accurate and fast transformations
between Geodetic (lon/lat/alt) and Cartesian (x/y/z) coordinate
expressions.

Quick Links:

* [Example Usage](#General-Use)
* [Quick Start](#Getting-Started)
* [MIT/X11 License](LICENSE)

General Content:

* [Key Points](#Key-Points)
* [Project Info](#Project-Info)
* [Getting Started](#Getting-Started)
* [General Use](#General-Use)
* [Technical Detail](#Technical-Detail)
* [Transformation Details](#Transformation-Details)
* [Technical Deep Dive](#Technical-Deep-Dive)



## Peridetic - Key Points <a id=Key-Points></a>

#### _Easy_:
Simplicity, simplicity, simplicity, lightweight, small and easy.

* Install (or simply clone or copy) the (two) header files from /include
	directory, compile and go.

#### _Freedom_:

* Peridetic is a liberally licensed (MIT/X11) for any uses including
	commercial ones.

#### _Useful_:

Provides the two most fundamental and ubiquitous geodetic/rectangular
coordinate transformations:

* Longitude/Parallel(latitude)/(ellipsoidal)Altitude from Cartesian (XYZ)

* Cartesian XYZ coordinates from Longitude/Parallel/Altitude

#### _Applications_:

* Land/Terrestrial - E.g. infrastructure, vehicles, pedestrians,
	mountain climbers, etc.

* Bathymetric - E.g. ships, submarines, undersea cables, buoys, etc.

* Aviation - E.g. aircraft, UAV/drones, balloons, instrumented birds, etc.

#### _Testing/Verification_:

* Optimal Domain is within +/-100[km] altitude from Earth's (ellipsoid)
	surface. (ref [Altitude Domain of Validity](#domain-of-validity)

* Precision: Better than ~7.6[nm] within the
	[optimal domain](#domain-of-validity) which corresponds
	with computational relative precision of approximately
	1.e-15 which approaches the limit of [64-bit IEEE-754 'double'
	type](https://en.wikipedia.org/wiki/Decimal64_floating-point_format)
	(with a touch of room for computation noise).

* Accuracy: Agrees with select globally distributed NOAA/NGS CORS stations
	within their published coordinate precision of 1[mm].

#### _Limitations and Cautions_:

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

* This implementation addresses standard geodetic coordinate
	conversions in which longitude/paralell angles are associated
	with an ellipsoidal Figure of Earth. By contrast, astrometric
	lon/par values are referenced to a Geoid (local level surface)
	instead.  Astrometric conversions are outside the scope of these
	Peridetic transformations.

* For detail on precision both inside and outside the optimal domain refer
	to the section on [Transformation Precision](#Transformation-Precision)


## Peridetic - Project Info <a id=Project-Info ></a>

### Project Motivation

GNSS (Global Navigation Satellite System) receivers are everywhere
(e.g. phones, cars, etc), and existing geodetic information is
ubiquitous (e.g. maps, infrastructure databases, etc). Therefore,
software operating with either/both types of data often needs to
convert coordinate locations expressed in one system into an equivalent
coordinates representation in the other.

Geodetic transformations provide the means to convert the Cartesian
coordinate data values (aka "ECEF", "XYZ") such as those from GNSS into
and/or from Geodetic coordinates associated with conventional terrestrial
locations expressed in geographic terms (aka "lon/lat" or "longitude,
latitude, altitude").

Software development projects, even simple ones, often need to perform
these transformations. However, available existing transformation software
is typically only available as a small part of a much larger code project.
This can be inconvenient by requiring a full installation processes,
dealing with the package learning curve, and, most critically, can create
the need to carry (potentially large) dependencies forward in the (often
small) project work at hand.

Peridetic addresses this need by offering the two essential geodetic
transformations without the overhead of installing or learning a new
environment and without creating any unnecessary dependencies for software
development or runtime environments.

If "All you want to do is convert XYZ into/from LPA values and vice versa",
and you're using C++, then this Peridetic package has been created just
for you.

### Project Concept

Peridetic provides (MIT license) C++ (header) code that performs
precise and accurate transformation from/into geodetic coordinates
(longitude, latitude, and ellipsoidal altitude) and Cartesian "XYZ"
(Earth Centered Earth Fixed - ECEF). The transformations are optimized
for locations "near" (within +/- 100[km]) of Earth's ellipsoidal surface.

### Project Name

The project name, Peridetic, is a catenation of "PERI" and "geoDETIC"
and is reflective of geodetic operations designed to be highly
performant within the important practical domain of operation "near"
(within approximately +/- 100[km]) Earth's surface (ref: [domain of
validity](#domain-of-validity)).

### Project Contributions

Thank you for your interest in Perdetic and for considering to provide
comments, suggestions, ideas and/or criticisms - all of which are
welcomed with open arms.

For quick comments, please send a short
[Email](mailto://peridetic@stellacore.com). If you wish to offer more
extensive specific change suggestions, please fork the repository and
issue a pull request from a rebased branch in the fork repo.

The "vision" behind Peridetic is "utility with simplicity"

* Description/Documentation:
	Feedback pertaining to further simplification of the code, its
	description and/or this github project structure are particularly
	welcome (especially to improve documentation).

* Correspoinding Coordinates:
	Another useful contribution would be identification and selection
	of known transformation pair data in the form of corresponding
	coordinates at other stations around the world. E.g. Station
	locations associated with different national systems and reference
	networks (ref [Accuracy section](#Transformation-Accuracy))

* Terminology:
	The Peridetic project (code and documentation) strives to use
	consistent accepted terminology. However, the initial version
	did not focus on this.	Identification of misleading or ambiguous
	terminology is most welcome. As a target, the Peridetic
	project terminology should be consistent with the [NGS
	glossary](https://www.ngs.noaa.gov/CORS-Proxy/Glossary/xml/NGS_Glossary.xml)


## Peridetic - Getting Started <a id=Getting-Started></a>

To get started refer to:

* Documentation:

	* Top-level README.md (this file)

	* TODO - Doxygen (github pages?)

* Installation:

	* Quick Hack: download or cut/paste the two header files
		(peridetic.h and periDetail.h) from this repo's "include"
		directory and put them into whatever directory you
		want (e.g. your own project or a local or system include
		directory as you like).

	* Auxiliary Packages:

		* TODO - header-only install package

		* TODO - examples package

		* TODO - documentation package

		* TODO - development and test environment package

* Usage:

	* Simple Examples - Ref [General Use](#General-Use) section below.

	* Detail Examples - ref Doxygen documents - TODO

	* Full Program Example -- TODO

* Questions and Feedback:

	* TODO - FAQ

	* Comments and questions via [Email](mailto://peridetic@stellacore.com).

	* TODO - Feature requests (ref email)



## Peridetic - General Use <a id=General-Use></a>

Peridetic transformations are very easy to use. Include the source
header file and call one of the transformation functions that are in the
"peri" namespace. Functions use standard C++ structures for argument and
return types. All data values are interpreted consistently in standard
units: _radians_ for angles; _meters_ for distances.

### Illustrative Description: <a id=Illustrative-Example></a>

In succinct terms, the core usage is:

	{
	using namespace peri;
	LPA const gotLPA{ lpaFromXyz( XYZ{ -1266326., -4725993., 4877986. }) };
	XYZ const gotXYZ{ xyzFromLpa( LPA{ .6981317, -1.8325957, 0. }) };
	}

With a bit more explanation:

	// include header file
    #include "peridetic.h" // indirectly includes implementation periDetail.h

	// standard C++ data types
	//	(generally, use peri::LPA and peri::XYZ from header, but for clarity)
	using LPA = std::array<double, 3u>;
	using XYZ = std::array<double, 3u>;

	// Start using the transformations (angle *Radians*, linear *meters*)
	// Note: these numeric values only approximately match each other
	//       (ref: Precision further below for full detail on precision)
	LPA const gotLPA{peri::lpaFromXyz(XYZ{-1266326., -4725993., 4877986.})};
	XYZ const gotXYZ{peri::xyzFromLpa(LPA{.6981317, -1.8325957, 0.})};

### Detailed Example Code

A complete collection of demonstration examples may be found here

	* TODO - [here](http://example.com).

	* TODO - link to doxygen generated outputs



## Peridetic - Technical Detail <a id=Technical-Detail></a>

### Terminology

A comment on terminology and notation:

* "XYZ" is used herein to denote classic Cartesian coordinates. In
	technical math speak - the three components of a vector as it is
	expressed with respect to a orthonormal dextral (right-handed)
	basis. For code, all data values (input arguments and return
	values) are interpreted as _meters_.

* "LPA" is used herein to denote geodetic coordinates. The letters stand
	for "Longitude", "Parallel" (of latitude), and (ellipsoidal)
	"Altitude".  The use of the "P" (instead of a second 'L')
	provides an easy way to distinguish the two geodetic location
	angles in single-letter notation.

Note that "Altitude" is used herein to mean "ellipsoidal height"
- I.e. the directed distance from the surface of the ellipsoid to
the point of interest.

The magnitude of the altitude value may be interpreted as the shortest
distance between the point of interest and any other point on the
ellipsoidal surface (singularities and multiple solution conditions near
Earth center are outside of the [design domain](#domain-of-validity).
The algebraic sign of altitude values is positive for point locations
outside the ellipsoid surface (locally upward) and is negative for points
of interest within the ellipsoid (locally downward)

### Basic Geodesy

The concept of geodetic location (the LP parts of LPA) is only meaningful
when the angles are interpreted in association with a specific
"Figure of Earth".

For standard geodetic coordinate conversions, the Figure of Earth, is
accepted to be an ellipsoid of revolution (specifically an oblate one
with equatorial radius larger then polar radius).

Peridetic header files include definitions for two shapes, the WGS84
and GRS80 ellipsoids, which are commonly encountered with GNSS data
and modern geographic data uses. Other standard and/or entirely custom
ellipsoids are easily created and can be used directly with the
transformation functions (TODO - ref peri::model namespace).

The transformations offer the WGS84 as a default ellipsoid so that
transformations are out-of-box compatible with most modern GNSS data.
Note that, in applications terms, the WGS84 ellipsoid is, for virtually
all practical purposes, simply a different way of specifying the same
shape as does the GRS80 ellipsoid (to within 0.1[mm] at the North pole).

### Peridetic - Software Considerations <a id=Software-Considerations></a>

At it's core, this Peridetic project comprises two source code header files:

* peridetic.h -- public interface (includes periDetail.h)

* periDetail.h -- implementation (inline functions) code

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

#### Data Preconditions

Function arguments have a generally valid interpretation if they contain
non-degenerate numeric values (i.e.  provided values are not NaN or
infinity values. Therefore, no input argument value testing is performed.

However, if input data contain any NaN values, those are propagated
through computations, then _all components_ of the return array are set
to NaN.

##### Crazy Values

Internal checking of values is limited to conditions that may occur
in the course of processing meaningful data.

There intentionally is no code (nor overhead) to test for infinity or
subnormal data values. Instead, such values are propagated through
computations based on the properties and characteristics of local
compiler, libraries and hardware.

If an ability to provide bad data values is important, then consider
wrapping the transforms functions inside a test along the lines of:

	double const wildValueComponent{ ... };
	// test for each component that might contain bad data
	if ((0. == wildValueComponent) || std::isnormal(wildValueComponent))
	{
		// Here, wildValueComponent is a valid LP or XYZ component
		// Altitude should still satisfy (MinAlt < Alt)
	}

NaN values are used internally to represent "null" (i.e. undefined or
uninitialized values). In general, these are not visible to consuming
code. However, if the calling code provides any NaN values as input
data then NaN populated data instances will be returned. Consuming code,
if desired, can rely on this for use as part of a "null object pattern"
paradigm.

E.g. consumer code can verify that return data values are valid (which they
should always been unless input values are out of control) by checking for
quite_NaN. For invalid return values, if any component is bad, then *all*
components are set to NaN values, so that only the first array element
need be tested. E.g.

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
(e.g. a numeric value of 90.0 is interpreted as an angle of ~1.7575743 not
the potentially (incorrectly)expected quarter turn of ~1.5707963 radians).

Overall, consuming code is expected to be nominally responsible
for itself in terms of providing values that are meaningful in a geodetic
sense.

For example, negative altitude values with magnitude greater than
the ellipsoid's polar radius are not valid geodetic coordinates.
To avoid the overhead in testing for such a silly condition, Peridetic
assumes the consuming code is sufficiently responsible to avoid this.

Overall, as long as the calling code is reasonably responsible handling its
own data values, then everything should be fine. If you are uncomfortable
with this level of responsibility, you might consider utilizing a few
utility functions from the project test environment:

	* TODO- periLocal.h dev/test package function access



## Peridetic - Transformation Details <a id=Transformation-Details></a>

### Transformations for Location Representations

* Geodetic Surface Location (LPA) from Earth Centered Cartesian (XYZ)

* Earth Centered Cartesian (XYZ) from Geodetic Surface Location (LPA)

The standard LPA and XYZ coordinate representations are both specified
as a triple of numeric values (std::array) but with dramatically
different interpretations. The specific interpretation of the
[XYZ](#XYZ-Coordinates) and [LPA](#LPA-Coordinates) coordinate systems
are described in detail further below.

It is important to note that the conversion between XYZ/LPA representations
is inextricably associated with a specific model for the Figure of Earth.
For Peridetic (as for geodetic coordinates in general), the Earth shape
model is an ellipsoid. (Note: astrometric lon/par values are referenced
to Geoid instead of an ellipsoid. Astrometric conversions are outside 
the scope of these Peridetic transformations).

#### Reference Ellipsoids

Reference Ellipsoids - are used to define the origin and alignment of
the underlying coordinate frames for both XYZ and LPA data values.  
Specific reference ellipsoid may be provided to each transformation. If
none is provided then the WGS84 ellipsoid is used.

There are many subtleties in exactly what constitutes a
"best-fit" ellipsoid means when it pertains to the
["Figure of the Earth"](https://en.wikipedia.org/wiki/Figure_of_the_Earth).
However most commonly used ellipsoid definitions have dimensions that
are within a few 100[m] of each other and the commonly used ellipsoids
differ from a pure spherical model by less than about +/- 11[km] between
equatorial and polar axes.

Peridetic provides two common Earth model specifications:

* GRS80 - ref peri::model::GRS80

* WGS84 - ref peri::model::WGS84

* Other standard as well as custom ellipsoids may be
	created easily via peri::Shape and peri::EarthModel data
	structures. Ref advanced documentation - TODO.

#### LPA <a id=LPA-Coordinates></a>

Geodetic Surface Location is expressed by three values denoted as
"LPA", where

* Origin.LP - is on the surface of ellipsoid at a point on its equator
	chosen (arbitrarily) to be identified with the prime meridian
	of particular convention (e.g. Greenwich, Paris, Mecca, [and
	others](https://en.wikipedia.org/wiki/Prime_meridian#List_of_prime_meridians_on_Earth)).

* L(ongitude) - is an azimuthal angle (in radians), positive Eastward
	from the prime meridian.

* P(arallel) - (of latitude) is an elevation angle (in radians),
	positive Northward from the equator.

* Origin.A - is dynamically defined to lie on the surface of the ellipsoid
	at a location specified by the LP coordinates.

* A(ltitude) - is the distance of a point from the Origin.A location.
	The "A" value is associated with a distance value that is
	"stationary" with respect to small changes in LP coordinates.
	In practice, dual solutions are resolved by selecting the "A"
	to be the shortest distance from the ellipsoid surface to the
	point of interest.

For explanation of Longitude and Latitude angles, ref: the NGS glossary
entry for
["coordinate, geocentric"](https://www.ngs.noaa.gov/PUBS_LIB/GEOID/Droman_papers/glossary/xml/C.xml).

Note that Peridetic uses the term "parallel" in place of latitude
(to facility unique naming notation) and uses the term "altitude"
to be clear that the interpretation is applied to an ellipsoid (in
effort to reduce confusion associated with the many uses of "height"
and "elevation").

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

#### XYZ <a id=XYZ-Coordinates></a>

The abbreviation "XYZ" is used to denote coordinates in a Cartesian
Coordinate system that are also known commonly as the "ECEF" (Earth
Centered, Earth Fixed) coordinates.

Ref: ["coordinate, Cartesian"](https://www.ngs.noaa.gov/PUBS_LIB/GEOID/Droman_papers/glossary/xml/C.xml)

##### Remarks on XYZ Coordinates

The XYZ (ECEF) Cartesian coordinate frame is often useful when working
with GNSS (Global Navigation Satellite System) observations and/or other
computations that are global in scope.

* Origin -  Is associated with the geometric center of the reference
	ellipsoid. For most ellipsoids, this is approximately at the
	centroid of Earth's mass distribution.

* Z - is the axis orthogonal to the equator (The rotation plane of
	symmetry associated with the ellipsoid under consideration). On
	Earth, points toward the North pole.

* X - is axis orthogonal to "Z" (in the equatorial plane) and directed
	toward the prime meridian. On Earth, for the Greenwich prime
	meridian, this points approximately toward Gulf of Guinea.

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

### Domain of Validity on Altitude <a id=domain-of-validity></a>

The quality of the results produced by the code in this project is
associated with a particular domain of validity. Results remain fairly useful
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
"outer space" while lower altitudes are physically inaccessible with
current technology.

The [Transformation Precision](#Transformation-Precision) section describes
what to expect for transformation of locations outside this optimum
domain, 

## Peridetic - Technical Deep Dive <a id=Technical-Deep-Dive></a>

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

Peridetic algorithms utilize a simple, fast and precise direct iterative
approach followed by the use of two trig function evaluations (of
std::atan2()) that are used to express the solution in angular units
for return of longitude/parallel coordinates as conventional angle values.

The XYZ from LPA algorithm utilizes a single std::sqrt() call.

The mathematical formulae and algorithm detail is described in:

* PerideticMath.{lyx,pdf} -- TODO

### Transformation Precision <a id=Transformation-Precision></a>

For purposes here, "precision" is defined loosely as "self-consistency",
"repeatability", "computational significance". As such, it is an intrinsic
metric that describes the quality of the computations more so than
the quality associated with data values.

Peridetic transformations are self-consistent with a precision on the
order of 7.6[nm] for locations within the operational optimal design
domain (i.e. within +/-100[km] from Earth surface). This level of
precision extends considerably farther in altitude although runtime
performance efficiency may drop slightly if operating outside the
optimal domain.

The precision (worst case for any coordinate component) changes as
function of point location altitudes (Alt values):

 * Alt < -6300[km]: -- Out of range! Don't do - expect garbage.
 
	* If you expect data in this range, then your application
		is outside the scope of validity for using Peridetic.

	* If you have a legitimate use-case for needing to do this, please
		describe briefly in an email, and the author will buy you a
		beverage.

 * -6300[km] <= Alt < -5800[km]: -- Reduced precision, < 100[um]

 * -5800[km] <= Alt < +11000[km]: -- Meets *design precision*, < 7.6[nm]

	* **-100[km] <= Alt <= 100[km]: -- is optimal design domain**

 * +11000[km] <= Alt < +405[Mm]: -- Reduced precision, < 0.2[um]

 * Beyond lunar distances, precision will continue to drop as
	altitude increases.  As an extreme example: at the altitude
	of the black hole, Sagittarius A-star, at the center of the
	Milkyway, the transform precision reduces to <100[km]. On the order
	of 1/4 of the way to Andromeda gallaxy, a 64-bit double completely
	loses all precision for expressing distances relative to the size
	of Earth.

The above precision estimates are created by testing round trip
transformations for various point locations well distributed with respect
ellipsoid. For testing within the optimal domain, results at approximately
1M point locations are evaluated.

For precision testing, the transformations are evaluated for
self-consistency.  Note that self-consistency does _not_ constitute a
proof of correctness, but only provides a measure of numeric/computational
noise involved.  Therefore, even if transformation computations are
precise, there is a completely independent question concerning how
accurate (correct) are the results.

The question of accuracy is addressed in the section
[Transform Accuracy](#Transformation-Accuracy) below.

### Transformation Accuracy <a id=Transformation-Accuracy></a>

To determine transformation accuracy it is necessary to compare results
with "absolute" or "known" values available from external sources.

This current implementation of Peridetic is evaluated using the
published data sources described in the following sections.

#### NOAA/NGS/CORS network comparisons:

A sampling of published CORS geodetic network reference station
locations is used to evaluate the accuracy of the Peridetic
coordinate expression transformations.

The National Geodetic Survey ([NGS](https://www.ngs.noaa.gov/)) is
a branch of the United States National Oceanic and Atmospheric Administration
[NOAA](https://www.noaa.gov/):

	"NOAA’s National Geodetic Survey (NGS) provides the framework
	for all positioning activities in the Nation."

The NGS manages the Continuously Operating Reference Stations
([CORS](https://www.ngs.noaa.gov/CORS/)) network of survey locations.
Although concentrated within the U.S. the network includes a number
of stations scattered around the globe.

The following stations were selected with an emphasis on variation 
of geodetic location (longitude and latitude) and elevation
(altitude).

	// PARAKOU (BJPA),  BORGOU
	// WESTEAST__AK2008 (AC60),  ALASKA
	// EAST DOCK (EDOC),  ALASKA
	// AMERICAN SAMOA (ASPA),  AMERICAN SAMOA
	// FORTALEZA 2005 (BRFT),  UNIDENTIFIED STATE OF BRAZIL
	// IRAQ SURVY BASRAH (ISBS),  IRAQ
	// ZANDERIJ (SRZN),  UNIDENTIFIED DISTRICT OF SURINAM
	// BOULDER (DSRC),  COLORADO
	// STEAMBOAT SPRINGS (STBT),  COLORADO
	// MEXICO CITY WAAS (MMX1),  DISTRITO FEDERAL

For each station, NGS publishes both Geodetic and Cartesian coordinates.
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
	and differences from zero are compared against tolerance values. The
	tolerance values used reflect the published data precision:

	* Angular tolerance of { 172. / 1024./1024./1024./1024. };
		This is <.16[nRad] (published angular values have resolution of
		approximately .049[nRad], but the published XYZ coordinates
		only have +/-1[mm] precision for comparison. A surface distance
		of 1[mm] corresponds with a larger angle ~.16[nRad] which is
		therefore used for testing.

	* Linear tolerance of { 1./1024. }; // CORS files XYZ published to [mm]

* All pairs of coordinates are required to pass this tolerance test (in
	both directions (lpaforXyz and xyzforLpa).

Overall, the Peridetic transforms are thought to be correct as described
above. The open source transparency of the algorithm and implementation
code along with opportunity for associated peer review should help
provides assurance on the quality. However...

...Please _note_ the _"AS IS"_ clause of the associated license.

### Transformation Quality and Testing

TODO

* Unit tests in ./test module (TODO integrate with project)

* Verification in ./verify module (TODO integrate with project)


