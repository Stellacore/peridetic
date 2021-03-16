
	MIT License
	Copyright (c) 2020 Stellacore Corporation.


# Peridetic - Easy and Effective Geodetic Transformations

Peridetic is an extremely lightweight, easy-to-use, simple C++
header file implementation providing precise, accurate and fast
transformations between Geodetic (lon/lat/alt) and Cartesian (x/y/z)
coordinate expressions.

Quick Links:

* [Example Usage](#General-Use)
* [Quick Start](#Getting-Started)
* [MIT/X11 License](LICENSE)
* [Context and Alternatives](#Context-and-Alternatives)

This Page Content:

* [Key Points](#Key-Points)
* [Project Info](#Project-Info)
* [Getting Started](#Getting-Started)
* [General Use](#General-Use)
	* [Using Installed Peridetic with CMake](#Using-Installed)
* [Technical Detail](#Technical-Detail)
	* [Software/Dev Considerations](#Software-Considerations)
* [Transformation Details](#Transformation-Details)
* [Technical Deep Dive](#Technical-Deep-Dive) -- math description in [./doc/perideticSummary.pdf](https://github.com/Stellacore/peridetic/tree/main/doc/perideticSummary.pdf)
* [Runtime Performance](#Runtime-Performance)



## Peridetic - Key Points <a id=Key-Points></a>

#### _Easy_:
Simplicity, simplicity, simplicity, lightweight, small and easy.

* Install (or simply copy the header file pair from /include directory),
	compile and go. (ref [Getting Started](#Getting-Started))

#### _Lightweight_:
Tiny code size with [extremely low memory/data usage](#Code-Overhead)
(on the order of 5kB).

* Ideal for low-power and low-performance devices (e.g. phones, raspberry-pi,
custom ARM processors, etc).

#### _Freedom_:

* Peridetic uses the permissive [MIT/X11 License](LICENSE) supporting
	any use cases including commercial ones.

#### _Useful_:

Provides the two most fundamental and ubiquitous geodetic/rectangular
coordinate transformations (ref [examples](#Definitive-Example-Code))
which can be used to obtain:

* Longitude/Latitude/(ellipsoidal)Altitude from Cartesian (XYZ)

* Cartesian XYZ coordinates from Longitude/Latitude/Altitude

#### _Applications_:

* Land/Terrestrial - E.g. infrastructure, vehicles, pedestrians,
	mountain climbing, etc.

* Bathymetric - E.g. ships, submarines, undersea cables, buoys, etc.

* Aviation - E.g. aircraft, UAV/drones, balloons, instrumented birds, etc.

* Space - E.g. near-Earth and geosynchronous satellites, but *slightly* less
	optimal performance.

#### _Testing/Verification_:

* Optimal domain is within +/-100[km] altitude from Earth's (ellipsoid)
	surface. (ref [Altitude Optimal Domain](#optimal-domain))

* [Precision](#Transformation-Precision):
	Better than ~7.6[nm] within the
	[optimal domain](#optimal-domain) which corresponds
	with computational relative precision of approximately
	1.e-15 which approaches the limit of [64-bit IEEE-754 'double'
	type](https://en.wikipedia.org/wiki/Decimal64_floating-point_format)
	(with a touch of room for computation noise).

* [Accuracy](#Transformation-Accuracy):
	Agrees with select globally distributed NOAA/NGS CORS stations
	within their published coordinate precision of 1[mm].

* [Speed](#Runtime-Performance):
	Comparable to fastest available algorithms in terms of floating point
	operations counts.

#### _Limitations and Cautions_:

* Less optimized for use in outer space (although
	internal precision remains under approximately 0.2[um] at lunar
	distances).

* Less optimized for use deep within Earth interior (although
	transformations still meet design precision limits until
	below approximately -5800[km] depths. 

	* Locations until altitudes below -6300[km], still transform
		with sub[mm] precision (roundtrip consistency better than
		100[um]). However, if depths below this region are important
		to you, then you are working in a extremely specialized
		and novel domain that is not within scope here.

* This implementation addresses standard *Geodetic* coordinate
	conversions in which longitude/paralell angles are associated
	to an ellipsoidal Figure of Earth. By contrast, astrodetic
	lon/par values are referenced to a Geoid (local level surface)
	instead.  Astrometric conversions are outside the scope of these
	Peridetic transformations.

	* [This web page](https://www.oc.nps.edu/oc2902w/coord/coord.pdf)
		provides a fairly clear and accessible description of the
		various types of coordinate systems and distinctions between
		them.

* The section, [Transformation Precision](#Transformation-Precision),
	provides detail on precision both inside and outside the optimal
	altitude domain.



## Peridetic - Project Info <a id=Project-Info ></a>

### Project Motivation

GNSS (Global Navigation Satellite System) receivers are everywhere
(e.g. phones, cars, etc), and existing geodetic information is
ubiquitous (e.g. maps, infrastructure databases, etc). Therefore,
software operating with either/both generally needs to
convert coordinate locations expressed in one system into an equivalent
representation in the other domain.

Geodetic transformations provide the means to convert the Cartesian
coordinate data values (aka "ECEF", "XYZ") into/from Geodetic coordinates
associated with conventional terrestrial locations expressed in geographic
terms (aka "lon/lat" or "longitude, latitude, altitude").

Software development projects, even simple ones, often need to perform
these transformations. However, available existing transformation software
is typically only available as a small portion of some much larger project.
This can be inconvenient in requiring a full installation processes,
dealing with the large package learning curve, and, most critically,
creates a need to carry (potentially large) dependencies forward in the
(often small) project work otherwise at hand.

Peridetic addresses this need by offering the two core geodetic
transformations without overhead of installing or learning a new
environment and without creating any unnecessary dependencies for software
development or runtime environments.

If "All you want to do is convert XYZ into/from LPA values and vice versa",
and you're using C++, then this Peridetic package has been created just
for you.
	
If you need more than these two most-simple of transformations, consider
[more extensive full-featured alternatives](#Context-and-Alternatives)

### Project Concept

Peridetic provides (MIT license) C++ header code that performs very
precise and accurate transformation from/into geodetic coordinates
(longitude, latitude, and ellipsoidal altitude) and Cartesian "XYZ" (Earth
Centered Earth Fixed - ECEF). The transformation performance is optimized
for locations "near" (within +/- 100[km]) of Earth's ellipsoidal surface.

### Project Name

The project name, Peridetic, is a catenation of "PERI" and [geo]"DETIC".
The components are associated with Greek words "peri", meaning near
and around (to Earth's surface), and "daiesthai" meaning "to divide"
(into measurable units). The name Peredetic is intended to be refelctive
of geodetic operations that are highly performant within the important
practical domain of operation within approximately +/- 100[km] Earth's
surface (ref: [optimal domain](#optimal-domain)).

### Project Contributions

If you're interested to contribute...

First of all, thank you for your interest in Perdetic and for considering
to provide comments, suggestions, ideas and/or criticisms - all of which
are welcomed with open arms.

For quick comments, you can send a short
[Email](mailto://peridetic@stellacore.com). If you wish to offer more
extensive and/or specific change suggestions, please fork the repository
and issue a pull request from a rebased branch in the fork repo.

The vision behind Peridetic is "utility with simplicity". Candidate areas
to consider contribution include:

* Description/Documentation:
	Feedback pertaining to further simplification of the code, its
	description and/or this github project structure are particularly
	welcome (especially to improve documentation). Perhaps starting
	with this "README" page.

* Correspoinding Coordinates:
	Another useful contribution would be identification and selection
	of known transformation pair data in the form of corresponding
	coordinates at other stations around the world. E.g. Station
	locations associated with different national systems and reference
	networks (ref [Accuracy section](#Transformation-Accuracy))
	and/or high-accuracy cooresponding coordinate pairs associated
	with other ellipsoid shapes.

* Terminology:
	The Peridetic project (code and documentation) strives to
	use consistent accepted terminology. However, some initial
	work did not focus on this and may be a touch sloppy with
	terminology.  Identification of any misleading or ambiguous
	terminology is most welcome. As a target, the Peridetic
	project terminology should be consistent with the
[NGS glossary](https://www.ngs.noaa.gov/CORS-Proxy/Glossary/xml/NGS_Glossary.xml)



## Peridetic - Getting Started <a id=Getting-Started></a>

To get started immediately:

* Grab and use the pair of header files - ref: [Grab-n-Go](#Grab-n-Go)
* -or- Build and install library [more formally](#CMake-Build)
* Use in your application as described in:
	* [Quick Illustrative Description](#Illustrative-Example) - nominal
	* [Definitive Example Code](#Definitive-Example-Code) - definitive

Documentation sources include:

* Top-level README.md (this file)


* Project documentation (doxygen pages) can be generated from a cloned
	repository copy (ref [CMake Build](#CMake-Build) section.
	
E.g. to only generate documentation pages:

	$ cd /tmp
	$ git clone https://github.com/Stellacore/peridetic.git
	$ mkdir /tmp/perideticBuild
	$ cd /tmp/perideticBuild
	$ cmake /tmp/peridetic
	$ make docs
	$ <favoriteBrowser> /tmp/perideticBuild/doc/html/index.html

### Installation

#### Grab-n-Go <a id=Grab-n-Go></a>

Simply download or cut/paste the two header files
[peridetic.h](https://github.com/Stellacore/peridetic/blob/main/include/peridetic.h)
and
[periDetail.h](https://github.com/Stellacore/peridetic/blob/main/include/periDetail.h)
		from this repo's
		[/include](https://github.com/Stellacore/peridetic/tree/main/include)
		directory and put them into your development environment wherever you
		want (e.g. your own project or a local or system include
		directory as you like).

#### CMake Build <a id=CMake-Build></a>

In general, for a formal development process, clone this project then
build and use in one of two ways:

* Library build: Since this is a header only implementation
	the 'build' process is not necessary in order to use the
	transformations. However it does produces the software API and
	implementation documentation (via doxygen) as well as the test
	and verification programs.

* Build a formal software distribution package, then install it onto local
	development environments (and/or into containers) for general use.

E.g. to build Perdietic in /tmp:

	$ mkdir /tmp/perideticWorkArea  # or wherevever you like
	$ cd /tmp/perideticWorkArea
	$ git clone https://github.com/Stellacore/peridetic.git
	$ mkdir tmpBuild  # i.e. /tmp/perideticWorkArea/tmpBuild
	$ cd tmpBuild
	$ cmake ../peridetic -DCMAKE_BUILD_TYPE=Release
		# -DCMAKE_INSTALL_PREFIX=/tmp # e.g. install location (here /tmp)
		# for other than default behavior
		# optionally add command line specifications
		# or edit CMakeCache file (e.g. cmake-gui) as you like
	$ make -j `nproc` # builds documentation with doxygen (also eval/tests)
		# provides peridetic{Targets,Config}.cmake files for cmake
	$ ctest # run test and verification programs
	$ cpack # creates install packages
		# e.g.: peridetic-*-Linux.deb

Then continue with platform/package specific below

#### For Debian packages (e.g. on Ubuntu/Debian linux)

By default (i.e. unless you provide cmake command options or changed
content of the CMakeCache.txt file), installation will include:

* /usr/local/include/peridetic/ -- containing the core header files
	(i.e. all that is needed for use in your application)

* /usr/local/share/doc/peridetic/ -- containing package documentation
	(i.e. API and software reference documentation). E.g. point brower
	to local html index:
	* /usr/local/share/doc/peridetic/html/index.html

* /usr/local/lib/cmake/peridetic/ -- containing cmake project files
	(e.g. to intergrate with other development work using cmake)

On any system that supports Debian formatted packages, continue with e.g.

	$ sudo apt-get install ./peridetic-<VER>-Linux.deb
		# use version ID in place of '<VER>'

Cleanup

	$ cd /tmp  # or somewhere else
	$ rm -rf /tmp/perideticWorkArea

To uninstall (e.g. some time later)

	$ sudu apt-get remove peridetic


Then proceed with use in your own work.

* Usage:

	* Simple Illustrative Examples -- description in text below.
		Ref [General Use](#General-Use) section below.

	* Definitive Examples -- stand-alone programs from project
		"/examples" directory.
		Ref [Detail Examples](#Definitive-Example-Code)

	* This project's [./periUse sub directory](#Using-Installed)
		contains an example mini-project (hello world style) using
		Peridetic as part of an external project built with CMake.

* Questions and Feedback:

	* Comments/requests and questions
		via [Email](mailto://peridetic@stellacore.com).


### Using Installed Peridetic <a id=Using-Installed></a>

For a concrete example of building an independent/stand-alone consuming
project against an already installed Peridetic, refer to the
[periUse sub directory](./periUse/README.md).
The "hello world" style  periUse 'mini-project' example also demonstrates
the use of CMake "find-package()" command to locate an already-installed
Peridetic resource.



## Peridetic - General Use <a id=General-Use></a>

Peridetic transformations are very easy to use. Include the source
header file and call one or both transformation functions that are in the
"peri" namespace. Functions use standard C++ structures for argument and
return types. All data values are interpreted consistently in standard
units: _radians_ for angles; _meters_ for distances.

### Illustrative Description: <a id=Illustrative-Example></a>

NOTE: this section contains generally illustrative text (missing a
few things and possibly including typos).
For definitive, compilable, and executable code, refer to the source file
links in [Definitive Example Code](#Definitive-Example-Code) section.

Using data type aliases from peridetic.h the core usage is:

	// include header for function declarations and definitions
    #include "peridetic.h" // indirectly includes implementation periDetail.h

	// Geodetic from Cartesian  (using type aliases for std::array returns)
	peri::XYZ const haveXyz{ -1266643.136, -4727176.539, 4079014.032 };
	peri::LPA const wantLPA{ peri::lpaForXyz(haveXyz) };

	// Cartesian from Geodetic
	peri::LPA const haveLpa{ -1.832595715, 0.698131701, 1600.000 };
	peri::XYZ const wantXYZ{ peri::xyzForLpa(haveLpa) };

In terms of explicit data types:

	// include header for function declarations and definitions
    #include "peridetic.h" // indirectly includes implementation periDetail.h

	// Geodetic from Cartesian
	std::array<double, 3u> const gotLPA
		{ peri::lpaForXyz
			( std::array<double, 3u>{ -1266643.136, -4727176.539, 4079014.032 }
			)
		};
	// Cartesian from Geodetic
	std::array<double, 3u> const gotXYZ
		{ peri::xyzForLpa
			( std::array<double, 3u>{ -1.832595715, 0.698131701, 1600.000 }
			)
		};

### Definitive Example Code <a id=Definitive-Example-Code></a>

The following main program source code files provide compilable
use-case examples and also can be used as stand-alone coordinate conversion
utilities (e.g. for one off conversions, call from scripts, etc).

Demonstration/utility example programs include:

* [lpaForXyz.cpp](https://github.com/Stellacore/peridetic/blob/main/examples/lpaForXyz.cpp)
	-- Report equivalent Geodetic coordinate values for three command line
	Cartesian XYZ coordinate values expressed in meters.

* [xyzForLpaRadians.cpp](https://github.com/Stellacore/peridetic/blob/main/examples/xyzForLpaRadians.cpp)
	-- Report equivalent Cartesian coordinates for three command line
	geodetic coordinate values with longitude/parallel(latitude)
	expressed in radians and altitude expressed in meters.

* [xyzForLpaDegrees.cpp](https://github.com/Stellacore/peridetic/blob/main/examples/xyzForLpaDegrees.cpp)
	-- Report equivalent Cartesian coordinates for three command line
	geodetic coordinate values with longitude/parallel(latitude)
	expressed in (non-standard)*degrees* and altitude expressed in meters.



## Peridetic - Context and Alternatives <a id =Context-and-Alternatives></a>

Perietic is extremely focused software capability that exists in the overall
context of general geospatial technologies. There are many software 
applications and development resources available in this domain.

* ... OSGeo -- A particularly useful archive of complementary geospatial
technologies may be found at the OSGeo Foundation website:

	* https://www.osgeo.org/


For Geodetic/Cartesian coordinate conversions, there are a number
of existing alternative software options for performing Geodetic
transformations. Most of these include many additional capabilities
and involve installing large software packages. However, if you
need more extensive geodetic capabilities and/or additional features
(e.g. cartography, magnetism, etc.), these are well worth consideration.

* ... GeographicLib -- Full C++ environment with *many* additional
capabilities (geodesic paths, cartographic projections, magnetism, Geoid, etc).

	* https://sourceforge.net/projects/geographiclib/

	* By comparison with Peridetic, requires installing large software
	package/dependencies and data files.

* ... PROJ -- Full blown mapping package with *many* additional
capabilities especially in relation to cartographic projections and datum
accommodations.

	* https://www.osgeo.org/projects/proj/

	* By comparison with Peridetic, requires installing a large
	software package/dependency and database files.

After creating the Peridetic code, a previously existing similarly capable
and lightweight package was discovered that was missed during initial searches
prior to developing Peridetic.

* ... ecef-geodetic -- provides a collection of multiple ECEF-to-geodetic
coordinate conversion functions.

	* https://github.com/planet36/ecef-geodetic

	* By comparison with Peridetic, this offers a selection of
	multiple algorithms whereas Peridetic provides only a single (but
	well tested, externally verified, and fast) computation algorithm.



## Peridetic - Technical Detail <a id=Technical-Detail></a>

A mathematical description of the equations and formulae involved
is presented in the
[technical note: ./doc/perideticSummary.pdf](https://github.com/Stellacore/peridetic/tree/main/doc/perideticSummary.pdf)

### Terminology

A comment on terminology and notation:

* "XYZ" is used herein to denote classic Cartesian coordinates. In
	technical math speak - the three components of a true 3D
	vector as it is expressed with respect to a orthonormal dextral
	(right-handed) basis. For code purposes, all data values (input
	arguments and return values) are interpreted as _meters_.

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
Earth center are outside of the [design domain](#optimal-domain).
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
transformation functions (Ref peri::model namespace classes in
periDetail.h implementation header).

The transformations offer the WGS84 as a default ellipsoid so that
transformations are out-of-box compatible with most modern GNSS data.
Note that, in applications terms, the WGS84 ellipsoid is, for virtually
all practical purposes, simply a different way of specifying the same
shape as does the GRS80 ellipsoid (to within 0.1[mm] at the North pole).

### Peridetic - Software Considerations <a id=Software-Considerations></a>

At it's core, this Peridetic project comprises two source code header files:

* peridetic.h -- public interface (which has internal #include periDetail.h)

* periDetail.h -- implementation (inline functions) code

#### Software Environment

Software development points include:

Language

* Standard C++11 for header use (i.e. public headers), C++17 for test code.

Compilers

* GCC - primary (9.3), (10.1)

* Clang - verified (7.0.1)

Build

* Via [CMake](https://cmake.org/)

Documentation

* Via ['doxygen'](https://www.doxygen.nl/index.html)
Processors

* Precision requires availability of 8-byte (64-bit) or larger "double" type.

* Tested on

	* x86 - (primary: AMD Ryzen 2700)

Optimization:

* NOTE! To obtain [reasonable performance](#Transformation-Speed)
	it is important to compile with
	_*optimization enabled*_. (Implementation involves many, many short
	inline functions and variable assignment operations that entirely
	disappear when compiled with optimization, but will represent a
	very large number of (slow) function calls and (unnecessary) copy
	assignments if not compiled with optimization
	
		If using CMake, note that the default CMake setup
		configuration, by default, creates build instructions
		for debug mode. To obtain reasonable transformation
		performance, be sure to override this behavior. (E.g.
		with "-DCMAKE_BUILD_TYPE=Release" option).

Error Handling:

* No exceptions are involved or utilized within the code
	(none used, none thrown)

* Internal code supports quite_NaN propagation
	(if any input data are NaN, then result is returned with all
	components set to NaN values)

Thread safety:

* Functions should be(+) both re-entrant and thread safe

	* (+) expected to be true, but has not been verified explicitly.

Executable Code Overhead: <a id=Code-Overhead></a>

* Peridetic adds very little code/data size to executable programs.
	Ref. the code size example in [periUse sub directory](./periUse/README.md).


#### Use and Integration

This project can be used in your own code in two ways:

* Copy the "peridetic.h" public header file, and the "periDetail.h"
	implementation detail header file wherever you wish and
	incorporate this location into build include path, then include
	the "peridetic.h" file in your own source code, compile and go.

* -or- Incorporate this project into another development effort
	using the ["cmake" paradigm](https://cmake.org/documentation/).
	If using CMake in your own project, you can include something 
	like the folloing in relevant CMakeLists.txt file.

Example CMakeLists.txt file syntax:

	# Find (previously) installed perdietic library
	find_package(peridetic REQUIRED NO_MODULE)
	message(Found: ${peridetic_FOUND})
	message(Version: ${peridetic_VERSION})

	# dependency for myTarget
	target_link_libraries(
		${myTarget}
		PRIVATE
			peridetic::peridetic
		)

#### Licensing

This code may be used for any purpose (including commercial) provided
the copy right notice is retained as described in the
[MIT (aka X11) License)](https://directory.fsf.org/wiki/License:X11)

#### Data Preconditions

Function arguments have a generally valid interpretation if they contain
non-degenerate numeric values (i.e.  provided values are not NaN or
infinity values. Therefore, no input argument value testing is performed
within transformation functions.

However, if input data contain any NaN values, then those are propagated
through computations, such that _all components_ of the return array are
set to NaN (e.g. supports null value propagation associated with
"null object" design pattern).

##### Crazy Values

Internal checking of values is limited to conditions that may reasonably
be expected to occur in the course of processing meaningful data.

It is intentional that there is no code (nor overhead) to test for poorly
formed data values such as infinity or subnormal values. Any such
values are propagated through computations based on the properties and
characteristics of local compiler, libraries and hardware.

If an handling of questionable data values is important, then consider
wrapping the transforms functions inside a data qualification/validity
test along the lines of:

	double const wildValueComponent{ ... };
	// test for each component that might contain bad data
	if ((0. == wildValueComponent) || std::isnormal(wildValueComponent))
	{
		// Here, wildValueComponent will be a computationally valid component.
		// -- application can introduce further restrictions as appropriate
		// -- e.g. to restrict angle values to principal domain, etc., etc.
	}

NaN values are used internally to represent "null" (i.e. undefined or
uninitialized values). In general, these are not visible to consuming
code. However, if the calling code provides any NaN values as input
data then all-NaN value data instances will be returned. Consuming code,
if desired, can rely on this for use as part of a "null object pattern"
paradigm.

E.g. consumer code can verify that return data values are valid (which they
should always been unless input values are out of control) by checking for
quite_NaN. For invalid return values, if any component is bad, then *all*
components are set to NaN values, so that only one of the array elements
need be tested. E.g.

	peri::XYZ const ecefLoc{ peri::xyzForLpa(...) };
	if (! std::isnan(ecefLoc[0]))
	{
		// Life is good
	}

##### Data Interpretation and Integrity

Large magnitude input angle values, will be wrapped into principal angle
domains, and interpreted accordingly.

Therefore, beware of the easy potential mistake of providing an angle
in units of (wrong)degrees (or gradians) instead of correct units of
radians. For example (e.g. a numeric value of 90.0 is interpreted as an
angle of ~=1.7575743 not the potentially (incorrectly)expected quarter
turn of ~=1.5707963 radians).

Overall, consuming code is expected to be nominally responsible
for itself in terms of providing values that are meaningful in a geodetic
sense.

For best numeric stability provide angle values within principle domain
between +/-pi in order to capture full precision of angle values.

For example, negative altitude values with magnitude greater than
the ellipsoid's polar radius are not valid geodetic coordinates.
To avoid the overhead in testing for such a silly condition, Peridetic
assumes the consuming code is sufficiently responsible to avoid this. If
not, add your own data validation wrapper guards around the Peridetic
function calls.

Overall, as long as the calling code is reasonably responsible handling its
own data values, then everything should be fine. If you are uncomfortable
with this level of responsibility, you might consider utilizing a few
utility functions from the project test environment:

* Ref the [periLocal.h] (https://github.com/Stellacore/peridetic/blob/main/tests/periLocal.h)
	header file in "/tests" subdirectory. This is a header file used
	in development/testing programs. It includes various functions
	that may be generally useful in the context of gedetic data value
	interpretations (e.g. isValid(), infoString(), principalAngle(),
	etc). For usage information, build the project and point browser at
	the generated doxygen html documentation.



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
model is an ellipsoid. (Note: astrodetic lon/par values are referenced
to an undulating Geoid instead of an ellipsoid. Astrometric conversions
are outside the scope of these Peridetic transformations).

#### Reference Ellipsoids

Reference Ellipsoids - are used to define the origin and alignment
of the underlying coordinate frames for both XYZ and LPA data values.
A specific and/or custom reference ellipsoid may be provided to each
transformation. If none is provided then the WGS84 ellipsoid is used.

There are many subtleties in exactly what constitutes a
"best-fit" ellipsoid as it relates to the physical
["Figure of the Earth"](https://en.wikipedia.org/wiki/Figure_of_the_Earth).
However most commonly used ellipsoid models have dimensions that are
within a few 100[m] of each other and the most commonly used ellipsoids
differ from a pure spherical model by less than about +/-11[km] between
equatorial and polar semi-axes.

Peridetic provides two common Earth model specifications:

* GRS80 - ref peri::model::GRS80

* WGS84 - ref peri::model::WGS84

* Other standard as well as custom ellipsoids may be
	created easily via peri::Shape and peri::EarthModel data
	structures. Ref doxygen-comment documentation in periDetail.h file
	(viz the peri::model namespace documentation).

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

* Origin.A - is dynamically defined to lie *on the surface* of the ellipsoid
	at the location specified by the LP ellispoid angle coordinates.

* A(ltitude) - is the distance of a point from the Origin.A location.
	The "A" value is associated with a distance value that is
	differentially "stationary" with respect to small changes in
	LP coordinates. In practice, multiple solutions are resolved by
	selecting the "A" to be at point on the ellipsoid that is closest
	to the point of interest. The altitude is positive for points
	outside the ellipsoid (e.g. "upward"), and negative for points
	inside it (e.g. "downward").

For an explanation of Longitude and Latitude angles, ref: the NGS glossary
entry for
["coordinate, geocentric"](https://www.ngs.noaa.gov/PUBS_LIB/GEOID/Droman_papers/glossary/xml/C.xml).

Note that Peridetic uses the term "parallel" in place of latitude
(to facility unique naming notations) and uses the term "altitude"
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
ellipsoid. I.e. the LPA coordinate values are determined _BOTH_ by
the point location _AND_ by the ellipsoid acting in combination with
each other.

To resolve mathematical ambiguities expressing "the" LPA representation of
a point in space, the specific LPA representation is the one associated
with the point on the ellipsoid surface that is *closest* to the point
of interest.

Note that, even given a specific fixed ellipsoid, the LPA values still
are *_NOT_ unique*. As just one illustrative example, a point at zero
longitude, on the equator and on surface of ellipsoid has conventional LPA
value of (0,0,0). A second, mathematically valid solution, is associated
with the antipodal point on the other side of Earth that is also on the
equator but for which the altitude is an extreme negative value equal to
twice the equatorial radius, i.e. mathematical LPA value of (pi,0,-2b)
where 'b' is an equatorial radius of the ellipsoid.

Every point in space has at least this kind of dual LPA
representation. Point locations near and below Earth's surface are
typically associated with an additional pair of candidate LPA solutions.
Some locations are associated with an infinity of potential LPA
expressions. E.g. the origin (center of Earth) has an infinite number
of LPA representations (both poles as well as anywhere on equator).

In mathematical terms, there are generally four separate and discrete LPA
representations for a given XYZ location (not near the very center of the
ellipsoid). In some cases, two of these solutions are mathematically
"imaginary" and have no physical correspondence.

Although LPA values are fundamentally defined for locations
arbitrarily far away from the ellipsoid (i.e. in deep space), it may be
questionable if it is correct to be using Geodetic LPA values (instead of
e.g. astrodetic ones). E.g. At geostationary satellite distances (altitude
of ~3e5[km]), a one arc-minute deflection of vertical (difference of
between Geodetic and Astrodetic locations) corresponds to "horizontal"
position difference approaching 10[m].

Peridetic's Geodetic transformations are optimized for locations "on" and
"near" the surface of Earth (ref: [optimal domain](#optimal-domain)).
In this case, the singularities near to the center of Earth are irrelevant
(ref: [special cases](#special-cases)).

LPA values are not true coordinates (in the sense of "co-important
ordinate values). Rather the LPA expression comprises a triple of
hetrogenous values for every point location in *combination with* two
global, and often implicit, values expressing the shape of an ellipsoid of
revolution. The individual components generally can not be meaningfully
combined with each other nor with the components of other LPA locations
without involving elaborate transformation operations (like first
converting into corresponding XYZ expressions;-).

#### XYZ <a id=XYZ-Coordinates></a>

The abbreviation "XYZ" is used to denote coordinates in a Cartesian
Coordinate system that are also known commonly as the "ECEF" (Earth
Centered, Earth Fixed) coordinates and/or as "Rectangular Geocentric
Coordinates".

Ref: ["coordinate, Cartesian"](https://www.ngs.noaa.gov/PUBS_LIB/GEOID/Droman_papers/glossary/xml/C.xml)

##### Remarks on XYZ Coordinates

The XYZ (ECEF) Cartesian coordinate frame is often useful when working
with GNSS (Global Navigation Satellite System) observations and/or other
computations that are global in scope.

* Origin -  Is associated with the geometric center of the reference
	ellipsoid. For modern geodetic ellipsoids, the ellipsoid center is
	associated with the centroid of Earth's overall mass distribution.
	For older geodetic ellipsoids, the center may be defined such that
	(only a portion of) the ellipsoid has a "best" fit for its surface
	over a specific region or country.

* Z - is the axis orthogonal to the equator (The rotation plane of
	symmetry associated with the ellipsoid under consideration). On
	Earth, the positive "Z" axis points toward the North pole.

* X - is axis orthogonal to "Z" (in the equatorial plane) and directed
	toward the prime meridian. On Earth, when using the Greenwich prime
	meridian(+), the positive "X" points approximately toward the Gulf of
	Guinea.

* Y - is axis mutually orthogonal to Z and X in "right-hand" sense. On
	Earth (for Greenwich prime), the positive "Y" axis points to a
	location south of the Bay of Bengal.

(+) Historical Note: the historically defined "Greenwich Meridian"
identified with the "Airy Transit Circle" is slightly West of the
WGS84 0-longitude meridian. Cf.
[Why the Greenwich meridian moved](https://link.springer.com/content/pdf/10.1007/s00190-015-0844-y.pdf)

The XYZ is a classic orthonormal rectangular coordinate system. It is
associated with three basis vectors (aka "axes"). Each basis vector
has unit magnitude (is "normalized"). The three axes are mutually
perpendicular (are "orthogonal") and have a dextral (aka "right-handed")
chirality interpretation (in the cyclic order: "X", "Y", "Z").

Using XYZ coordinates, distances and angles can be computed directly
from the coordinate component values (e.g. via Pythagorean theorem,
the law of cosines, etc). The XYZ (vector) coordinates can be added,
subtracted, and otherwise manipulated *as* vectors (because they *are*
vectors).

### Optimal Domain on Altitude <a id=optimal-domain></a>

The good performance of computations and high quality of results produced
by the code in this project are associated with a particular optimal
domain. Results remain fairly useful outside of this domain, although
the quality and performance may be less than when operating inside the
optimum domain.

Practical applications are typically concerned with locations and motions
within the atmosphere, at sea, and on, or just below, the surface
of Earth.  A useful (conservative) definition of this region may be
specified arbitrarily as the volume of space within +/-100[km] of the
ellipsoidal shape being used to approximate the Earth surface.

There is nothing magic about this arbitrary +/-100[km] threshold. However
it is a useful limit at which to evaluate and express transformation
precision and accuracy since this region encompasses the overwhelming
majority of practical use-cases. Higher altitudes are conventionally
associated with "outer space" while lower altitudes are physically
inaccessible with current technology.

The [Transformation Precision](#Transformation-Precision) section describes
what to expect for transformation of locations outside this optimum
domain.



## Peridetic - Technical Deep Dive <a id=Technical-Deep-Dive></a>

### Transformation Algorithms

In general, geodetic transformations are non-linear operations for which
computation can be somewhat complex and involved. However, most practical
applications of Geodesy involve spatial locations that are fairly "near"
to Earth's surface (e.g. within +/-100[km]). Within this domain, the
transformations can exploit techniques that remain highly accurate and
precise while being simpler to implement and faster to run than various
"theoretically-elegant" analytical solutions.

There are a number of papers that espouse the use of "closed form"
solutions for the LPA from XYZ conversion. Certainly these approaches
could be utilized for implementation of Peridetic. However, those
solutions (generally based on solving a quartic equation) generally
involve a fair number of steps, sometimes also testing and switching
logic, and they require evaluation of n-th roots computations (i.e. cube
roots and square roots).

Numerically iterative approaches are sometimes considered inelegant.
However, consider that the n-th root computations required to evaluate
"closed form" geodetic solutions are not elementary algebraic
operations. Rather the root extraction computations are themselves
evaluated by iterative algorithms. These root computations may be
implemented within processor silicon and by very fast. However, for
some architectures, the iterative root algorithms are implemented in
processor firmware or in math library source code which may be much
slower than evaluating the simple algebraic instructions employed by
the direct numeric iteration in Peridetic.

Peridetic algorithms utilize a simple, [fast](#Runtime-Performance) and
[precise](#Transformation-Precision) direct numeric iteration approach
followed by the use of two trig function evaluations (for std::atan2())
used to return longitude/parallel solution as conventional angle values.

The mathematical formulae and algorithm detail is described in:

* PerideticMath.{lyx,pdf} -- TODO (cleanup/post)

### Transformation Precision <a id=Transformation-Precision></a>

For purposes here, "precision" is defined loosely as "self-consistency",
"repeatability", "computational significance". As such, it is an intrinsic
metric that describes the quality of the *computations* more so than
the quality associated with data values. (Quality of data *values* is
addressed in following section on accuracy).

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
		describe it briefly in an [email](mailto://peridetic@stellacore.com),
		and the project author will buy you a beverage.

 * -6300[km] <= Alt < -5800[km]: -- Reduced precision, < 100[um]

 * -5800[km] <= Alt < +11000[km]: -- Meets *design precision*, < 7.6[nm]

	* **-100[km] <= Alt <= 100[km]: -- is optimal design domain**

 * +11000[km] <= Alt < +405[Mm]: -- Reduced precision, < 0.2[um]

 * Beyond lunar distances, precision will continue to drop as
	altitude increases.  As an extreme example: at the altitude of
	Sagittarius A-star, the Milkyway's black hole, transform
	precision reduces to <100[km]. On the order of 1/4 of the way to
	Andromeda gallaxy, a 64-bit double completely loses all precision
	for expressing distances relative to the size of Earth.

The above precision reports are created by testing round trip
transformations for various point locations well distributed with respect
ellipsoid. For testing within the optimal domain, results for approximately
1M point locations have been evaluated.

For precision testing, the transformations are evaluated for
self-consistency.  Note that self-consistency does _not_ constitute a
proof of correctness, but only provides a measure of numeric/computational
noise involved.  Therefore, even if transformation computations are
precise, there is the completely independent question concerning how
accurate (correct) are the results.

The question of accuracy is addressed in the section
[Transform Accuracy](#Transformation-Accuracy) below.

### Transformation Accuracy <a id=Transformation-Accuracy></a>

To determine transformation accuracy it is necessary to compare results
with "absolute" or "known" values available from external sources.

This current implementation of Peridetic is evaluated using
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
([CORS](https://www.ngs.noaa.gov/CORS/)) network of station location
observations.  Although the network is concentrated within the U.S. it
includes a number of stations scattered around the globe.

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

* Cut-n-paste data from the NGS interactive map site into a header
	file (in the peridetic development and testing codebase). Ref
	[corsDataPairs.h](https://github.com/Stellacore/peridetic/blob/main/tests/corsDataPairs.h)
	header in /tests directory.

* Unit test code reads data from this header file and populates program
	variables with "expected" data values. E.g.:

	* expXYZ

	* expLPA -- after decoding d-m-s values into Radians.

* The expected values are transformed in each direction to obtain "got"
	values. E.g.

	* gotXYZ=xyzForLpa(expLPA)

	* gotLPA=lpaForXYZ(expXYZ).

* The "got" values are compared with their corresponding "expected" values
	and differences from zero are compared against tolerance values. The
	tolerance values used reflect the published data precision:

	* Angular tolerance of { 172. / 1024./1024./1024./1024. };
		This is <.16[nRad] (published angular values have resolution of
		approximately .049[nRad], but the published XYZ coordinates
		only have +/-1[mm] precision for comparison. A surface distance
		of 1[mm] corresponds with a larger angle ~.16[nRad] which is
		therefore used for testing.

	* Linear tolerance of { 1./1024. }; // CORS file XYZ published to [mm]

* All pairs of coordinates are required to pass this tolerance test (in
	both directions (lpaForXyz and xyzForLpa).

Overall, the Peridetic transforms are thought to be correct as described
above. The open source transparency of the algorithm and implementation
code along with opportunity for associated peer review should help
provides assurance on the quality. However...

...Please _note_ the _"AS IS"_ clause of the associated license.

### Transformation Quality and Testing

The quality metrics described above are checked by the unit and verification
test programs contained in the project "/tests" directory and may be
compiled and run using the CMake/CTest paradigm.

* Test programs in
	[/tests](https://github.com/Stellacore/peridetic/tree/main/tests)
	projects directory.

* Each test is an independent program.
	Diagnostic and error messages are sent to std::cout and std::cerr
	respectively.

* Tests may be built as described in
	[CMake Build](#CMake-Build)
	section

* The full test suite may run using CTest as described in
	[CMake Build](#CMake-Build)
	Individual tests may be run independently (e.g. from command line
	or desktop).


### Runtime Performance<a id=Runtime-Performance></a>

The Peridetic code implementation emphasizes clarity and portability.
Even though there are no special optimizations, the computational
performance is reasonably good. The following presents simple
(wall-clock-style) timing results obtained by transforming a reasonably
sized sample of points.

#### Example Timing Results

Configuration:

* Hardware: AMD Ryzen 7 2700  (using only *one* single processor core)

* Compiled with GCC 9.3

* Generated sample point distribution spanning:

	* All longitude values in the range [-pi:pi) radians

	* All parallel (latitude) values in the range [-pi/2:pi/2] radians

	* Altitudes in the range [-100,000:+100,000] meters

Each test includes transformation of the entire collection of samples. Each
transformation involves fetching a (pre-computed) data value from memory,
applying a specific transformation operation, and then putting the result
into a (pre-allocated) work space.

The test transformation operations include:

* "copy": Simply copy sample data from source into work space (no computation).
	This should provide an indication of data handling overhead which is
	present in *all* transformation operations

* "multiply": Multiply each component of the 3D data input by a constant
	value (the non-binary values of .1, .3, .7).

* "sqrt(abs)": Compute the square root of absolute value of each component.
	Note on this AMD processor, sqrt() operation is about as fast as 
	a multiply.

* "xyzForLpa": Full geodetic transformation: compute Cartesian "X,Y,Z"
	coordinates from Geodetic "Lon,Lat,Alt" values.

* "lpaForXyz": Full geodetic transformation: compute Geodetic "Lon,Lat,Alt"
	values from Cartesisn "X/Y/Z" coordinates.

Results:

Time values include all computation infrastructure and overhead as well
as specific computation times. The 'copy' transform is included to provide
a ballpark estimate of time associated with data handling overhead.

	# Number samples tested: 17508141

	# Absolute times per test
	# -- time values are 'wall-clock' elapsed [in sec]
	# -- absolute total and 'per-each' times

	0.385952462     0.000000022  : Reference evaluation - copy:
	0.476933213     0.000000027  : Reference evaluation - multiply:
	0.482271823     0.000000028  : Reference evaluation - sqrt(abs()):
	1.196965868     0.000000068  : Cartesian from Geodetic - xyzForLpa():
	2.642841151     0.000000151  : Geodetic from Cartesian - lpaForXyz():

	# Relative test times
	# -- test times with respect to each other [ratio]
	# -- column order matches row order

	1.00   0.81   0.80   0.32   0.15  : Reference evaluation - copy:
	1.24   1.00   0.99   0.40   0.18  : Reference evaluation - multiply:
	1.25   1.01   1.00   0.40   0.18  : Reference evaluation - sqrt(abs()):
	3.10   2.51   2.48   1.00   0.45  : Cartesian from Geodetic - xyzForLpa():
	6.85   5.54   5.48   2.21   1.00  : Geodetic from Cartesian - lpaForXyz():

Note the timing values fluctuate by a percent or two from run to run, but
this provides a general idea of what to expect (at least this class
of processor).

For this test, the xyzForLpa() computation is about two and a half times
more expensive than simple multiplication of all three coordinate values.

Computation of geodetic "Lon,Lat,Alt" values from Cartesian "X,Y,Z"
coordinates is approximately five and half times more expensive than simple
multiplication of all three coordinate values.

Note that these results apply to point locations within the optimum design
domain (locations within approximately +/-100[km] of Earth surface). For
points outside this range (e.g. geosynchronous orbits, etc), computations
require slightly (but only slightly) more time for the lpaForXyz()
conversion (approximately 50%, or possibly 100% longer).

