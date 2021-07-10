

// include header for function declarations and definitions
#include "peridetic.h" // indirectly includes implementation periDetail.h

int
main
	()
{

	// Geodetic from Cartesian
	std::array<double, 3u> const gotLPA
		{ peri::lpaForXyz
			( std::array<double, 3u> { -1266643.136, -4727176.539, 4079014.032 }
			)
		};
	// Cartesian from Geodetic
	std::array<double, 3u> const gotXYZ
		{ peri::xyzForLpa
			( std::array<double, 3u> { -1.832595715, 0.698131701, 1600.000 }
			)
		};

	return 0;
}

