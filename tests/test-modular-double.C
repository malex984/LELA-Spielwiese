/* tests/test-modular-double.C
 * Brazenly stolen by bds from 
 * tests/test-modular-short.C
 * Brazenly stolen by Zhendong Wan (Copyright 2003) from 
 * tests/test-modular.C
 * Copyright 2001, 2002 Bradford Hovinen,
 * Copyright 2002 Dave Saunders
 *
 * Written by Bradford Hovinen <hovinen@gmail.com>,
 *            Dave Saunders <saunders@cis.udel.edu>
 *
 * ------------------------------------
 * 
 * This file is part of LELA, licensed under the GNU General Public
 * License version 3. See COPYING for more information.
 */

#include "lela/lela-config.h"

#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <queue>

#include "lela/ring/mymodular.h"

#include "test-common.h"
#include "test-ring.h"
#include "test-blas-level1.h"

using namespace LELA;

int main (int argc, char **argv)
{
	static int iterations = 1;

	static Argument args[] = {
		{ 'i', "-i I", "Perform each test for I iterations.", TYPE_INT,     &iterations },
		{ '\0' }
	};

	parseArguments (argc, argv, args);

	// Make sure some more detailed messages get printed
	commentator.getMessageClass (INTERNAL_DESCRIPTION).setMaxDepth (6);
	commentator.getMessageClass (INTERNAL_DESCRIPTION).setMaxDetailLevel (Commentator::LEVEL_UNIMPORTANT);
//	commentator.setReportStream(cout); // For debug...
	
	std::ostream& report = commentator.report();

	commentator.start("MyModular<double> ring test suite", "MyModular<double>");
	bool pass = true;


	// MyModular<double> F2 (2); 
	// report << "Ring F2" << std::endl;
	// if (!runRingTests (F2,  "MyModular<double>",  iterations,
	// false)) pass = false;
	
	MyModular<double> F3 (3); 
	report << "Ring F3" << std::endl;
	if (!runRingTests (F3,  "MyModular<double>",  iterations, false)) pass = false;
	
	MyModular<double> F5 (5); 
	report << "Ring F5" << std::endl;
	if (!runRingTests (F5,  "MyModular<double>",  iterations, false)) pass = false;
	
	MyModular<double> F7 (7); 
	report << "Ring F7" << std::endl;
	if (!runRingTests (F7,  "MyModular<double>",  iterations, false)) pass = false;

	MyModular<double> F11 (11); 
	report << "Ring F11" << std::endl;
	if (!runRingTests (F11,  "MyModular<double>",  iterations, false)) pass = false;
	
	MyModular<double> F (32749); 
	report << "Ring F" << std::endl;
	if (!runRingTests (F,  "MyModular<double>",  iterations, false)) pass = false;

	MyModular<double> G (65521); 
	report << "Ring G" << std::endl;
	if (!runRingTests (G,  "MyModular<double>",  iterations, false)) pass = false;


	//MyModular<double> H (1099511627689); 

	commentator.stop(MSG_STATUS (pass));
	return pass ? 0 : -1;
}
/* -*- mode: C++; tab-width: 8; indent-tabs-mode: t; c-basic-offset: 8 -*- */
// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s:syntax=cpp.doxygen:foldmethod=syntax
