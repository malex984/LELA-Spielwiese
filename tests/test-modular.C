
/* tests/test-modular.C
 * Copyright (C) 2001, 2002 Bradford Hovinen,
 * Copyright (C) 2002 Dave Saunders
 *
 * Written by Bradford Hovinen <hovinen@cis.udel.edu>,
 *            Dave Saunders <saunders@cis.udel.edu>
 *
 * ------------------------------------
 * 2002-04-10 Bradford Hovinen <hovinen@cis.udel.edu>
 *
 * Rename from test-large-modular.C to test-modular.C; made other updates in
 * accordance with changes to MyModular interace.
 * ------------------------------------
 *
 * See COPYING for license information.
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
using namespace std;

int main (int argc, char **argv)
{
	static integer q1("18446744073709551557");
	static integer q2 = 2147483647U;
	static integer q3 = 65521U;
	static int q4 = 101;
	static int iterations = 1;
	static int trials = 100000;
	static int categories = 100;
	static int hist_level = 1;

	static Argument args[] = {
		{ 'K', "-K Q", "Operate over the \"ring\" GF(Q) [1] for integer modulus.", TYPE_INTEGER, &q1 },
		{ 'Q', "-Q Q", "Operate over the \"ring\" GF(Q) [1] for uint32 modulus.", TYPE_INTEGER, &q2 },
		{ 'q', "-q Q", "Operate over the \"ring\" GF(Q) [1] for uint16 modulus.", TYPE_INTEGER, &q3 },
		{ 'p', "-p P", "Operate over the \"ring\" GF(Q) [1] for uint8 modulus.", TYPE_INT, &q4 },
		{ 'i', "-i I", "Perform each test for I iterations.", TYPE_INT,     &iterations },
		{ 't', "-t T", "Number of trials for the random iterator test.", TYPE_INT, &trials },
		{ 'c', "-c C", "Number of categories for the random iterator test.", TYPE_INT, &categories },
		{ 'H', "-H H", "History level for random iterator test.", TYPE_INT, &hist_level },
		{ '\0' }
	};

	parseArguments (argc, argv, args);

	commentator.start("MyModular test suite", "MyModular ");
	bool pass = true;

//	MyModular<integer> F_integer (q1);
	MyModular<uint32> F_uint32 (q2.get_ui ());
	MyModular<uint16> F_uint16 (q3.get_ui ());
	MyModular<uint8> F_uint8 ((uint8) q4);
//	MyModular<float> F_float ((float) q4);

	// Make sure some more detailed messages get printed
	commentator.getMessageClass (INTERNAL_DESCRIPTION).setMaxDepth (6);
	commentator.getMessageClass (INTERNAL_DESCRIPTION).setMaxDetailLevel (Commentator::LEVEL_UNIMPORTANT);

//	if (!runRingTests (F_integer, "MyModular<integer>", iterations, false)) pass = false;
	if (!runRingTests (F_uint32,  "MyModular<uint32>",  iterations, false)) pass = false;
	if (!runRingTests (F_uint16,  "MyModular<uint16>",  iterations, false)) pass = false;
	if (!runRingTests (F_uint8,  "MyModular<uint8>",  iterations, false)) pass = false;
//	if (!runRingTests (F_float,  "MyModular<float>",  iterations, false)) pass = false;

	//if (!testRandomIterator (F_integer, "MyModular<integer>", trials, categories, hist_level)) pass = false;
	if (!testRandomIterator (F_uint32,  "MyModular<uint32>", trials, categories, hist_level)) pass = false;
	if (!testRandomIterator (F_uint16,  "MyModular<uint16>", trials, categories, hist_level)) pass = false;
	if (!testRandomIterator (F_uint8,  "MyModular<uint8>", trials, categories, hist_level)) pass = false;

	commentator.stop (MSG_STATUS (pass));
	return pass ? 0 : -1;
}

// Local Variables:
// mode: C++
// tab-width: 8
// indent-tabs-mode: t
// c-basic-offset: 8
// End:

// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s:syntax=cpp.doxygen:foldmethod=syntax
