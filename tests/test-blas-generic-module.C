/* tests/test-matrix-domain.C
 * Copyright 2001, 2002 Bradford Hovinen
 *
 * Written by Bradford Hovinen <bghovinen@math.uwaterloo.ca>
 *
 * ---------------------------------------------------------
 *
 * See COPYING for license information.
 *
 * Test suite for MatrixDomain
 */

#include "linbox/util/commentator.h"
#include "linbox/blas/context.h"
#include "linbox/ring/modular.h"
#include "linbox/matrix/dense.h"
#include "linbox/matrix/sparse.h"
#include "linbox/vector/stream.h"
#include "linbox/matrix/transpose.h"

#include "test-common.h"
#include "test-blas-level1.h"
#include "test-blas-level2.h"
#include "test-blas-level3.h"

using namespace LinBox;

int main (int argc, char **argv)
{
	bool pass = true;

	static long l = 50;
	static long m = 50;
	static long n = 40;
	static long p = 30;
	static long k = 10;
	static integer q = 2147483647U;
	static int iterations = 1;

	static Argument args[] = {
		{ 'l', "-l L", "Set row-dimension of matrix A to L.", TYPE_INT, &l },
		{ 'm', "-m M", "Set row-dimension of matrix B and column-dimension of A to M.", TYPE_INT, &m },
		{ 'n', "-n N", "Set row-dimension of matrix C and column-dimension of B to N.", TYPE_INT, &n },
		{ 'p', "-p P", "Set column-dimension of matrix C to P.", TYPE_INT, &p },
		{ 'k', "-k K", "K nonzero elements per row/column in sparse matrices.", TYPE_INT, &k },
		{ 'q', "-q Q", "Operate over the \"field\" GF(Q) [1] for uint32 modulus.", TYPE_INTEGER, &q },
		{ 'i', "-i I", "Perform each test for I iterations.", TYPE_INT, &iterations },
		{ '\0' }
	};

	parseArguments (argc, argv, args);

	typedef Modular<uint32> Ring;
	typedef Ring::Element Element;

	Ring F (q);
	Context<Ring, GenericModule<Ring> > ctx (F);

	commentator.setBriefReportParameters (Commentator::OUTPUT_CONSOLE, false, false, false);
	commentator.getMessageClass (INTERNAL_DESCRIPTION).setMaxDepth (5);
	commentator.getMessageClass (INTERNAL_DESCRIPTION).setMaxDetailLevel (Commentator::LEVEL_UNIMPORTANT);
	commentator.getMessageClass (TIMING_MEASURE).setMaxDepth (3);

	commentator.start ("BLAS GenericModule test suite", "BLASGenericModule");

	if (!testBLAS1 (ctx, "Modular <uint32>", l, iterations)) pass = false;

	RandomDenseStream<Ring, Vector<Ring>::Dense> stream_v1 (F, l, 1);
	RandomDenseStream<Ring, Vector<Ring>::Dense> stream_v2 (F, m, 1);
	RandomDenseStream<Ring, Vector<Ring>::Dense> stream_v3 (F, n, 1);
	RandomDenseStream<Ring, Vector<Ring>::Dense> stream_v4 (F, p, 1);

	Vector<Ring>::Dense v1 (l), v2 (m), v3 (n), v4 (p);
	stream_v1 >> v1;
	stream_v2 >> v2;
	stream_v3 >> v3;
	stream_v4 >> v4;

	RandomDenseStream<Ring, DenseMatrix<Element>::Row> stream11 (F, m, l);
	RandomDenseStream<Ring, DenseMatrix<Element>::Row> stream12 (F, n, m);
	RandomDenseStream<Ring, DenseMatrix<Element>::Row> stream13 (F, p, n);

	DenseMatrix<Element> M1 (stream11);
	DenseMatrix<Element> M2 (stream12);
	DenseMatrix<Element> M3 (stream13);

	if (!testBLAS2 (ctx, "dense", M1, M2, v1, v2, iterations,
			DenseMatrix<Element>::IteratorType ()))
		pass = false;
	if (!testBLAS3 (ctx, "dense", M1, M2, M3, iterations,
			DenseMatrix<Element>::IteratorType ()))
		pass = false;

	RandomSparseStream<Ring, SparseMatrix<Element>::Row> stream21 (F, (double) k / (double) m, m, l);
	RandomSparseStream<Ring, SparseMatrix<Element>::Row> stream22 (F, (double) k / (double) n, n, m);
	RandomSparseStream<Ring, SparseMatrix<Element>::Row> stream23 (F, (double) k / (double) p, p, n);

	SparseMatrix<Element> M4 (stream21);
	SparseMatrix<Element> M5 (stream22);
	SparseMatrix<Element> M6 (stream23);

	if (!testBLAS2 (ctx, "sparse row-wise", M4, M5, v1, v2, iterations,
			SparseMatrix<Element>::IteratorType ()))
		pass = false;
	if (!testBLAS3 (ctx, "sparse row-wise", M4, M5, M6, iterations,
			SparseMatrix<Element>::IteratorType ()))
		pass = false;

	TransposeMatrix<SparseMatrix<Element> > M7 (M6);
	TransposeMatrix<SparseMatrix<Element> > M8 (M5);
	TransposeMatrix<SparseMatrix<Element> > M9 (M4);

	if (!testBLAS2 (ctx, "sparse column-wise", M7, M8, v4, v3, iterations,
			TransposeMatrix<SparseMatrix<Element> >::IteratorType ()))
		pass = false;
	if (!testBLAS3 (ctx, "sparse column-wise", M7, M8, M9, iterations,
			TransposeMatrix<SparseMatrix<Element> >::IteratorType ()))
		pass = false;

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