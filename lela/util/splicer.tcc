/* lela/util/splicer.tcc
 * Copyright 2011 Bradford Hovinen <hovinen@gmail.com>
 *
 * Written by Bradford Hovinen <hovinen@gmail.com>
 *
 * Class to chop matrices into pieces and splice them together
 *
 * ------------------------------------
 * 
 * This file is part of LELA, licensed under the GNU General Public
 * License version 3. See COPYING for more information.
 */

#ifndef __LELA_UTIL_SPLICER_TCC
#define __LELA_UTIL_SPLICER_TCC

#include "lela/util/splicer.h"
#include "lela/util/commentator.h"
#include "lela/util/debug.h"
#include "lela/blas/context.h"
#include "lela/blas/level1.h"
#include "lela/vector/bit-vector.h"
#include "lela/vector/bit-subvector.h"
#include "lela/vector/sparse-subvector.h"
#include "lela/matrix/submatrix.h"

namespace LELA
{

template <class Ring, class Vector1, class Vector2>
void Splicer::attach_block_specialised (const Ring &F, Vector1 &out, const Vector2 &in, size_t src_idx, size_t dest_idx, size_t size,
					VectorRepresentationTypes::Sparse, VectorRepresentationTypes::Sparse)
{
	typename Vector2::const_iterator i = std::lower_bound (in.begin (), in.end (), src_idx, VectorUtils::FindSparseEntryLB ());
	typename Vector2::const_iterator i_end = std::lower_bound (in.begin (), in.end (), src_idx + size, VectorUtils::FindSparseEntryLB ());

	for (; i != i_end; ++i)
		out.push_back (typename Vector1::value_type (i->first - src_idx + dest_idx, i->second));
}

template <class Ring, class Vector1, class Vector2>
void Splicer::attach_block_specialised (const Ring &F, Vector1 &out, const Vector2 &in, size_t src_idx, size_t dest_idx, size_t size,
					VectorRepresentationTypes::Sparse, VectorRepresentationTypes::Dense)
{
	typename Vector2::const_iterator i = in.begin () + src_idx;
	typename Vector2::const_iterator i_end = in.begin () + (src_idx + size);

	for (; i != i_end; ++i, ++dest_idx)
		if (!F.isZero (*i))
			out.push_back (typename Vector1::value_type (dest_idx, *i));
}

template <class Ring, class Vector1, class Vector2>
void Splicer::attach_block_specialised (const Ring &F, Vector1 &out, const Vector2 &in, size_t src_idx, size_t dest_idx, size_t size,
					VectorRepresentationTypes::Dense, VectorRepresentationTypes::Sparse)
{
	typename Vector2::const_iterator i = std::lower_bound (in.begin (), in.end (), src_idx, VectorUtils::FindSparseEntryLB ());
	typename Vector2::const_iterator i_end = std::lower_bound (in.begin (), in.end (), src_idx + size, VectorUtils::FindSparseEntryLB ());

	Context<Ring> ctx (F);
	typename Vector1::SubvectorType out_sub (out, dest_idx, dest_idx + size);
	BLAS1::scal (ctx, ctx.F.zero (), out_sub);

	for (; i != i_end; ++i)
		out[i->first - src_idx + dest_idx] = i->second;
}

template <class Ring, class Vector1, class Vector2>
void Splicer::attach_block_specialised (const Ring &F, Vector1 &out, const Vector2 &in, size_t src_idx, size_t dest_idx, size_t size,
					VectorRepresentationTypes::Dense01, VectorRepresentationTypes::Dense01)
{
	typedef BitSubvector<typename Vector1::iterator, typename Vector1::const_iterator> DestSubvector;
	typedef BitSubvector<typename Vector2::const_iterator, typename Vector2::const_iterator> SourceSubvector;

	SourceSubvector v1 (in.begin () + src_idx, in.begin () + (src_idx + size));
	DestSubvector v2 (out.begin () + dest_idx, out.begin () + (dest_idx + size));

	Context<Ring> ctx (F);

	BLAS1::copy (ctx, v1, v2);
}

template <class Ring, class Vector1, class Vector2>
void Splicer::attach_block_specialised (const Ring &F, Vector1 &out, const Vector2 &in, size_t src_idx, size_t dest_idx, size_t size,
					VectorRepresentationTypes::Sparse01, VectorRepresentationTypes::Sparse01)
{
	typename Vector2::const_iterator i = std::lower_bound (in.begin (), in.end (), src_idx);
	typename Vector2::const_iterator i_end = std::lower_bound (in.begin (), in.end (), src_idx + size);

	for (; i != i_end; ++i)
		out.push_back (typename Vector1::value_type (*i - src_idx + dest_idx));
}

template <class Vector>
void Splicer::append_word (Vector &v, size_t index, typename Vector::word_type word)
{
	typedef WordTraits<typename Vector::word_type> WT;

	typename Vector::Endianness::word_pair w;
	w.parts.low = word;
	w.parts.high = 0ULL;
	w.full = Vector::Endianness::shift_right (w.full, index & WT::pos_mask);

	if (v.empty () || v.back ().first != index >> WT::logof_size)
		v.push_back (typename Vector::value_type (index >> WT::logof_size, w.parts.low));
	else
		v.back ().second |= w.parts.low;

	if (w.parts.high != 0)
		v.push_back (typename Vector::value_type ((index >> WT::logof_size) + 1, w.parts.high));
}

template <class Ring, class Vector1, class Vector2>
void Splicer::attach_block_specialised (const Ring &F, Vector1 &out, const Vector2 &in, size_t src_idx, size_t dest_idx, size_t size,
					VectorRepresentationTypes::Hybrid01, VectorRepresentationTypes::Hybrid01)
{
	SparseSubvector<const Vector2, VectorRepresentationTypes::Hybrid01> v1 (in, src_idx, src_idx + size);
	typename SparseSubvector<const Vector2, VectorRepresentationTypes::Hybrid01>::const_iterator i_v1;

	for (i_v1 = v1.begin (); i_v1 != v1.end (); ++i_v1)
		append_word (out, dest_idx + (i_v1->first << WordTraits<typename Vector2::word_type>::logof_size), i_v1->second);
}

template <class Ring, class Vector1, class Vector2>
void Splicer::attach_block_specialised (const Ring &F, Vector1 &out, const Vector2 &in, size_t src_idx, size_t dest_idx, size_t size,
					VectorRepresentationTypes::Dense01, VectorRepresentationTypes::Sparse01)
{
	SparseSubvector<const Vector2, typename VectorTraits<Ring, Vector2>::RepresentationType> v1 (in, src_idx, src_idx + size);
	BitSubvector<typename Vector1::iterator, typename Vector1::const_iterator> v2 (out.begin () + dest_idx, out.begin () + (dest_idx + size));

	Context<Ring> ctx (F);

	BLAS1::copy (ctx, v1, v2);
}

template <class Ring, class Vector1, class Vector2>
void Splicer::attach_block_specialised (const Ring &F, Vector1 &out, const Vector2 &in, size_t src_idx, size_t dest_idx, size_t size,
					VectorRepresentationTypes::Sparse01, VectorRepresentationTypes::Dense01)
{
	typename Vector2::const_iterator i = in.begin () + src_idx;
	typename Vector2::const_iterator i_end = in.begin () + (src_idx + size);
	size_t idx = dest_idx;

	for (; i != i_end; ++i, ++idx)
		if (*i)
			out.push_back (idx);
}

template <class Ring, class Vector1, class Vector2>
void Splicer::attach_block_specialised (const Ring &F, Vector1 &out, const Vector2 &in, size_t src_idx, size_t dest_idx, size_t size,
					VectorRepresentationTypes::Sparse01, VectorRepresentationTypes::Hybrid01)
{
	typedef SparseSubvector<const Vector2, typename VectorTraits<Ring, Vector2>::RepresentationType> Subvector2;

	Subvector2 v1 (in, src_idx, src_idx + size);

	typename Subvector2::const_iterator i;

	for (i = v1.begin (); i != v1.end (); ++i) {
		typename Subvector2::word_type t;
		size_t idx = (i->first << WordTraits<typename Subvector2::word_type>::logof_size) + dest_idx;

		for (t = Subvector2::Endianness::e_0; t != 0; t = Subvector2::Endianness::shift_right (t, 1), ++idx)
			if (i->second & t)
				out.push_back (idx);
	}
}

template <class Ring, class Vector1, class Vector2>
void Splicer::attach_block_specialised (const Ring &F, Vector1 &out, const Vector2 &in, size_t src_idx, size_t dest_idx, size_t size,
					VectorRepresentationTypes::Hybrid01, VectorRepresentationTypes::Dense01)
{
	typedef BitVector<typename Vector1::Endianness> TmpVector;
	typedef BitSubvector<typename TmpVector::iterator, typename TmpVector::const_iterator> TmpSubvector;
	typedef BitSubvector<typename Vector2::const_iterator, typename Vector2::const_iterator> ConstTmpSubvector;

	size_t dest_idx_offset = dest_idx & WordTraits<typename Vector2::word_type>::pos_mask;

	TmpVector tmp (round_up (size + dest_idx_offset, WordTraits<typename Vector2::word_type>::bits));
	ConstTmpSubvector v1 (in.begin () + src_idx, in.begin () + (src_idx + size));
	TmpSubvector v2 (tmp.begin () + dest_idx_offset, tmp.begin () + (dest_idx_offset + size));

	Context<Ring> ctx (F);

	BLAS1::copy (ctx, v1, v2);

	typename TmpVector::word_iterator i = tmp.word_begin ();
	size_t idx = dest_idx >> WordTraits<typename Vector1::word_type>::logof_size;

	if (i == tmp.word_end ()) {
		if (!out.empty () && out.back ().first == idx)
			out.back ().second |= tmp.back_word ();
		else if (tmp.back_word ())
			out.push_back (typename Vector1::value_type (idx, tmp.back_word ()));
	} else {
		if (!out.empty () && out.back ().first == idx) {
			out.back ().second |= *i;
			++i;
			++idx;
		}

		for (; i != tmp.word_end (); ++i, ++idx)
			if (*i)
				out.push_back (typename Vector1::value_type (idx, *i));

		if (tmp.back_word ())
			out.push_back (typename Vector1::value_type (idx, tmp.back_word ()));
	}
}

template <class Ring, class Vector1, class Vector2>
void Splicer::attach_block_specialised (const Ring &F, Vector1 &out, const Vector2 &in, size_t src_idx, size_t dest_idx, size_t size,
					VectorRepresentationTypes::Hybrid01, VectorRepresentationTypes::Sparse01)
{
	SparseSubvector<const Vector2, typename VectorTraits<Ring, Vector2>::RepresentationType> v1 (in, src_idx, src_idx + size);

	typename SparseSubvector<const Vector2, typename VectorTraits<Ring, Vector2>::RepresentationType>::const_iterator i_v1;

	for (i_v1 = v1.begin (); i_v1 != v1.end (); ++i_v1) {
		typename Vector1::word_type m = Vector1::Endianness::e_j (*i_v1 & WordTraits<typename Vector1::word_type>::pos_mask);

		if (*i_v1 >> WordTraits<typename Vector1::word_type>::logof_size == out.back ().first)
			out.back ().second |= m;
		else
			out.push_back (typename Vector1::value_type (*i_v1 >> WordTraits<typename Vector1::word_type>::logof_size, m));
	}
}

template <class Ring, class Vector1, class Vector2>
void Splicer::attach_block_specialised (const Ring &F, Vector1 &out, const Vector2 &in, size_t src_idx, size_t dest_idx, size_t size,
					VectorRepresentationTypes::Dense01, VectorRepresentationTypes::Hybrid01)
{
	SparseSubvector<const Vector2, typename VectorTraits<Ring, Vector2>::RepresentationType> v1 (in, src_idx, src_idx + size);
	BitSubvector<typename Vector1::iterator, typename Vector1::const_iterator> v2 (out.begin () + dest_idx, out.begin () + (dest_idx + size));

	Context<Ring> ctx (F);

	BLAS1::copy (ctx, v1, v2);
}

template <class Ring, class Matrix1, class Matrix2>
void Splicer::copyBlockSpecialised (const Ring &R, const Matrix1 &source, Matrix2 &dest, const Block &horiz_block, const Block &vert_block,
				    MatrixIteratorTypes::Row, MatrixIteratorTypes::Row)
{
	typename Matrix2::RowIterator i_A = dest.rowBegin () + horiz_block.destIndex ();
	typename Matrix2::RowIterator i_A_end = dest.rowBegin () + (horiz_block.destIndex () + horiz_block.size ());
	typename Matrix1::ConstRowIterator i_S = source.rowBegin () + horiz_block.sourceIndex ();

	for (; i_A != i_A_end; ++i_A, ++i_S)
		attach_block (R, *i_A, *i_S, vert_block.sourceIndex (), vert_block.destIndex (), vert_block.size ());
}

template <class Ring, class Matrix1, class Matrix2>
void Splicer::copyBlockSpecialised (const Ring &R, const Matrix1 &source, Matrix2 &dest, const Block &horiz_block, const Block &vert_block,
				    MatrixIteratorTypes::Col, MatrixIteratorTypes::Col)
{
	typename Matrix2::ColIterator i_A = dest.colBegin () + vert_block.destIndex ();
	typename Matrix2::ColIterator i_A_end = dest.colBegin () + (vert_block.destIndex () + vert_block.size ());
	typename Matrix1::ConstColIterator i_S = source.colBegin () + vert_block.sourceIndex ();

	for (; i_A != i_A_end; ++i_A, ++i_S)
		attach_block (R, *i_A, *i_S, horiz_block.sourceIndex (), horiz_block.destIndex (), horiz_block.size ());
}

template <class Ring, class Matrix>
void Splicer::copyIdentitySpecialised (const Ring &R, Matrix &dest, const Block &horiz_block, const Block &vert_block,
				       MatrixIteratorTypes::Row)
{
	typename Matrix::RowIterator i_A = dest.rowBegin () + horiz_block.destIndex ();
	typename Matrix::RowIterator i_A_end = dest.rowBegin () + (horiz_block.destIndex () + horiz_block.size ());
	size_t idx = horiz_block.sourceIndex ();

	for (; i_A != i_A_end; ++i_A, ++idx)
		if (vert_block.isSourceIndexInBlock (idx))
			Splicer::attach_e_i (R, *i_A, vert_block.sourceToDestIndex (idx));
}

template <class Ring, class Matrix>
void Splicer::copyIdentitySpecialised (const Ring &R, Matrix &dest, const Block &horiz_block, const Block &vert_block,
				       MatrixIteratorTypes::Col)
{
	typename Matrix::ColIterator i_A = dest.colBegin () + vert_block.destIndex ();
	typename Matrix::ColIterator i_A_end = dest.colBegin () + (vert_block.destIndex () + vert_block.size ());
	size_t idx = vert_block.sourceIndex ();

	for (; i_A != i_A_end; ++i_A, ++idx)
		if (vert_block.isSourceIndexInBlock (idx))
			Splicer::attach_e_i (R, *i_A, horiz_block.sourceToDestIndex (idx));
}

template <class Grid>
void Splicer::splice (Grid grid) const
{
	lela_check (check ());

	commentator.start ("Splicing matrices", __FUNCTION__);

	typename std::vector<Block>::const_iterator horiz_block, vert_block;

	for (horiz_block = _horiz_blocks.begin (); horiz_block != _horiz_blocks.end (); ++horiz_block)
		for (vert_block = _vert_blocks.begin (); vert_block != _vert_blocks.end (); ++vert_block)
			grid (*horiz_block, *vert_block);

	commentator.stop (MSG_DONE);
}

} // namespace LELA

#endif // __LELA_UTIL_SPLICER_TCC

// Local Variables:
// mode: C++
// tab-width: 8
// indent-tabs-mode: t
// c-basic-offset: 8
// End:

// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s:syntax=cpp.doxygen:foldmethod=syntax
