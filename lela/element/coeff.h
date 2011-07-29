// -*- Mode: C++; tab-width: 2; indent-tabs-mode: nil; c-basic-offset: 2 -*-
/** @headerfile coeff.h lela/element/coeff.h
 * 
 * This file containg reference counting wrapper for Singular number.
 *
 * Copyright (C) 2011 Oleksandr Motsak
 *
 * Written by @author Oleksand Motsak <http://goo.gl/mcpzY>
 **/
/*****************************************************************************/

#ifndef __SINGULAR_element_coeff_H
#define __SINGULAR_element_coeff_H

#include <omalloc/omalloc.h>
#include <misc/auxiliary.h>

#include <coeffs/coeffs.h>
#include <coeffs/numbers.h>

#include <cassert>

/// @class ReferenceCountedElement
/// This structure is needed for Number: the smart pointer stores all
/// the data in here and completely manages it.
struct ReferenceCountedElement
{
  /// This omalloc Bin will contain all the instances of this structure
  static omBin lela_reference_counted_numbers_bin;

  typedef number SingularNumber;
  typedef coeffs SingularRing;

  /// Create _undefined_ number representation...
  ReferenceCountedElement(): m_counter(0), m_element(0), m_coeffs(NULL) {}

  /// Specify the number representation: Singular number and coeffs.
  ReferenceCountedElement(SingularNumber c, SingularRing R): m_counter(0), m_element(c), m_coeffs(R) {}  

  /// Use omalloc for allocating all such structures
  void* operator new ( std::size_t size )
  {
    //omTypeAlloc(void*, addr, size);
    return omAlloc0Bin(ReferenceCountedElement::lela_reference_counted_numbers_bin);
  }

  /// Use omalloc for deallocating all such structures
  void operator delete ( void* block )
  { //omfree( block );
    omFreeBin(block, ReferenceCountedElement::lela_reference_counted_numbers_bin);
  }

  /// Standard intrusive_ptr_add_ref for this structure
  static inline void intrusive_ptr_add_ref(ReferenceCountedElement* p)
  {
    assert( p!= NULL );
    assert( p->m_counter >= 0 );

    ++(p->m_counter); 
  }

  /// Standard intrusive_ptr_release for this structure, which handles
  /// the destruction of the underlying number as well
  static inline void intrusive_ptr_release(ReferenceCountedElement* p)
  {
    assert( p!= NULL );
    assert( p->m_counter > 0 );

    if( --(p->m_counter) == 0 )
    {
      assume( n_Test(p->m_element, p->m_coeffs) );
      n_Delete(& (p->m_element), p->m_coeffs); // Remove number,
      delete p; /// TODO: Use Omalloc to delete *m_pointer
    }
  }

  long           m_counter; ///< Reference Counter
  SingularNumber m_element; ///< Referenced Element
  SingularRing   m_coeffs;  ///< Element's Coeff. Domain, needed for delete

}; // Use omalloc for allocating/deallocating...?

omBin ReferenceCountedElement::lela_reference_counted_numbers_bin = omGetSpecBin(sizeof(ReferenceCountedElement));

/// @class Number This is a simple smart pointer guarding Singular numbers.
///
/// This class is needed Ring-Elements are expected to be able to
/// construct and destruct themselves.
/// 
/// Moreover its purpose is to count references and destroy the
/// underlying number when it's not references anymore.
///
/// The use of underlying Singular number before assiging or
/// initializing a Number is a bug and will trigger an assert failure
class Number
{
  public:
    typedef typename ReferenceCountedElement::SingularNumber SingularNumber;
    typedef typename ReferenceCountedElement::SingularRing   SingularRing;

    Number () : px (NULL) {}

    explicit Number( const int i, SingularRing R )
    {
      assert( R != 0 );
      
      SingularNumber c = n_Init(i, R);
      assume( n_Test(c, R) );

      px = new ReferenceCountedElement(c, R); /// TODO: Use Omalloc to create it!
      assume( px != 0 );

      ReferenceCountedElement::intrusive_ptr_add_ref( px );
    }

    explicit Number( SingularNumber c, SingularRing R, bool add_ref = true )
    {
      assert( R != 0 );
      
      assume( n_Test(c, R) );

      px = new ReferenceCountedElement(c, R); /// TODO: Use Omalloc to create it!
      assume( px != 0 );

      if( add_ref ) ReferenceCountedElement::intrusive_ptr_add_ref( px );
    }

    Number(Number const & rhs): px( rhs.px )
    {
      if( px != 0 ) ReferenceCountedElement::intrusive_ptr_add_ref( px );
    }

    ~Number()
    {
      if( px != 0 ) ReferenceCountedElement::intrusive_ptr_release( px );
      px = 0;      
    }

    Number & operator=(Number const & rhs)
    {
      Number(rhs).swap(*this);
      return *this;
    }

    /// Reset this instance to an undefined state
    inline void reset()
    {
      Number().swap( *this );
    }

    /// Reset this instance to a different state (c in R)
    inline void reset( SingularNumber c, SingularRing R )
    {
      Number( c, R ).swap( *this );
    }

    inline operator SingularNumber() const
    {
      return get_number();
    }

    inline bool test() const
    {
      if( px == 0 )
        return false;

      if (get_ring() == 0)
        return false;

      return n_Test( get_number(), get_ring() );
    }

    inline operator bool () const
    {
      return test();
    }

    // operator! is redundant, but some compilers need it
    inline bool operator! () const // never throws
    {
      return !test();
    }    

    inline bool belongs_to(SingularRing R) const 
    {
      if (px != 0)
      {
        assert( R != 0 );
        assert( R == get_ring() );
        return (R == get_ring());
      }

      return (R==0);
    }

    inline void swap(Number & rhs)
    {
      ReferenceCountedElement* tmp = px;
      px = rhs.px;
      rhs.px = tmp;
    }

    inline Number copy() const
    {
      const SingularRing R = get_ring();
      
      Number cpy( n_Copy(get_number(), R), R );

      assume( cpy.test() );
      assume( cpy.belongs_to(R) );

      assume( cpy == (*this) );
      assume( &cpy != this ); // Real Copy!
      assume( cpy.get_counter() == 1 );
      
      return cpy;
    }


    inline long get_counter() const // TODO: can be used in all the in-place ring arithmetic routines of CoeffDomain!
    {
      assert( px != 0 );
      const long c = px->m_counter;
      assert( c >= 0 );
      return c;
    }

    inline const SingularRing& get_ring() const
    {
      assert( px != 0 );      
      return px->m_coeffs;
    }

    /// Addition of Numbers
    Number operator+(const Number& right) const
    {
      const SingularRing R = get_ring();

      assume( right.test() );
      assume( test() );
      assume( right.belongs_to( R ) );
      
      return Number( n_Add( get_number(), right.get_number(), R ), R );
    }

    /// Subtraction of Numbers
    Number operator-(const Number& right) const
    {
      const SingularRing R = get_ring();
      
      assume( right.test() );
      assume( test() );
      assume( right.belongs_to( R ) );

      return Number( n_Sub( get_number(), right.get_number(), R ), R );
    }

    /// Multiplication of Numbers
    Number operator*(const Number& right) const
    {
      const SingularRing R = get_ring();
      
      assume( right.test() );
      assume( test() );
      assume( right.belongs_to( R ) );

      return Number( n_Mult( get_number(), right.get_number(), R ), R );
    }


    /// Division of Numbers
    Number operator/(const Number& right) const
    {
      const SingularRing R = get_ring();

      assume( right.test() );
      assume( test() );
      assume( right.belongs_to( R ) );

      return Number( n_Div( get_number(), right.get_number(), R ), R );
    }
    


    /// unary Additive Inverse (Negation) operator
    Number operator-() const
    {
      assume( test() );
      
      const SingularRing R = get_ring();

      // Note: n_Neg creates NO NEW number!
      return Number( n_Neg( n_Copy( get_number(), R ), R ), R );
    }

    /// Test whether this number is divisible by another (right)
    inline bool div_by(const Number& right) const
    {
      const SingularRing R = get_ring();

      assume( right.test() );
      assume( test() );
      assume( right.belongs_to( R ) );

      return n_DivBy(get_number(), right.get_number(), R);
    }
    
    inline bool operator==(const Number& right) const
    {
      const SingularRing R = get_ring();

      assume( right.test() );
      assume( test() );
      assume( right.belongs_to( R ) );

      // TODO: No additionaly actions if equal?
      // Proposal to test: remove that method's const modifier
      // and in equality case: assign *this = right
      // counter: sum of counters?
      // Problem: right is immutable... :(
      // Do we want to hack around that?
      return n_Equal(get_number(), right.get_number(), R);
    }


    Number& operator*=(const Number& right)
    {
      const SingularRing R = get_ring();

      assume( right.test() );
      assume( test() );
      assume( right.belongs_to( R ) );

      if( get_counter() == 1 )
      {
        SingularNumber& n = get_number();
        n_InpMult( n, right, R );
      }
      else
      {
        (*this) = (*this) * right;
      }
      return *this;
    }


    Number& negin()
    {
      const SingularRing R = get_ring();

      assume( test() );

      if( get_counter() == 1 )
      {
        SingularNumber& n = get_number();
        n = n_Neg( n, R ); // Note: n_Neg creates NO NEW number!

        assume( test() );
      }
      else
      {
        (*this) = -(*this);
      }
      return *this;     
    }
                                  
    SingularNumber get_number() const
    {
      assert( px != 0 );      
      assert( n_Test( px->m_element, get_ring() ) );
      return px->m_element;
    }
    
  private:
    SingularNumber& get_number()
    {
      assert( px != 0 );      
      assert( n_Test( px->m_element, get_ring() ) );
      return px->m_element;
    }

    ReferenceCountedElement* px; ///< contained pointer, with counter
};


void swap(Number & lhs, Number & rhs)
{
  lhs.swap(rhs);
}


#endif /* __SINGULAR_element_coeff_H */

