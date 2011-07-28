// -*- Mode: C++; tab-width: 2; indent-tabs-mode: nil; c-basic-offset: 2 -*-
/** @headerfile coeffs.h lela/ring/coeffs.h
 * 
 * This file containg the LELA ring interface wrapper for Number
 *
 * ABSTRACT: This wrapper are needed in order to use the
 * Singular numbers inside the LELA.
 *
 * Copyright (C) 2011 Oleksandr Motsak
 *
 * Written by @author Oleksand Motsak <http://goo.gl/mcpzY>
 **/
/*****************************************************************************/

#ifndef __SINGULAR_ring_coeff_H
#define __SINGULAR_ring_coeff_H

#include <lela/lela-config.h>

#ifdef __LELA_HAVE_LIBPOLYS

#include <lela/ring/interface.h>
#include <lela/util/error.h>

#include <lela/element/coeff.h>

#include <omalloc/omalloc.h>
#include <misc/auxiliary.h>

#include <coeffs/coeffs.h>
#include <coeffs/numbers.h>

#include <reporter/reporter.h>
#include <resources/feResource.h>

#include <cassert>


/** Singular coefficient domain
 *
 * This class defines the ring-interface for general Singular coefficients.
 */
class CoeffDomain
{
  private:
    CoeffDomain& operator=(const CoeffDomain&);

  protected:
    typedef typename Number::SingularNumber BaseElement;
    typedef typename Number::SingularRing  BaseCoeffs;

    void CleanMe()
    {
      if (_coeffs != NULL)
      {
        _zero.reset();
        _one.reset();
        _minus_one.reset();

        if (_clean_coeffs)
          nKillChar(_coeffs);
        
        _type = n_unknown;
        _param = NULL;        
      }
    }

    void InitMe(BaseCoeffs R, bool cleanup)
    {
      CleanMe();
      
      assume (R != 0);

      if( R == 0 )
        throw LELA::LELAError ("Could not create the needed Singular coeffs!");        

      assert (R != 0);
      
      _coeffs = R;
      _clean_coeffs = cleanup;

      init( _zero, (int)0 );
      init( _one,  (int)1 );
      init( _minus_one, (int)-1 );
    }

    
    void InitMe(n_coeffType type, void* p)
    {
      InitMe( nInitChar( type, p ), true);
        
      _type = type;
      _param = p;
    }
    
  public:
    
    CoeffDomain(): _coeffs( NULL ), _clean_coeffs(false) {}

    
    CoeffDomain(n_coeffType type, void* p = NULL): _coeffs( NULL ), _clean_coeffs(false)
    {
      InitMe(type, p);
    }

    CoeffDomain(const CoeffDomain& D): _coeffs( NULL ), _clean_coeffs(false)
    {
      InitMe(D._type, D._param);
    }
    

    ~CoeffDomain()
    {
      CleanMe();
    }

//    CoeffDomain(coeffs &r): _coeffs(r), _clean_coeffs(false)    { InitMe();    }
    
    /** @name Singular coeffs Interface as a LELA Ring.
     * These methods are required of all \ref{LELA} rings.
     */
    //@{

    /// the type in which ring elements are represented.
    typedef Number Element;
    
//TODO!//    /// An object of this type is a generator of random ring elements.
//TODO!//    typedef RandIterInterface RandIter;

    /// @name Object Management
    //@{    
    Element &init (Element &x) const
    {
      x = zero();
      
      assume( x );
      assume( x.belongs_to( _coeffs ) );

      return x;
    }
    
    inline Element &init(Element &x, BaseElement y) const
    {
      assume( n_Test(y, _coeffs ) );

      x = Number(y, _coeffs);

      assume( x );
      assume( x.belongs_to( _coeffs ) );
      assume( x.get_counter() == 1 );

      return x;
      
    }

    /** \brief Initialization of ring element from an integer.
     *
     * x becomes the image of n under the natural map from the
     * integers to the ring. The element x need not have been
     * previously initialised.
     *
     * @return reference to x.
     * @param x output ring element.
     * @param n input integer.
     */
    Element &init (Element &x, const int i) const
    {
      x = Number(i, _coeffs);

      assume( x );
      assume( x.belongs_to( _coeffs ) );
      assume( x.get_counter() == 1 );

      return x;
    }

    Element &init (Element &x, const unsigned int &y) const
    { 
      return init(x, (int)(y));
    }

    // TODO: would not it fail on machines with int == long?
    Element &init (Element &x, const unsigned long &y) const
    { 
      return init(x, (int)(y));
    }

    Element &init (Element &x, const LELA::integer &y) const
    {
      const int n = y.get_ui();
      return init(x, n );
    }
    
    Element &init (Element &x, const double &y) const
    { 
      return init(x, (int)(y+.5));
    }
    
    Element &init (Element &x, const float &y) const
    { 
      return init(x, (int)(y+.5));
    }

    /// Version of init which takes a Property rather than an element.
    ///
    /// It should do exactly what the version above taking an element does.
    template <class Iterator, class Accessor>
        Element &init (LELA::Property<Iterator, Accessor> x, const LELA::integer &n) const
    {
      init(x.ref (), n);
      
      assume( x.ref () );
      assume( x.ref ().belongs_to( _coeffs ) );

      return x.ref ();
      
    }

    template <class Iterator, class Accessor>
        Element &init (LELA::Property<Iterator, Accessor> x) const
    {
      x.ref() = zero();
      
      assume( x.ref () );
      assume( x.ref ().belongs_to( _coeffs ) );

      return x.ref ();
    }

  /** \brief  Copy one ring element to another.
	 *
	 * This function makes a deep copy of the element y into the
	 * element x (as opposed to the assignment-operator, which
	 * makes only a shallow copy). The element x need not have
	 * been previously initialised.
	 *
	 * @return reference to x
	 * @param  x destination ring element.
	 * @param  y source ring element.
	 */
    Element &copy (Element &x, const Element &y) const
    {      
      assume( y );
      assume( y.belongs_to( _coeffs ) );
      
      x = y.copy();
      
      assume( x );
      assume( x.belongs_to( _coeffs ) );
      assume( x.get_counter() == 1 );

      return x;
      
    }

    /// Version of copy which takes a Property rather than an element.
    ///
    /// It should do exactly what the version above taking an element does.
    template <class Iterator, class Accessor>
        Element &copy (LELA::Property<Iterator, Accessor> x, const Element &y) const
    {
      assume( y );
      assume( y.belongs_to( _coeffs ) );
      
      return copy (x.ref (), y);
    }

    /** \brief Cardinality.
     *
     * Return c, integer representing cardinality of the ring.  c
     * becomes a positive integer for all rings with finite
     * cardinality, and 0 to signify a ring of infinite
     * cardinality.
     */
    virtual int &cardinality (int &c) const
    {
      return characteristic(c); // TODO: Overwrite it in speciall cases...?
    }
    
    LELA::integer &cardinality (LELA::integer &c) const
    {
      int _c; cardinality(_c);
      return c = LELA::integer(_c);
    }
    

    /** \brief Characteristic.
     *
     * Return c, integer representing characteristic of the ring
     * (the least positive n such that the sum of n copies of x is
     * 0 for all ring elements x).  c becomes a positive integer
     * for all rings with finite characteristic, and 0 to signify
     * a ring of infinite characteristic.
     */
    LELA::integer &characteristic (LELA::integer &c) const
    {
      int _c; characteristic(_c);
      return c = LELA::integer(_c);
    }
    
    int &characteristic (int &c) const
    {
      c = n_GetChar(_coeffs);
      return c;      
    }
    


    //@} Object Management

    /** @name Arithmetic Operations 
     * x <- y op z; x <- op y
     * These operations require all elements, including x, to be initialized
     * before the operation is called.  Uninitialized ring elements will
     * give undefined results.
     */
    //@{

    /** \brief Equality of two elements.
     *
     * This function assumes both ring elements have already been 
     * constructed and initialized.
     *
     * @return boolean true if equal, false if not.
     * @param  x ring element
     * @param  y ring element
     */
    bool areEqual (const Element &x, const Element &y) const
    {
      assume( x );
      assume( x.belongs_to( _coeffs ) );
      assume( y );
      assume( y.belongs_to( _coeffs ) );
      
      return n_Equal(x, y, _coeffs);
    }

    /** \brief Addition, x <-- y + z.
     *
     * This function assumes all the ring elements have already been 
     * constructed and initialized.
     *
     * @return reference to x.
     */
    Element &add (Element &x, const Element &y, const Element &z) const
    {
      assume( y );
      assume( y.belongs_to( _coeffs ) );
      assume( z );
      assume( z.belongs_to( _coeffs ) );
      
      x = y + z;
      
      assume( x );
      assume( x.belongs_to( _coeffs ) );

      return x;
    }

    /** \brief Subtraction, x <-- y - z.
     *
     * This function assumes all the ring elements have already been 
     * constructed and initialized.
     *
     * @return reference to x.
     */
    Element &sub (Element &x, const Element &y, const Element &z) const
    {
      assume( y );
      assume( y.belongs_to( _coeffs ) );
      assume( z );
      assume( z.belongs_to( _coeffs ) );

      x = y - z;
      
      assume( x );
      assume( x.belongs_to( _coeffs ) );

      return x;
    }

    /** \brief Multiplication, x <-- y * z.
     *
     * This function assumes all the ring elements have already been 
     * constructed and initialized.
     *
     * @return reference to x.
     */
    Element &mul (Element &x, const Element &y, const Element &z) const
    {
      assume( y );
      assume( y.belongs_to( _coeffs ) );
      assume( z );
      assume( z.belongs_to( _coeffs ) );

      x = y * z;
          
      assume( x );
      assume( x.belongs_to( _coeffs ) );

      return x;
    }

    /** Division, x <-- y / z.
     *
     * This function attempts to compute the quotient of y and z
     * and, if it exists, stores the result in x. If the quotient
     * does not exist in the ring, the value of x is left
     * unchangd.
     *
     * This function assumes all the ring elements have already been 
     * constructed and initialized.
     *
     * @return true if the quotient exists in the ring, false if not
     */
    bool div (Element &x, const Element &y, const Element &z) const
    {
      assume( y );
      assume( y.belongs_to( _coeffs ) );
      assume( z );
      assume( z.belongs_to( _coeffs ) );

      if( !n_DivBy(y, z, _coeffs) )
        return false;

      init(x, n_Div(y, z, _coeffs));

      assume( x );
      assume( x.belongs_to( _coeffs ) );
      
      return true;
    }

    /** \brief Additive Inverse (Negation), x <-- - y.
     *
     * This function assumes both ring elements have already been 
     * constructed and initialized.
     *
     * @return reference to x.
     */
    Element &neg (Element &x, const Element &y) const
    {
      assume( y );
      assume( y.belongs_to( _coeffs ) );

      init(x, n_Neg( n_Copy(y, _coeffs), _coeffs) ); // Note: n_Neg creates NO NEW number!
      
      assume( x );
      assume( x.belongs_to( _coeffs ) );

      return x;

    }

    /** \brief Multiplicative Inverse, x <-- 1 / y.
     *
     * This computes the multiplicative inverse if possible and if
     * it exists, stores it in x. If not, then x is left
     * unchanged.
     *
     * This function assumes both ring elements have already been 
     * constructed and initialized.
     *
     * @return true if the inverse exists in the ring, false if not
     */
    bool inv (Element &x, const Element &y) const
    {
      assume( y );
      assume( y.belongs_to( _coeffs ) );

      if( !n_DivBy(one(), y, _coeffs) )
        return false;

      init(x, n_Invers(y, _coeffs));

      assume( x );
      assume( x.belongs_to( _coeffs ) );
      
      return true;      
      
    }

    /** \brief Ring element AXPY, r  <-- a * x + y.
     *
     * This function assumes all ring elements have already been 
     * constructed and initialized.
     *
     * @return reference to r.
     */
    Element &axpy (Element &r, const Element &a, const Element &x, const Element &y) const
    {
      assume( a );
      assume( a.belongs_to( _coeffs ) );
      assume( x );
      assume( x.belongs_to( _coeffs ) );
      assume( y );
      assume( y.belongs_to( _coeffs ) );
      
      return addin(mul(r, a, x), y);
      
    }

    //@} Arithmetic Operations

    /** @name Predicates
     */
    //@{
    /** Equality with the zero-element
     *
     * Test if ring element is equal to zero.
     * This function assumes the ring element has already been 
     * constructed and initialized.
     *
     * @return boolean true if equals zero, false if not.
     * @param  x ring element.
     */
    bool isZero (const Element &x) const
    {
      assume( x );
      assume( x.belongs_to( _coeffs ) );
      
      return n_IsZero(x, _coeffs);
    }

    /** Equality with the one-element
     *
     * Test if ring element is equal to one.
     * This function assumes the ring element has already been 
     * constructed and initialized.
     *
     * @return boolean true if equals one, false if not.
     * @param  x ring element.
     */
    bool isOne (const Element &x) const
    {
      assume( x );
      assume( x.belongs_to( _coeffs ) );

      return n_IsOne(x, _coeffs);
    }
    //@}

    /** @name Inplace Arithmetic Operations 
     * x <- x op y; x <- op x
     *
     * These operations require all elements, including x, to be initialized
     * before the operation is called.  Uninitialized ring elements will
     * give undefined results.
     */
    //@{

    /** Inplace Addition; x += y
     *
     * This function assumes both ring elements have already been 
     * constructed and initialized.
     *
     * @return reference to x.
     * @param  x ring element (reference returned).
     * @param  y ring element.
     */
    Element &addin (Element &x, const Element &y) const
    {
      assume( x );
      assume( x.belongs_to( _coeffs ) );
      assume( y );
      assume( y.belongs_to( _coeffs ) );
      
      return add(x, x, y); // No inplace addition -> TODOs!
    }

    /** Inplace Subtraction; x -= y
     *
     * This function assumes both ring elements have already been 
     * constructed and initialized.
     *
     * @return reference to x.
     * @param  x ring element (reference returned).
     * @param  y ring element.
     */
    Element &subin (Element &x, const Element &y) const
    {
      assume( x );
      assume( x.belongs_to( _coeffs ) );
      assume( y );
      assume( y.belongs_to( _coeffs ) );

      return sub(x, x, y);  // No inplace subtraction!!!?!
    }

    /** Inplace Multiplication; x *= y
     *
     * This function assumes both ring elements have already been 
     * constructed and initialized.
     *
     * @return reference to x.
     * @param  x ring element (reference returned).
     * @param  y ring element.
     */
    Element &mulin (Element &x, const Element &y) const
    {
      assume( x );
      assume( x.belongs_to( _coeffs ) );
      assume( y );
      assume( y.belongs_to( _coeffs ) );
      
      return mul(x, x, y); // TODO: n_InpMult?
    }

   
    /** Inplace Multiplication using a property; x *= y
     *
     * This function does the same thing as above but takes a
     * reference to a Property rather than to an element. This is
     * required for that the scal operation on sparse matrices
     * work.
     *
     * @return reference to x.
     * @param  x Property (reference returned).
     * @param  y ring element.
     */
      template <class Iterator, class Accessor>
          Element &mulin (LELA::Property<Iterator, Accessor> x, const Element &y) const
      {
        assume( x.ref () );
        assume( x.ref ().belongs_to( _coeffs ) );
        assume( y );
        assume( y.belongs_to( _coeffs ) );
        
        return mulin (x.ref (), y);
      }


    /** Inplace Division; x /= y
     *
     * This function attempts to compute the quotient of x and y
     * and, if it exists, stores the result in x. If the quotient
     * does not exist in the ring, the value of x is left
     * unchanged.
     *
     * This function assumes both ring elements have already been 
     * constructed and initialized.
     *
     * @return true if the quotient exists in the ring, false otherwise
     * @param  x ring element (reference returned).
     * @param  y ring element.
     */
    bool divin (Element &x, const Element &y) const
    {
      assume( x );
      assume( x.belongs_to( _coeffs ) );
      assume( y );
      assume( y.belongs_to( _coeffs ) );


      if( !n_DivBy(x, y, _coeffs) )
        return false;
      
      init(x, n_Div(x, y, _coeffs));

      assume( x );
      assume( x.belongs_to( _coeffs ) );
      
      return true;
    }

    /** Inplace Additive Inverse (Inplace Negation).
     * x = - x
     * This function assumes the ring element has already been 
     * constructed and initialized.
     *
     * @return reference to x.
     * @param  x ring element (reference returned).
     */
    Element &negin (Element &x) const
    {
      assume( x );
      assume( x.belongs_to( _coeffs ) );

      if( x.get_counter() == 1 )
      {
        x.raw_set( n_Neg(x, _coeffs) );  // Note: n_Neg creates NO NEW number!
        return x;
      } else
      {
        Element y; copy(y, x);  assume( &y != &x ); // real copy!
        y.raw_set( n_Neg(y, _coeffs) );  // Note: n_Neg creates NO NEW number!
        assume( y );
        assume( y.belongs_to( _coeffs ) );
        return x = y;
      }
    }

    /** Inplace Multiplicative Inverse.
     * x = 1 / x
     * This function assumes the ring element has already been 
     * constructed and initialized.
     *
     * @return true if the inverse of x exists in the ring, false if not
     * @param  x ring element.
     */
    bool invin (Element &x) const
    {
      assume( x );
      assume( x.belongs_to( _coeffs ) );

      Element y;
      
      if( inv(y, x) )
      {
        assume( y );
        assume( y.belongs_to( _coeffs ) );
        
        return x = y;
      }
      
      return false;
      
    }

    /** Inplace AXPY.
     * r  += a * x
     * This function assumes all ring elements have already been 
     * constructed and initialized.
     * @return reference to r.
     * @param  r ring element (reference returned).
     * @param  a ring element.
     * @param  x ring element.
     */
    Element &axpyin (Element &r, const Element &a, const Element &x) const
    {
      assume( r );
      assume( r.belongs_to( _coeffs ) );
      assume( a );
      assume( a.belongs_to( _coeffs ) );
      assume( x );
      assume( x.belongs_to( _coeffs ) );
      
      Element t;
      return addin(r, mul(t, a, x));
      
    }

    //@} Inplace Arithmetic Operations

    /** @name Input/Output Operations */
    //@{

    /** Print ring.
     * @return output stream to which ring is written.
     * @param  os  output stream to which ring is written.
     */
    std::ostream &write (std::ostream &os) const
    {
      StringSetS("");
      n_CoeffWrite(_coeffs);
      const char* s = StringAppendS(""); // TODO!!!

      return os << s;
    }

    /** Read ring.
     * @return input stream from which ring is read.
     * @param  is  input stream from which ring is read.
     */
    std::istream &read (std::istream &is)
    {
      assume(FALSE); // not implemented yet!
      return is;
    }

    /** Print ring element.
     * This function assumes the ring element has already been 
     * constructed and initialized.
     *
     * In this implementation, this means for the <tt>
     * _elem_ptr</tt> for x exists and does not point to
     * null.
     *
     * @return output stream to which ring element is written.
     * @param  os  output stream to which ring element is written.
     * @param  x   ring element.
     */
    std::ostream &write (std::ostream &os, const Element &x) const
    {
      assume( x );
      assume( x.belongs_to( _coeffs ) );
      
      BaseElement t = n_Copy(x, _coeffs); // TODO: this is _REALLY_ ugly!
      assume( n_Test(t, _coeffs) );
      
      StringSetS("");
      n_Write(t, _coeffs); 
      const char* s = StringAppendS("");

      assume( n_Test(t, _coeffs) );
      n_Delete(&t, _coeffs);

      return os << s;
    }

    /** Read ring element.
     * This function assumes the ring element has already been 
     * constructed and initialized.
     *
     * In this implementation, this means for the <tt>
     * _elem_ptr</tt> for x exists and does not point to
     * null.
     *
     * @return input stream from which ring element is read.
     * @param  is  input stream from which ring element is read.
     * @param  x   ring element.
     */
    std::istream &read (std::istream &is, Element &x) const
    {
      int i;
      
      if ( is >> i ) // TODO: this is a workaround... since n_Read needs a string :(
        init(x, i); // and those matrices are with integers anyway!
      else
        assume(FALSE); // not implemented yet!
        
      return is;
    }

    /** Obtain the width in characters of a typical element
     *
     * This can be used to format the output of a matrix in a
     * readable way.
     */
    virtual size_t elementWidth () const
        { LELA::integer c; return (cardinality (c) == 0) ? 10 : (size_t) ceil (log (c.get_d ()) / M_LN10); }
    

    //@} Input/Output Operations
    /** @name Standard elements
     */
    //@{

    /// Return a reference to the zero-element of the ring
    const Element &zero () const
    {
      assume( _zero );      
      assume( _zero.belongs_to( _coeffs ) );
      assume( n_IsZero(_zero, _coeffs) );

      return _zero;
    }

    /// Return a reference to the one-element of the ring
    const Element &one () const
    {
      assume( _one );
      assume( _one.belongs_to( _coeffs ) );      
      assume( n_IsOne(_one, _coeffs) );
      
      return _one;      
    }

    /// Return a reference to the negative of the one-element of the ring
    const Element &minusOne () const
    {
      assume( _minus_one );
      assume( _minus_one.belongs_to( _coeffs ) );
      assume( n_IsMOne(_minus_one, _coeffs) );
      
      return _minus_one;      
    }

      //@}
      //@} Common Object Interface

    /** \brief Conversion of ring element to an integer.  The
     * meaning of conversion is specific to each ring
     * class. However, if x is in the image of the integers in the
     * ring, the integer n returned is such that an init from n
     * will reproduce x. Most often 0 &leq; n &lt; characteristic.
     *
     * @return reference to n.
     * @param n output integer.
     * @param x input ring element.
     */
    int &convert (int &n, const Element &y) const
    {
      assume( y );
      assume( y.belongs_to( _coeffs ) );

      BaseElement t = n_Copy(y, _coeffs);
      assume( n_Test(t, _coeffs) );
      n = n_Int(t, _coeffs);
      assume( n_Test(t, _coeffs) );
      n_Delete(&t, _coeffs);      
      return n;
    }

    LELA::integer &convert (LELA::integer &x, const Element &y) const 
    {
      int n; convert (n, y); return x = LELA::integer(n);      
    }

    double &convert (double &x, const Element &y) const
    {
      int n; convert (n, y); return  x = (double) n;
    }

    float &convert (float &x, const Element &y) const
    {
      int n; convert (n, y); return  x = (float) n;
    }

    inline void cleanup(Element& y) const
    {
      assume( y );
      assume( y.belongs_to( _coeffs ) );

      y.reset();      
    }

    inline bool test(Element& y) const
    {
      assume( y.belongs_to( _coeffs ) );
      return bool(y);
    }

  private:
    
    BaseCoeffs _coeffs;
    bool _clean_coeffs;

    n_coeffType _type; ///< needed for the copy constructor CoeffDomain(CoeffDomain&)
    void* _param;      ///< what about using reference counting of coeffs (themselves)
    
    Element _one; // = n_Init(1, _coeffs);
    Element _zero; // = n_Init(0, _coeffs); // When do we clean these???
    Element _minus_one; // = n_Init(-1, _coeffs);


    
}; // class CoeffDomain


#ifndef HAVE_initializeGMP
#define HAVE_initializeGMP
//#ifdef HAVE_FACTORY
int initializeGMP(){ return 1; }
//#endif
#endif 


#include "lela/blas/level1-generic.h"
#include "lela/blas/level2-generic.h"
#include "lela/blas/level3-generic.h"
#include "lela/blas/level3-sw.h"


#endif /* __LELA_HAVE_LIBPOLYS */

#endif /* __SINGULAR_ring_coeff_H */

