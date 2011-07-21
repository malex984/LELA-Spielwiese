/* lela/ring/coeffs.h
 * Copyright (C) 2011 Oleksandr Motsak
 *
 * Written by Oleksand Motsak <http://goo.gl/mcpzY>
 *
 * ------------------------------------
 *
 * BSD
 */

#ifndef __SINGULAR_ring_coeff_H
#define __SINGULAR_ring_coeff_H

#include <lela/lela-config.h>
#include <lela/ring/interface.h>

#ifdef __LELA_HAVE_LIBPOLYS

#include <omalloc/omalloc.h>
#include <misc/auxiliary.h>

#include <coeffs/coeffs.h>
#include <coeffs/numbers.h>

#include <reporter/reporter.h>
#include <resources/feResource.h>

/** Singular coefficient domain
 *
 * This class defines the ring-interface for general Singular coefficients.
 */

class CoeffDomain
{
  private:
    CoeffDomain& operator=(const CoeffDomain&);

  protected:

    void CleanMe()
    {
      if (_coeffs != NULL)
      {
        assume( n_Test(_zero, _coeffs) );
        assume( n_IsZero(_zero, _coeffs) );

        assume( n_Test(_one, _coeffs) );
        assume( n_IsOne(_one, _coeffs) );

        assume( n_Test(_minus_one, _coeffs) );
        assume( n_IsMOne(_minus_one, _coeffs) );

        n_Delete(&_zero, _coeffs);
        n_Delete(&_one, _coeffs);
        n_Delete(&_minus_one, _coeffs);

        if (_clean_coeffs)
          nKillChar(_coeffs);
        
        _type = n_unknown;
        _param = NULL;        
      }
    }

    void InitMe(coeffs R, bool cleanup)
    {
      CleanMe();
      
      assume (R != NULL);

      _coeffs = R;

      _zero = n_Init(0, _coeffs);
      _one  = n_Init(1, _coeffs);
      _minus_one = n_Init(-1, _coeffs);

      _clean_coeffs = cleanup;
      

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
    typedef number Element;

//TODO!//    /// An object of this type is a generator of random ring elements.
//TODO!//    typedef RandIterInterface RandIter;

    /// @name Object Management
    //@{

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
    Element &init (Element &x, const int i = 0) const
    {
      x = n_Init(i, _coeffs);
//      assume( x != (Element)NULL);
      assume( n_Test(x, _coeffs) );
      return x;
    }

    Element &init (Element &x, const unsigned int &y) const
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
    template <class Iterator, class Accessor, class T>
        Element &init (LELA::Property<Iterator, Accessor> x, const T &n = 0) const
        { return init (x.ref (), n); }


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
      return x = n_Copy(y, _coeffs);
    }

    /// Version of copy which takes a Property rather than an element.
    ///
    /// It should do exactly what the version above taking an element does.
    template <class Iterator, class Accessor>
        Element &copy (LELA::Property<Iterator, Accessor> x, const Element &y) const
        { return copy (x.ref (), y); }

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
      assume( n_Test(x, _coeffs) );
      assume( n_Test(y, _coeffs) );
      
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
      assume( n_Test(y, _coeffs) );
      assume( n_Test(z, _coeffs) );
      
      x = n_Add(y, z, _coeffs);
      assume( n_Test(x, _coeffs) );
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
      assume( n_Test(y, _coeffs) );
      assume( n_Test(z, _coeffs) );

      x = n_Sub(y, z, _coeffs);
      assume( n_Test(x, _coeffs) );
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
      assume( n_Test(y, _coeffs) );
      assume( n_Test(z, _coeffs) );

      x = n_Mult(y, z, _coeffs);
      assume( n_Test(x, _coeffs) );
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
      assume( n_Test(y, _coeffs) );
      assume( n_Test(z, _coeffs) );
      
      if( !n_DivBy(y, z, _coeffs) )
        return false;

      x = n_Div(y, z, _coeffs);
      assume( n_Test(x, _coeffs) );
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
      assume( n_Test(y, _coeffs) );
      copy(x, y);
      return negin(x);      
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
      assume( n_Test(y, _coeffs) );     

      if( !n_DivBy(one(), y, _coeffs) )
        return false;

      x = n_Invers(y, _coeffs);
      assume( n_Test(x, _coeffs) );      
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
      assume( n_Test(a, _coeffs) );     
      assume( n_Test(x, _coeffs) );     
      assume( n_Test(y, _coeffs) );

      Element t = n_Copy(a, _coeffs);
      n_InpMult(t, x, _coeffs);
      assume( n_Test(t, _coeffs) );
      
      r = n_Add(t, y, _coeffs);

      n_Delete(&t, _coeffs);

      assume( n_Test(r, _coeffs) );
      return r;
      
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
      assume( n_Test(x, _coeffs) );     
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
      assume( n_Test(x, _coeffs) );     
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
      assume( n_Test(y, _coeffs) );     
      Element z = n_Add(x, y, _coeffs); // No inplace addition -> TODOs!
      n_Delete(&x, _coeffs);
      assume( n_Test(z, _coeffs) );     
      return x = z;
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
      assume( n_Test(y, _coeffs) );     
      Element z = n_Sub(x, y, _coeffs); // No inplace subtraction!!!?!
      n_Delete(&x, _coeffs);
      assume( n_Test(z, _coeffs) );     
      return x = z;
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
      assume( n_Test(x, _coeffs) );     
      assume( n_Test(y, _coeffs) );
      
      n_InpMult(x, y, _coeffs);
      assume( n_Test(x, _coeffs) );     
      
      return x;
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
          { return mulin (x.ref (), y); }


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
      assume( n_Test(x, _coeffs) );     
      assume( n_Test(y, _coeffs) );

      if( !n_DivBy(x, y, _coeffs) )
        return false;
      
      Element r = n_Div(x, y, _coeffs);
      assume( n_Test(r, _coeffs) );
      
      n_Delete(&x, _coeffs);      
      x = r;

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
      assume( n_Test(x, _coeffs) );
      
      x = n_Neg(x, _coeffs);
      assume( n_Test(x, _coeffs) );     

      return x;
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
      assume( n_Test(x, _coeffs) );
      
      if( !n_DivBy(one(), x, _coeffs) )
        return false;

      Element r = n_Div(one(), x, _coeffs);
      assume( n_Test(r, _coeffs) );

      n_Delete(&x, _coeffs);      
      x = r;
      return true;
      
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
      assume( n_Test(r, _coeffs) );     
      assume( n_Test(a, _coeffs) );
      assume( n_Test(x, _coeffs) );

      Element t = n_Copy(a, _coeffs);
      assume( n_Test(t, _coeffs) );
      
      n_InpMult(t, x, _coeffs);
      assume( n_Test(t, _coeffs) );
      Element s = n_Add(r, t, _coeffs);
      assume( n_Test(s, _coeffs) );
      
      n_Delete(&t, _coeffs);
      
      n_Delete(&r, _coeffs);      
      return r = s;
      
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
//      assume( x != (Element)NULL ); // TODO: writing an undefined element / QQ!!! :(((
      assume( n_Test(x, _coeffs) ); // undefined => error!

      // TODO: this is _REALLY_ ugly!
      Element t = n_Copy(x, _coeffs);
//      assume( t != (Element)NULL ); // yeah, copy of NULL is NULL :(
      assume( n_Test(t, _coeffs) ); // thus: undefined!
      
      StringSetS("");
      n_Write(t, _coeffs); 
      const char* s = StringAppendS("");

      assume( n_Test(t, _coeffs) ); // and, yes this will fail as well 
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
      assume( n_Test(_zero, _coeffs) );
      assume( n_IsZero(_zero, _coeffs) );

      return _zero;
    }

    /// Return a reference to the one-element of the ring
    const Element &one () const
    {
      assume( n_Test(_one, _coeffs) );
      assume( n_IsOne(_one, _coeffs) );
      
      return _one;      
    }

    /// Return a reference to the negative of the one-element of the ring
    const Element &minusOne () const
    {
      assume( n_Test(_minus_one, _coeffs) );
      assume( n_IsMOne(_minus_one, _coeffs) );
      
      return _minus_one;      
    }

    /** @name Reference-counting of elements
     *
     * These functions need not be implemented, but are useful for
     * memory-management with elements.
     */
    //@{

    virtual void ref (Element &x) const {}
    virtual void unref (Element &x) const {}
    
    //@}

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
        assume( n_Test(y, _coeffs) );
        Element t = n_Copy(y, _coeffs);
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

    //@}
    //@} Common Object Interface

      void cleanup(Element& v) const
      {
        assume( n_Test(v, _coeffs) );
        n_Delete(&v, _coeffs);
      }


      bool test(Element& v) const
      {
        return n_Test(v, _coeffs);
      }
    


  private:
    coeffs _coeffs;
    bool _clean_coeffs;

    n_coeffType _type;
    void* _param;
    
    Element _one; // = n_Init(1, _coeffs);
    Element _zero; // = n_Init(0, _coeffs); // When do we clean these???
    Element _minus_one; // = n_Init(-1, _coeffs);


    
}; // class RingInterface




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

