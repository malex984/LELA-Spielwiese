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

#ifdef __LELA_HAVE_LIBPOLYS

#include <lela/ring/interface.h>

#include <omalloc/omalloc.h>
#include <misc/auxiliary.h>

#include <coeffs/coeffs.h>
#include <coeffs/numbers.h>

#include <reporter/reporter.h>
#include <resources/feResource.h>

#include <cassert>



struct ReferenceCountedElement
{
  typedef number SingularNumber;
  typedef coeffs SingularRing;

  ReferenceCountedElement(): m_counter(0), m_element(0), m_coeffs(NULL) {}
  ReferenceCountedElement(SingularNumber c, SingularRing R): m_counter(0), m_element(c), m_coeffs(R) {}  

  long           m_counter; ///< Reference Counter
  SingularNumber m_element; ///< Referenced Element
  SingularRing   m_coeffs;  ///< Element's Coeff. Domain, needed for delete


  static inline void intrusive_ptr_add_ref(ReferenceCountedElement* p)
  {
    assert( p!= NULL );
    assert( p->m_counter >= 0 );

    ++(p->m_counter); 
  }

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

}; // Use omalloc for allocating/deallocating...?


/// @class Number simple smart pointer guarding numbers.
///
/// Its purpose is to count references (in assigns etc.)
/// and destroy the underlying number when it's not needed anymore!
class Number
{

  
  public:
    typedef typename ReferenceCountedElement::SingularNumber SingularNumber;
    typedef typename ReferenceCountedElement::SingularRing   SingularRing;

    Number () : px (NULL) {}

    explicit Number( SingularNumber c, SingularRing R, bool add_ref = true )
    {
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

    void reset()
    {
      Number().swap( *this );
    }

    void reset( SingularNumber c, SingularRing R )
    {
      Number( c, R ).swap( *this );
    }

    // Note: this looks like a workaround... TODO?
    void raw_set( SingularNumber c)
    {
      assert( px != 0 );
      px->m_element = c;
    }    

/*
    SingularNumber const & operator*() const
    {
      assert( px != 0 );
      assume( n_Test(px->m_element, px->m_coeffs) );
      return px->m_element;
    }
    SingularNumber & operator*()
    {
      assert( px != 0 );
//      assume( n_Test(px->m_element, px->m_coeffs) );
      return px->m_element;
    } */
    
    operator SingularNumber() const
    {
      return get();
    }

    operator bool () const
    {
      if( px == 0 )
        return false;
      
      if (px->m_coeffs == 0)
        return false;
      
      return n_Test(px->m_element, px->m_coeffs);
    }

    inline bool belongs_to(SingularRing R) const 
    {
      if (px != 0)
          return (R == px->m_coeffs);
      
      return (R==0);
    }

    void swap(Number & rhs)
    {
      ReferenceCountedElement* tmp = px;
      px = rhs.px;
      rhs.px = tmp;
    }

  private:
    inline SingularRing base_ring() const 
    {
      if (px != 0)

      return (0);
    }

    SingularNumber get() const
    {
      assert( px != 0 );      
      assert( n_Test(px->m_element, px->m_coeffs) );
      return px->m_element;
    }
    

  private:
    ReferenceCountedElement* px; ///< contained pointer, with counter
};

void swap(Number & lhs, Number & rhs)
{
  lhs.swap(rhs);
}


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
      
      assume (R != NULL);

      _coeffs = R;
      _clean_coeffs = cleanup;

      init( _zero, 0 );
      init( _one,  1 );
      init( _minus_one, -1 );
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
    inline Element &init(Element &x, BaseElement y) const
    {
      assume( n_Test(y, _coeffs ) );

      return x = Number(y, _coeffs);
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
    Element &init (Element &x, const int i = 0) const
    {
      return init(x, n_Init(i, _coeffs));
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
    template <class Iterator, class Accessor, class T>
        Element &init (LELA::Property<Iterator, Accessor> x, const T &n) const
    {
      return init(x.ref (), n);
    }

    template <class Iterator, class Accessor>
        Element &init (LELA::Property<Iterator, Accessor> x) const
    {
      return x.ref() = zero();
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
      assume( y.belongs_to( _coeffs ) );
      assume( y );
      
      return init(x, n_Copy(y, _coeffs));
    }

    /// Version of copy which takes a Property rather than an element.
    ///
    /// It should do exactly what the version above taking an element does.
    template <class Iterator, class Accessor>
        Element &copy (LELA::Property<Iterator, Accessor> x, const Element &y) const
    {
      assume( y.belongs_to( _coeffs ) );
      assume( y );
      
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
      assume( x.belongs_to( _coeffs ) );
      assume( x );
      assume( y.belongs_to( _coeffs ) );
      assume( y );
      
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
      assume( y.belongs_to( _coeffs ) );
      assume( z.belongs_to( _coeffs ) );
      assume( y );
      assume( z );
      
      return init(x, n_Add(y, z, _coeffs));
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
      assume( y.belongs_to( _coeffs ) );
      assume( z.belongs_to( _coeffs ) );
      assume( y );
      assume( z );

      return init(x, n_Sub(y, z, _coeffs));
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
      assume( y.belongs_to( _coeffs ) );
      assume( z.belongs_to( _coeffs ) );
      assume( y );
      assume( z );

      return init(x, n_Mult(y, z, _coeffs));
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
      assume( y.belongs_to( _coeffs ) );
      assume( z.belongs_to( _coeffs ) );
      assume( y );
      assume( z );

      if( !n_DivBy(y, z, _coeffs) )
        return false;

      init(x, n_Div(y, z, _coeffs));
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
      assume( y.belongs_to( _coeffs ) );
      assume( y );

      copy(x, y); assume( &y != &x ); // real copy!
      x.raw_set( n_Neg(x, _coeffs) ); // Note: n_Neg creates NO NEW number!
      
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
      assume( y.belongs_to( _coeffs ) );
      assume( y );

      if( !n_DivBy(one(), y, _coeffs) )
        return false;

      init(x, n_Invers(y, _coeffs));
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
      assume( a.belongs_to( _coeffs ) );
      assume( a );
      assume( x.belongs_to( _coeffs ) );
      assume( x );
      assume( y.belongs_to( _coeffs ) );
      assume( y );
      
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
      assume( x.belongs_to( _coeffs ) );
      assume( x );
      
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
      assume( x.belongs_to( _coeffs ) );
      assume( x );

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
      assume( x.belongs_to( _coeffs ) );
      assume( x );
      assume( y.belongs_to( _coeffs ) );
      assume( y );


      if( !n_DivBy(x, y, _coeffs) )
        return false;
      
      init(x, n_Div(x, y, _coeffs));

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
      assume( x.belongs_to( _coeffs ) );
      assume( x );
      
      Element y; copy(y, x);  assume( &y != &x ); // real copy!
      y.raw_set( n_Neg(y, _coeffs) );  // Note: n_Neg creates NO NEW number!
      
      return x = y;
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
      assume( x.belongs_to( _coeffs ) );
      assume( x );

      Element y;
      
      if( inv(y, x) )
        return x = y;
      
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
      assume( x.belongs_to( _coeffs ) );
      assume( x );
      
      // TODO: writing an undefined element / QQ!!! :(((

      // TODO: this is _REALLY_ ugly!
      BaseElement t = n_Copy(x, _coeffs);
      assume( n_Test(t, _coeffs) );
      
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
      assume( _zero.belongs_to( _coeffs ) );
      assume( _zero );      
      assume( n_IsZero(_zero, _coeffs) );

      return _zero;
    }

    /// Return a reference to the one-element of the ring
    const Element &one () const
    {
      assume( _one.belongs_to( _coeffs ) );
      assume( _one );
      
      assume( n_IsOne(_one, _coeffs) );
      
      return _one;      
    }

    /// Return a reference to the negative of the one-element of the ring
    const Element &minusOne () const
    {
      assume( _minus_one.belongs_to( _coeffs ) );
      assume( _minus_one );

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
      assume( y.belongs_to( _coeffs ) );
      assume( y );

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
      assume( y.belongs_to( _coeffs ) );
      assume( y );

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

