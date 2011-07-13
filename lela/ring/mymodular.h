/* lela/ring/mymodular.h
 * Copyright (C) 2011 Oleksandr Motsak
 *
 * Written by Oleksand Motsak <http://goo.gl/mcpzY>
 *
 * ------------------------------------
 *
 * BSD
 */

#ifndef __SINGULAR_ring_modular_H
#define __SINGULAR_ring_modular_H


#include <lela/lela-config.h>
#include <lela/ring/coeffs.h>

#ifdef __LELA_HAVE_LIBPOLYS

template <class _Element> // Just ignore this template parameter
class MyModular: public CoeffDomain
{
  typedef CoeffDomain Base;
  
  void InitMe(unsigned long modulus)
  {
    if( _coeffs != NULL )
      if( _clean_coeffs )
        nKillChar(_coeffs);

    _coeffs = nInitChar( n_Zp, (void*) modulus );
    assume( _coeffs != NULL);
    _clean_coeffs = true;
  }

  public:    
    typedef typename Base::Element Element;

    class RandIter;
    
    MyModular (unsigned long _modulus): Base()
    {
      InitMe( _modulus );
    }

    MyModular (const LELA::integer &modulus): Base()
    {
      const unsigned long _modulus = modulus.get_ui();
      InitMe( _modulus );
    }

    MyModular (const MyModular<Element> &F): Base()
    {
      int _modulus;
      F.cardinality(_modulus);
      InitMe(_modulus);
    }
};

// #include <lela/randiter/modular.h>

#include <iostream>
#include <vector>

#include <lela/randiter/mersenne-twister.h>
#include <lela/util/commentator.h>

template <class E>
class MyModularRandIter
{
  public:
    typedef typename MyModular<E>::Element Element;

  /** Constructor from field, sampling size, and seed.
	 * The random field element iterator works in the field F, is seeded
	 * by seed, and it returns any one element with probability no more
	 * than 1/min (size, F.cardinality (c)).
	 * A sampling size of zero means to sample from the entire field.
	 * A seed of zero means to use some arbitrary seed for the generator.
	 * Purely virtual.
	 * @param F LinBox field archetype object in which to do arithmetic
	 * @param size constant integer reference of sample size from which to 
	 *             sample (default = modulus of field)
	 * @param seed constant integer reference from which to seed random number
	 *             generator (default = 0)
	 */
    MyModularRandIter (const MyModular<E> &F, 
                     const LELA::integer &size = 0, 
                     const LELA::integer &seed = 0)
        : _MT (seed.get_ui ()), _F (F), _size (size), _seed (seed.get_ui ())
    {
      LELA::integer cardinality;

      F.cardinality (cardinality);

      if ((_size == 0) || (_size > cardinality))
        _size = cardinality;

      LELA::commentator.report (10, INTERNAL_DESCRIPTION)
          << "Created random generator with size " << _size 
          << " and seed " << _seed << std::endl;
    }

  /** Copy constructor.
	 * Constructs ModularRandIter object by copying the random field
	 * element generator.
	 * This is required to allow generator objects to be passed by value
	 * into functions.
	 * @param  R ModularRandIter object.
	 */
    MyModularRandIter (const MyModularRandIter<E> &R) 
        : _F (R._F), _size (R._size), _seed (R._seed) {}

  /** Destructor.
	 * This destructs the random field element generator object.
	 */
    ~MyModularRandIter () {}

  /** Assignment operator.
	 * Assigns ModularRandIter object R to generator.
	 * @param  R ModularRandIter object.
	 */
    MyModularRandIter<E> &operator = (const MyModularRandIter<E> &R)
    {
      if (this != &R) { // guard against self-assignment
        _size = R._size;
        _seed = R._seed;
        _MT.setSeed (_seed);
      }

      return *this;
    }

  /** Random field element creator.
	 * This returns a random field element from the information supplied
	 * at the creation of the generator.
	 * Required by abstract base class.
	 * @return reference to random field element
	 */
    unsigned int &random (unsigned int &a) const
    {
      return a = _MT.randomIntRange (0, _size.get_ui ());
    }
    Element &random (Element &a) const
    {
      unsigned int n; random(n);
      return _F.init (a, (int)n);
    }

  private:
    LELA::MersenneTwister _MT;

  /// Field in which arithmetic is done
    MyModular<E> _F;

  /// Sampling size
    LELA::integer _size;

  /// Seed
    long _seed;

}; // class ModularRandIter


template <class E>
class MyModular<E>::RandIter
{
  MyModularRandIter<E> _r;

  public:
    RandIter (const MyModular<E> &F, const LELA::integer &size = 0, const LELA::integer &seed = 0)
        : _r (F, size, seed) {}
    RandIter (const MyModular<E>::RandIter &r)
        : _r (r._r) {}

    ~RandIter () {}
    RandIter &operator= (const RandIter &r)
                        { _r = r._r; return *this; }

    typedef typename MyModular<E>::Element Element;

    // Workaround?
    unsigned int &random (unsigned int &a) const
        { return _r.random (a); }
    

    Element &random (Element &a) const
        { return _r.random (a); }
};


#else

// The following should be the only inclusion of the old stuff!
#include <lela/ring/old.modular.h>
   
template <class _Element> // Just ignore this template parameter
class MyModular: public LELA::Modular<_Element>
{
  typedef LELA::Modular<_Element> Base;
  
public:
  typedef typename Base::Element Element;

  MyModular (unsigned long modulus): Base(modulus) {}

  MyModular (const LELA::integer &modulus): Base(modulus) {}

  MyModular (const MyModular<Element> &F): Base(F) {}
};


#endif /* __LELA_HAVE_LIBPOLYS */

#endif /* __SINGULAR_ring_modular_H */
