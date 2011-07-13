/* tests/test-coeffs.C
 * Copyright 2011 Oleksandr Motsak
 * Written by Oleksand Motsak <http://goo.gl/mcpzY>
 *
 * ---------------------------------------------------------
 *
 * BSD
 *
 * Test for Singular Coefficient Domains
 */

// #define DISABLE_LIBPNG 1
// #define DISABLE_COMMENTATOR 1

#include "test-common.h"

#include <lela/ring/coeffs.h>
#include <lela/ring/mymodular.h>

#include <lela/util/commentator.h>

#include <lela/vector/sparse.h>
#include <lela/vector/stream.h>

#include <lela/matrix/dense.h>
#include <lela/matrix/sparse.h>

#include <lela/blas/context.h>

#include <lela/blas/level1.h>
#include <lela/blas/level2.h>
#include <lela/blas/level3.h>

#include <lela/solutions/echelon-form.h>

using namespace LELA;

#ifndef __LELA_HAVE_LIBPOLYS
#error Cannot be used without libpolys
#endif 



#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <cstdio>
#include <cmath>

using namespace std;



// #include <polys/monomials/ring.h>
// #include <polys/monomials/p_polys.h>

template <class Ring>
bool  MyArithTest (const Ring &R) //  , typename Ring::Element &v  
{
  commentator.start ("MyArithTest<Ring>(R)", __FUNCTION__);

  ostream &report = commentator.report (Commentator::LEVEL_NORMAL, INTERNAL_DESCRIPTION);  // TDO: ASK BRAD!!!! 
  ostream &error = commentator.report (Commentator::LEVEL_IMPORTANT, INTERNAL_ERROR);


  Context<Ring> ctx(R);

  typedef typename Ring::Element MYElement; // :((((

  report << (":::::::::::::::: TESTING LB :::::::::::::::!");
  report << endl;
  R.write(report << "Ring: " << endl ) << endl << flush ;

  integer i(0);  
  report << "Ring'  char: " << R.characteristic(i).get_ui() << endl << flush;
  report << "Ring'  card: " << (int)R.cardinality(i).get_ui() << endl << flush;

  R.write(report << "Ring'  0: ", R.zero()) << endl << flush;
  R.write(report << "Ring' +1: ", R.one()) << endl << flush;
  R.write(report << "Ring' -1: ", R.minusOne()) << endl << flush;

  typename Ring::Element a, b, sum, mul, div_a, div_b;

  const int _a = 33;
  const int _b = 103;

  R.add (sum, R.init (a, _a), R.init (b, _b));
  R.mul (mul, a, b);

  report << "a: " << flush; R.write (report, a) << " (" << a << ") + " << flush;
  report << "b: " << flush; R.write (report, b) << " (" << b << ") --->>>>> " << flush;

  report << "sum: " << flush; R.write (report, sum) << ", " << sum << endl << flush;
  report << "mul: " << flush; R.write (report, mul) << ", " << mul << endl << flush;

  bool _test_divide = false;
  if( R.isZero(a) || R.isZero(b) )
  {
    assume( R.isZero(mul) );    
  } else
  {
    _test_divide = R.div (div_a, mul, a); assume(_test_divide);
    report << "div: mul/a: " << flush; R.write (report, div_a) << ", " << div_a << endl << flush;

    _test_divide = R.div (div_b, mul, b); assume(_test_divide);
    report << "div: mul/b: " << flush; R.write (report, div_b) << ", " << div_b << endl << flush;

  }


  commentator.start ("Testing Addition: ");

  report << "" << flush; R.write (report, R.init (a, _a + _b) ) << endl << flush;
  report << "ZeroDifference?: "  << flush; R.write (report, R.subin (a, sum) )  << endl << flush;

  report << ("Sum Test: ");

  commentator.stop (MSG_STATUS (R.isZero(a)));

  commentator.start ("Testing Mult: ");

  report << "Mult: " << flush; R.write (report, R.init (a, _a * _b) ) << endl << flush;
  report << "ZeroDifference?: "  << flush; R.write (report, R.subin (a, mul) )  << endl << flush;

  report << ("Mul Test: ");
  commentator.stop (MSG_STATUS (R.isZero(a)));

  if( _test_divide )
  {
    commentator.start ("Testing Div_a: ");

    report << "Division by a: " << flush; R.write (report, R.init (a, (_a * _b) / _a ) ) << endl << flush;
    report << "ZeroDifference?: "  << flush; R.write (report, R.subin (a, div_a) )  << endl << flush;

    report << ("Div_a Test: ");
    commentator.stop (MSG_STATUS (R.isZero(a)));


    commentator.start ("Testing Div_b: ");

    report << "Division by b: " << flush; R.write (report, R.init (b, (_a * _b) / _b ) ) << endl << flush;
    report << "ZeroDifference?: "  << flush; R.write (report, R.subin (b, div_b) )  << endl << flush;

    report << ("Div_b Test: ");
    commentator.stop (MSG_STATUS (R.isZero(b)));
  }


  if(true)
  {
    commentator.start ("Testing LA: ");
    size_t N = 5;

    // This business the 'typename Ring::Element' is just ugly and unavoidable :( 
    typedef DenseMatrix<typename Ring::Element> MDense;
    typedef SparseMatrix<typename Ring::Element> MSparse;

    typedef typename Vector<Ring>::Dense  VDense;
    typedef typename Vector<Ring>::Sparse VSparse;
    //    typedef typename Vector<Ring>::Hybrid VHybrid; // No Hybrid per default :(

    typedef typename VSparse::value_type  VValue;

    VSparse w; VDense  v(2*N);

    MDense  D(2*N, 2*N); MSparse S(2*N, 2*N);

    for ( size_t i = 0; i < N; i++ )
    {
      typename Ring::Element t;
      R.init(t, (const int)(17 - 3 * i));
      w.push_back( VValue(i, t) ); 
      v[i] = t; // OK? TODO: ASK Brad!

      //      D.setEntry ( (11 + 719*i) % (2*N), (71*i + 213) % (2*N), t);
      S.setEntry ( i/2, i, t);
      // No cleanup for t?
    }

    // And backwards...:
    for ( size_t i = 2*N-1; i >= N; i-- )
    {
      typename Ring::Element t;
      R.init(t, (const int)(7 * i - 41));

      R.init(v[i], (const int)(37 * i));
      //      assert v[i] = t; // OK? TODO: ASK Brad!

      w.push_back( VValue(i, t) ); // NO WAY!!! -> ASK Brad!


      //      D.setEntry ( (11 - 719*i) % (2*N), (113 - 17*i) % (2*N), t);
      S.setEntry ( i/2, i, t);

      // No cleanup for t?
    }

    BLAS1::write(ctx, report << "v: ", v) << endl << flush;
    BLAS1::write(ctx, report << "w: ", w) << endl << flush;



    report << "Matrices: "<< endl << flush;

    report << "D: " << D.coldim () << " x " << D.rowdim () << endl << flush;
    report << "S: " << S.coldim () << " x " << S.rowdim () << endl << flush;

    // The following need libpng???!!! :(((
    BLAS3::write(ctx, report << "S: " << endl, S, FORMAT_SAGE) << endl << flush;


    //    MDense  DD(2*N, 2*N); MSparse SS(2*N, 2*N);


    // The following is buggy for GF5^2 ('1' for zeroes) and somewhat for QQ ('o'!?)
    BLAS3::write(ctx, report << "D = S: " << endl, BLAS3::copy(ctx, S, D), FORMAT_SAGE) << endl << flush;
    //    BLAS3::write(ctx, report << "SS = D: ", BLAS3::copy(ctx, D, SS), FORMAT_SAGE) << endl << flush;




    VSparse ww; BLAS1::write(ctx, report << "ww = v: ", BLAS1::copy(ctx, v, ww)) << endl << flush;

    //:(//    VDense  vv(2*N); BLAS1::write(ctx, report << "vv = w: ", BLAS1::copy(ctx, w, vv)) << endl << flush;

    EchelonForm<Ring> EF(ctx);

    // The following is Ok (where D is ~fine), and Seg.faults with QQ
    // :(
    try
    {
      BLAS3::write(ctx, report << "red RowEchelonForm(D): " << endl, EF.RowEchelonForm(D, true), FORMAT_SAGE);
    }
    catch (LinboxError e)
    {
      commentator.report (Commentator::LEVEL_IMPORTANT, INTERNAL_ERROR) << e;
    }
    catch (...)
    {
      commentator.report (Commentator::LEVEL_IMPORTANT, INTERNAL_ERROR)
          << "No idea.... :(" << std::endl;
    }

    try
    {
      BLAS3::write(ctx, report << "red RowEchelonForm(S): " << endl, EF.RowEchelonForm(S, true), FORMAT_SAGE); 
    }
    catch (LinboxError e)
    {
      commentator.report (Commentator::LEVEL_IMPORTANT, INTERNAL_ERROR) << e;
    }
    catch (...)
    {
      commentator.report (Commentator::LEVEL_IMPORTANT, INTERNAL_ERROR)
          << "No idea.... :(" << std::endl;
    }


    commentator.stop (MSG_STATUS (TRUE));

  }


  //  report << flush;

  commentator.stop ("done", (const char *) 0, "MyArithTest");

  return TRUE;

}



bool TestLinbox(coeffs r)
{
  return MyArithTest(CoeffDomain(r)); 
}


BOOLEAN Test(const n_coeffType type, void* param)
{
  BOOLEAN ret = TRUE;

  const coeffs r = nInitChar( type, param );

  if( r == NULL )
  {
    PrintS( "Test: could not get this coeff. domain" ); PrintLn();
    return FALSE;
  }

  assume( r->cfCoeffWrite != NULLp );

  if( r->cfCoeffWrite != NULL )
  {
    PrintS( "Coeff-domain: " ); n_CoeffWrite(r); PrintLn();
  }

  number _zero = n_Init(0, r);
  assume( n_Test(_zero, r) );
  assume( n_IsZero(_zero, r) );
  n_Delete(&_zero, r);


  number _one = n_Init(1, r);
  assume( n_Test(_one, r) );
  assume( n_IsOne(_one, r) );
  n_Delete(&_one, r);


  number _minus_one = n_Init(-1, r);
  assume( n_Test(_minus_one, r) );
  assume( n_IsMOne(_minus_one, r) );
  n_Delete(&_minus_one, r);





  number t = n_Init(1, r);  
  ndInpAdd(t, t, r);  

  number two = n_Init(2, r);
  assume(n_Equal(two, t, r));
  ret = ret && n_Equal(two, t, r);
  n_Delete(&t, r);
  n_Delete(&two, r);




  if( !TestLinbox(r) )
    return FALSE; 

/*

  const int N = 2; // number of vars
  char* n[N] = {"x", "y"}; // teir names

  ring R = rDefault( r, N, n);  // now r belongs to R!

  if( R == NULL )
  {
    PrintS( "Test: could not get a polynomial ring over this coeff. domain" ); PrintLn();
    nKillChar( r );
    return FALSE;
  }


  rWrite(R); PrintLn();
  #define RDEBUG
  #ifdef  RDEBUG
//    rDebugPrint(R); PrintLn();
  #endif

  const int exp[N] = {5, 8};

  poly p = p_ISet(1, R);
  assume( p != NULL );
  assume(pNext(p) == NULL);
  ret = ret && (pNext(p) == NULL);

  p_SetExp(p,1,exp[0],R);
  p_SetExp(p,2,exp[1],R);
  p_Setm(p, R);

  assume( p_GetExp(p, 1, R) == exp[0] );
  assume( p_GetExp(p, 2, R) == exp[1] );

  p = p_Add_q(p_Copy(p, R), p, R);  

  ret = ret && p_EqualPolys(p, p, R);


  poly p2 = p_ISet(2, R);
  assume( p2 != NULL );
  assume(pNext(p2) == NULL);
  ret = ret && (pNext(p2) == NULL);

  p_SetExp(p2,1,exp[0],R);
  p_SetExp(p2,2,exp[1],R);
  p_Setm(p, R);

  assume( p_GetExp(p2, 1, R) == exp[0] );
  assume( p_GetExp(p2, 2, R) == exp[1] );

  number s = p_GetCoeff(p, R);
  two = n_Init(2, r);
  assume(n_Equal(two, s, r));
  ret = ret && n_Equal(s, two, r);  
  n_Delete(&two, r);

  assume(p_EqualPolys(p2, p, R));
  ret = ret && p_EqualPolys(p, p, R);

  p_Delete(&p, R);

  p_Delete(&p2, R);




  rDelete(R);
*/

  return ret;
}


BOOLEAN simple(const string tn, const n_coeffType _type, void* param = NULL)
{
  bool pass = true;
  string title = "MyArithTest ( CoeffDomain(_type, param) ) for [" + tn + "]";
  
  commentator.start (title.c_str(), __FUNCTION__);
  try
  {
    pass = Test(_type, param);
  }
  catch (LinboxError e)
  {
    commentator.report (Commentator::LEVEL_IMPORTANT, INTERNAL_ERROR) << e;
  }
  catch (...)
  {
    commentator.report (Commentator::LEVEL_IMPORTANT, INTERNAL_ERROR)
        << "No idea.... :(" << std::endl;
  }
  commentator.stop (MSG_STATUS (pass)); 
  
  return pass;
}

bool TestLinbox()
{
  bool pass = true;
   
  typedef MyModular<uint8> Ring;
  typedef Ring::Element Element;
   
  Ring F(5);

  commentator.start ("MyArithTest( MyModular<uint8>(5) )", __FUNCTION__);
  try
  {
    pass = MyArithTest(F);
  }
  catch (LinboxError e)
  {
    commentator.report (Commentator::LEVEL_IMPORTANT, INTERNAL_ERROR) << e;
  }
  catch (...)
  {
    commentator.report (Commentator::LEVEL_IMPORTANT, INTERNAL_ERROR)
        << "No idea.... :(" << std::endl;
  }
    
  commentator.stop (MSG_STATUS (pass)); 

  return pass;
}









//#ifdef HAVE_FACTORY
//int initializeGMP(){ return 1; }
//#endif

int main (int argc, char **argv)
{
  bool pass = true;

  static Argument args[] = {
    { '\0' }
  };

  parseArguments (argc, argv, args);

  commentator.setBriefReportParameters (Commentator::OUTPUT_CONSOLE, false, false, false);
  commentator.getMessageClass (INTERNAL_DESCRIPTION).setMaxDepth (5);
  commentator.getMessageClass (INTERNAL_DESCRIPTION).setMaxDetailLevel (Commentator::LEVEL_UNIMPORTANT);
  commentator.getMessageClass (TIMING_MEASURE).setMaxDepth (3);

//  commentator.setBriefReportStream(cout);
//  commentator.setDefaultReportFile("ForBrad.log");

  commentator.start ("Singular Coeffs", "");

//  pass = testAdd () && pass;

  feInitResources(argv[0]);

  StringSetS("ressources in use (as reported by feStringAppendResources(0):\n");
  feStringAppendResources(0);
  PrintS(StringAppendS("\n"));


  pass = TestLinbox() && pass;

  // modulop
  pass = simple("F7", n_Zp, (void*)7) && pass;

  // due to coeffs/ffields.h
  // TODO: remedy this!
  struct 
  {
    int GFChar;
    int GFDegree;
    char* GFPar_name;
  } param;

  param.GFChar= 5;
  param.GFDegree= 2;
  param.GFPar_name= "Q";

  pass = simple("GF", n_GF, (void*)&param) && pass;
  
  // longrat
  pass = simple("QQ", n_Q) && pass;

  commentator.stop (MSG_STATUS (pass));

  return pass ? 0 : -1;
}
