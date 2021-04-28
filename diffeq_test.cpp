//  file: diffeq_test.cpp
// 
//  Program to study the error in differential equation algorithms
//
//  Programmer:  Jiaqi Zhang zhang.9203@osu.edu
//
//  Revision history:
//      02/07/04  original version, translated from diffeq_test.c
//      01/30/05  comments improved and function names changed
//      01/22/06  improved code organization
//      04/28/20  add rk2 and rk45 for final project
//
//
//********************************************************************

// include files
#include <iostream>		// note that .h is omitted
#include <iomanip>		// note that .h is omitted
#include <fstream>		// note that .h is omitted
using namespace std;		// we need this when .h is omitted
#include <cmath>
#include "diffeq_routines.h"	// diffeq routine prototypes

// function prototypes
double rhs (double t, double y[], int i, void *params_ptr);
double exact_answer (double t, void *params_ptr);

// structures
typedef struct			// define a type to hold parameters 
{
  double alpha;
  double beta;
}
f_parameters;			// now we can define a structure of this type
				//   using the keyword "f_parameters" 

//************************** main program ****************************
int
main (void)
{
  void *params_ptr;		   // void pointer passed to functions 
  f_parameters funct_parameters;   // parameters for the function 

  const int N = 1;		// size of arrays of y functions
  double y_euler[N],y_rk2[N], y_rk4[N],y_rk45[N];	// arrays of y functions

  ofstream out1 ("diffeq_test1.dat");	// open the output file
    ofstream out ("diffeq_test.dat");    // open the output file
    ofstream out2 ("diffeq_test2.dat");    // open the output file
  funct_parameters.alpha = 1.;	// function parameter to be passed 
  funct_parameters.beta = 1.;	// function parameter to be passed
  params_ptr = &funct_parameters;	// structure to pass to function 
    out << "#      t           y_euler(t)     y_rk2(t)     y_rk4(t)      y_rk45(t)\n";
  double tmin = 0.;		// starting t value 
  double tmax = 3.;		// last t value


  // print out a header line and the first set of points 
  out1 << "#     t         y_euler(t)    y_rk2(t)     y_rk4(t)       y_rk45(t)       y_exact(t)         \n";
  out1 << scientific << setprecision (16)
    << tmin << "  "
    << y_euler[0]  << "  "
    << y_rk2[0]  << "  "
    << y_rk4[0]<<"    "<<y_rk45[0]  << "  " << exact_answer (tmin, params_ptr) << endl;



  double h = 0.1;		// initialize mesh spacing
//    double hmax = h*64.;
//    double hmin = h/64.;
      // initial condition for y(t)
    while (h >= 0.00001)
    {
        y_euler[0] = 1.0;        // initial condition for y(t)
          y_rk2[0] = 1.0;        // initial condition for y(t)
          y_rk4[0] = 1.0;        // initial condition for y(t)
          y_rk45[0] = 1.0;        // initial condition for y(t)
    for (double t = tmin; t <= tmax; t += h)
  {

    // find y(t+h) and output vs. t+h1 
    euler (N, t, y_euler, h, rhs, params_ptr);	// Euler's algorithm
    runge2 (N, t, y_rk2, h, rhs, params_ptr);            // 2th order R-K
    runge4 (N, t, y_rk4, h, rhs, params_ptr);	        // 4th order R-K
      runge45 (N, t, y_rk45, h, rhs, params_ptr);            // R-K-F
    out1 << scientific << setprecision (16)
      << t + h << "  "
      << y_euler[0] << "  "
      << y_rk2[0] << "  "
      << y_rk4[0] << "  "<< y_rk45[0] << "  " << exact_answer (t+h, params_ptr) << endl;

  

            out << scientific << setprecision (16)
                          << t + h << "   "
                          << fabs ((y_euler[0] - exact_answer (t + h, params_ptr)) / exact_answer (t + h, params_ptr)) << "   "
                          << fabs ((y_rk2[0] - exact_answer (t + h, params_ptr)) / exact_answer (t + h, params_ptr)) << "   "
            << fabs ((y_rk4[0] - exact_answer (t + h, params_ptr)) / exact_answer (t + h, params_ptr))<< "   "
                          << fabs ((y_rk45[0] - exact_answer (t + h, params_ptr)) / exact_answer (t + h, params_ptr))<< endl;
//
//
//  }
//            if(fabs(t - 2.) < 0.000000001)
//            {
//                out2 << scientific << setprecision (16)
//                              << t << "   "
//                              << fabs ((y_euler[0] - exact_answer (t + h, params_ptr)) / exact_answer (t + h, params_ptr)) << "   "
//                              << fabs ((y_rk2[0] - exact_answer (t + h, params_ptr)) / exact_answer (t + h, params_ptr)) << "   "
//                << fabs ((y_rk4[0] - exact_answer (t + h, params_ptr)) / exact_answer (t + h, params_ptr))<< "   "
//                              << fabs ((y_rk45[0] - exact_answer (t + h, params_ptr)) / exact_answer (t + h, params_ptr))<< endl;
//
//            }
//      }

    


      // print relative errors to output file
//      out2 << scientific << setprecision (20)
//        << log10 (h) << "   "
//        << log10 (fabs ((y_euler[0] - exact_answer (1. + h, params_ptr)) / exact_answer (1 + h, params_ptr))) << "   "
//        << log10 (fabs ((y_rk2[0] - exact_answer (1. + h, params_ptr)) / exact_answer (1 + h, params_ptr))) << "   "
//        << log10 (fabs ((y_rk4[0] - exact_answer (1. + h, params_ptr)) / exact_answer (1 + h, params_ptr)))<< endl;
//              << h << "   "
//              << fabs ((y_euler[0] - exact_answer (2. + h, params_ptr)) / exact_answer (2. + h, params_ptr)) << "   "
//              << fabs ((y_rk2[0] - exact_answer (2. + h, params_ptr)) / exact_answer (2. + h, params_ptr)) << "   "
//              << fabs ((y_rk4[0] - exact_answer (2. + h, params_ptr)) / exact_answer (2. + h, params_ptr))<< endl;
                    if(fabs(t - 2.) < 0.00000001)
                    {
                                out2<< log10 (h) << "   "
                                << log10 (fabs ((y_euler[0] - exact_answer (t + h, params_ptr)) / exact_answer (t + h, params_ptr))) << "   "
                                << log10 (fabs ((y_rk2[0] - exact_answer (t + h, params_ptr)) / exact_answer (t + h, params_ptr))) << "   "
                                << log10 (fabs ((y_rk4[0] - exact_answer (t + h, params_ptr)) / exact_answer (t + h, params_ptr)))<< "   "
                        <<log10 (fabs ((y_rk45[0] - exact_answer (t + h, params_ptr)) / exact_answer (t + h, params_ptr)))<<endl;
                    }

      
    }
        h /= 2;        // reduce mesh by 2
  }

    

  cout << "data stored in diffeq_test1.dat\n";
  out1.close ();			// close the output file 

  return (0);			// successful completion 
}


//************************** rhs ************************************
//
//  * This is the function defining the right hand side of the diffeq
//
//*********************************************************************
double
rhs (double t, double y[], int i, void *params_ptr)
{
  double a = ((f_parameters *) params_ptr)->alpha;
//   double b = ((f_parameters *) params_ptr)->beta;

  if (i == 0)
  {
    return (-a *  t * y[0]);
  }

  return (1);			// something's wrong if we get here 
}

//********************** exact_answer **************************
//
//  * This is the exact answer for y(t)
//
//******************************************************************
double
exact_answer (double t, void *params_ptr)
{
  // recover a and b from the void pointer params_ptr
  double a = ((f_parameters *) params_ptr)->alpha;
  double b = ((f_parameters *) params_ptr)->beta;

  return (b * exp (-a * t * t / 2.));
}
