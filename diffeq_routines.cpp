//  file: diffeq_routines.cpp  
//
//  Routines for Euler's and 4th order Runge-Kutta diff. eq. routines 
//                                                                    
//  Programmer:  Dick Furnstahl  furnstahl.1@osu.edu
//                                                                     
//  Revision history:                                                  
//   30-Jan-2004 --- original version, translated from diffeq_routines.c,
//                    which was based on rk4.c from             
//                    "Computational Physics" by Landau and Paez       
//   30-Jan-2005 --- comments improved and function names changed
//   22-Jan-2006 --- made i local to loops
//                                                                     
//   * Based originally on the discussion of differential equations 
//      in Chap. 9 of "Computational Physics" by Landau and Paez
//   * See the 780.20 Session 6 background notes for formulas 
//   * As a convention (advocated in "Practical C++"), we'll append
//      "_ptr" to all pointers.
//
//  To do:
//
//************************************************************************

// include files
#include <iostream>		// note that .h is omitted
#include <iomanip>		// note that .h is omitted
#include <fstream>		// note that .h is omitted
#include <cmath>
#include "diffeq_routines.h"	// diffeq routine prototypes 

const int NMAX=5;

//************************************************************************ 
//  
//   Euler's Algorithm Differential Equation Solver
//
// This routine takes all of the y's one step, from t to t+h.  
//  The original values of y[0], y[1], etc. are lost.
//
// inputs:
//   N --- number of y(t)'s
//   t --- independent variable
//   y[] --- vector of y(t)'s
//   h --- step size
//   f --- function for the right hand sides
//   *params_ptr --- pointer to parameters for rhs function f 
//
// outputs:
//   y[] --- predictions for the values of y(t+h)  
//
//
//************************************************************************
int
euler (const int N, double t, double y[], double h,
       double (*f) (double t, double y[], int i, void *params_ptr),
       void *params_ptr)
{
  for (int i = 0; i < N; i++)
  {
    y[i] += h * f (t, y, i, params_ptr);   // Eq.(6.45) in Session 6 notes 
  }

  return (0);			// successful completion 
}


//************************************************************************
//
//   2nd Order Runge-Kutta Differential Equation Solver
//
// This routine takes all of the y's one step, from t to t+h.
//  The original values of y[0], y[1], etc. are lost.
//
// inputs:
//   N --- number of y(t)'s
//   t --- independent variable
//   y[] --- vector of y(t)'s
//   h --- step size
//   f --- function for the right hand sides
//   *params_ptr --- pointer to parameters for rhs function f
//
// outputs:
//   y[] --- predictions for the values of y(t+h)
//
//
//***********************************************************************
int
runge2 (const int N, double t, double y[], double h,
    double (*f) (double t, double y[], int i, void *params_ptr),
    void *params_ptr)
{
  double y1[NMAX];    // intermediate y values
  double k1[NMAX], k2[NMAX]; // Runge-Kutta notation

  for (int i = 0; i < N; i++)
  {
    k1[i] = h * f (t, y, i, params_ptr);
    y1[i] = y[i] + k1[i] / 2.;    // argument for k2
  }

  for (int i = 0; i < N; i++)
  {
    k2[i] = h * f (t + h / 2., y1, i, params_ptr);
  }


  for (int i = 0; i < N; i++)
  {
    y[i] += k2[i];
  }

  return (0);            // successful completion
}


//************************************************************************ 
//  
//   4th Order Runge-Kutta Differential Equation Solver
//
// This routine takes all of the y's one step, from t to t+h.  
//  The original values of y[0], y[1], etc. are lost.
//
// inputs:
//   N --- number of y(t)'s
//   t --- independent variable
//   y[] --- vector of y(t)'s
//   h --- step size
//   f --- function for the right hand sides
//   *params_ptr --- pointer to parameters for rhs function f 
//
// outputs:
//   y[] --- predictions for the values of y(t+h)  
//
// Notes:
//   * The algorithm is from Eqs. (6.47)-(6.48) in the Session 6 notes.
//
//***********************************************************************
int
runge4 (const int N, double t, double y[], double h,
	double (*f) (double t, double y[], int i, void *params_ptr),
	void *params_ptr)
{
  double y1[NMAX], y2[NMAX], y3[NMAX];	// intermediate y values 
  double k1[NMAX], k2[NMAX], k3[NMAX], k4[NMAX]; // Runge-Kutta notation  

  for (int i = 0; i < N; i++)
  {
    k1[i] = h * f (t, y, i, params_ptr);
    y1[i] = y[i] + k1[i] / 2.;	// argument for k2 
  }

  for (int i = 0; i < N; i++)
  {
    k2[i] = h * f (t + h / 2., y1, i, params_ptr);
    y2[i] = y[i] + k2[i] / 2.;	// argument for k3 
  }

  for (int i = 0; i < N; i++)
  {
    k3[i] = h * f (t + h / 2., y2, i, params_ptr);
    y3[i] = y[i] + k3[i];	// argument for k4 
  }

  for (int i = 0; i < N; i++)
  {
    k4[i] = h * f (t + h, y3, i, params_ptr);
  }

  for (int i = 0; i < N; i++)
  {
    y[i] += (k1[i] + 2. * k2[i] + 2. * k3[i] + k4[i]) / 6.0;
  }

  return (0);			// successful completion 
}

//************************************************************************
//
//   Runge-Kutta-Fehlberg Differential Equation Solver
//
// This routine takes all of the y's one step, from t to t+h.
//  The original values of y[0], y[1], etc. are lost.
//
// inputs:
//   N --- number of y(t)'s
//   t --- independent variable
//   y[] --- vector of y(t)'s
//   h --- step size
//   f --- function for the right hand sides
//   *params_ptr --- pointer to parameters for rhs function f
//
// outputs:
//   y[] --- predictions for the values of y(t+h)
//
// Notes:
//   * The algorithm is from Eqs. (6.47)-(6.48) in the Session 6 notes.
//
//***********************************************************************
int
runge45 (const int N, double t, double y[], double h,
    double (*f) (double t, double y[], int i, void *params_ptr),
    void *params_ptr)
{
    double y1[NMAX], y2[NMAX], y3[NMAX], y4[NMAX], y5[NMAX], err[NMAX];   // intermediate y values
  double k1[NMAX], k2[NMAX], k3[NMAX], k4[NMAX],k5[NMAX],k6[NMAX]; // Runge-Kutta notation
    double s, hmin,hmax,tol;
    hmin = h /64.;
    hmax = h *64.;
    tol = 10E-4;
    err[0] = 1.;
  while(err[0]>tol)
  {

  for (int i = 0; i < N; i++)
  {
    k1[i] = h * f (t, y, i, params_ptr);
    y1[i] = y[i] + k1[i]/4.;    // argument for k2
  }

  for (int i = 0; i < N; i++)
  {
    k2[i] = h * f (t + h/4., y1, i, params_ptr);
    y2[i] = y[i] + 3.* k1[i] / 32. + 9.*k2[i]/32.;    // argument for k3
  }

  for (int i = 0; i < N; i++)
  {
    k3[i] = h * f (t + 3.*h / 8., y2, i, params_ptr);
    y3[i] = y[i] + 1932.*k1[i]/2197.-7200.*k2[i]/2197.+ 7296*k3[i]/2197.;    // argument for k4
  }
    for (int i = 0; i < N; i++)
    {
      k4[i] = h * f (t + 12.*h / 13., y3, i, params_ptr);
      y4[i] = y[i] + 439.*k1[i]/216.-8.*k2[i]+ 3680*k3[i]/513.-845.*k4[i]/4104;    // argument for k5
    }
    for (int i = 0; i < N; i++)
    {
      k5[i] = h * f (t + h, y4, i, params_ptr);
      y5[i] = y[i] - 8.*k1[i]/27. + 2.* k2[i]- 3544.*k3[i] / 2565. + 1859*k4[i]/4104 - 11.*k5[i]/40.;    // argument for k6
    }
  for (int i = 0; i < N; i++)
  {
    k6[i] = h * f (t + h/2, y5, i, params_ptr);
  }
    for (int j = 0; j < N; j++)
    {
      err[j] = fabs(k1[j]/360 - 128. * k3[j]/4275 - 2197. * k4[j]/75240. + k5[j]/50+2*k6[j]/55);

        if(err[j] < tol|| h <= 2*hmin)
        {
            for (int i = 0; i < N; i++)
            {
//                y[i] += (25.*k1[i]/216 + 1408. * k3[i]/2565 + 2197. * k4[i]/4104. - k5[i]/5);
        y[i] += (16.*k1[i]/135 + 6656. * k3[i]/12825. + 28561. * k4[i]/56430. - 9. * k5[i]/50 +2.*k6[i] /55.);
            }

            return(0);// successful completion


        }else{
            s =0.84*pow(tol*h/err[j],0.25);

        if (s<0.75 and h > 2*hmin){
            h /=2.;
        }else if (s> 1.5 and 2*h < hmax)
                {
                    h*=2.;
                }
             }
    
    }

  }
  return (1);

}
