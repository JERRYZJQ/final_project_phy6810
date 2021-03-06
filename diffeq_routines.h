//  file: diffeq_routines.h
// 
//  Header file for diffeq_routines.cpp
//
//
//  Programmer:  Dick Furnstahl  furnstahl.1@osu.edu
//
//  Revision History:
//    02/07/04 --- original version converted from C version, 
//                  based on rk4.cpp from "Computational
//                  Physics" and derivative_test.cpp
//
//  To do:
//
//************************************************************************

// function prototypes 
extern int euler ( const int N, double t, double y[], double h,
	    double (*f) (double t, double y[], int i, void *params_ptr), 
            void *params_ptr );
extern int runge2 ( const int N, double t, double y[], double h,
        double (*f) (double t, double y[], int i, void *params_ptr),
            void *params_ptr );
extern int runge4 ( const int N, double t, double y[], double h,
	    double (*f) (double t, double y[], int i, void *params_ptr), 
            void *params_ptr );
 
extern int runge45 ( const int N, double t, double y[], double h,
        double (*f) (double t, double y[], int i, void *params_ptr),
            void *params_ptr );
