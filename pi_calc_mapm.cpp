//Compile with: g++ pi_calc_mapm.cpp libmapm.a -s -o piapproximate -lm
//or
//Compile with: g++ pi_calc_mapm.cpp MAPM/libmapm.a -s -o piapproximate -lm

//This program computes Pi using Mikes Arbirary precison library
//download the library from http://www.tc.umn.edu/~ringx004/mapm-main.html
//Machin-like formula is used to calculate Pi (and this is nice 
//because the algorithm uses integers only :) )

#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <fstream>
#include "MAPM/m_apm.h"
#include <ctime>
#include <limits>
using namespace std;

#if !defined (CLOCKS_PER_SEC) && defined (CLK_TCK)
#define CLOCKS_PER_SEC CLK_TCK
#endif

void printMAPM(const char *caption, MAPM m, long decimalplaces);
MAPM series_ArcTan(MAPM denominator, long iterations);

int main( void )
{
    MAPM piApproximation = 0;
    MAPM denom1 = 57, denom2 = 239, denom3 = 682, denom4 = 12943;
    long iterations = 10;
    long decimalplaces = 100;
    time_t tt;
    clock_t tc;
    cout << "Machin-like formula will be used to calculate PiApprox.\n";
    cout<<"Pi/4 = 44ArcTan(1/57) + 7ArcTan(1/239) - 12ArcTan(1/682) + 24ArcTan(12943)\n\n";
    while( true )
    {
        cout << "How many series terms would you like to approximate ArcTan?: ";
        cin >> iterations;
        cout << "How decimal places of precision should be used for PiApprox? : ";
        cin >> decimalplaces;
        if ( isnan(decimalplaces) || isnan(iterations) ){
            cout << "Error!\nBoth inputs were not number\n";
            cin.clear();
            cin.ignore(std::numeric_limits<std::streamsize>::max(),'\n');}
        else          
            break;
    }    
    m_apm_cpp_precision(decimalplaces);   //Sets the MAPM variable type minimal precision     
    tt = time (0);
    tc = clock ();
     
    piApproximation =   44.0*series_ArcTan( denom1, iterations ) +
                       7.0*series_ArcTan( denom2, iterations ) +
                        -12.0*series_ArcTan( denom3, iterations ) +
                        24*series_ArcTan( denom4, iterations );
    piApproximation *= 4.0;
    cout << "\n\n";
    printMAPM("PiApprox is:..  ", piApproximation, decimalplaces);
    cout << "\nCalculation took " << (clock () - tc) / (double) CLOCKS_PER_SEC << " seconds CPU time, elapsed time " << difftime (time (0), tt) << " seconds" << endl;
    cout << "\n\n\n";
    cout << "PI is stored locally in a file (correct to 1 MIllion decimal places)\n";
    cout << "would you like to compute the difference between local PI and PiApprox? (0/1)?: ";
    cin >> iterations;
    if (iterations == 1)
    {
        ifstream filePiDigits("million_pi_digits.txt");
        char *pidigits = new char[10000100];
        if( filePiDigits.is_open() && pidigits != NULL)
        {    
            filePiDigits >> pidigits;
            denom1 = pidigits;
            filePiDigits.close();
            denom1 = denom1 - piApproximation;
            cout << "\n\n";
            printMAPM("Difference between piApprox and local PI is  ", denom1, 10);
        }
        else
        {
            cout <<"File or String allocation error ... program is now exiting\n";
        }
        delete [] pidigits;    
    }
    return(0);
}


void printMAPM(const char *caption, MAPM m, long decimalplaces)
{
	char *mBuf = new char[decimalplaces+10];
 if( mBuf == NULL)
    return;
	m.toString( mBuf, decimalplaces);
	printf("%s %s\n", caption, mBuf);
 delete [] mBuf;
}

MAPM series_ArcTan(MAPM denominator, long iterations)
{
    MAPM sum = 0;
    int i;
    for( i = 0; i < iterations; i++)
    {
        sum += pow(-1, i)/((2*i + 1)*pow(denominator, 2*i + 1));
//      SEE NOTE BELOW   
//      sum = sum.round( decimalplaces + 10);   //Rounding off stops MAPMs growing without bound!
    }
    return(sum); 
}


//NOTE!: Some real life use with the C++ wrapper has revealed a certain
//       tendency for a program to become quite slow after many iterations
//       (like in a for/while loop).  After a little debug, the reason
//       became clear. Remember that multiplication will maintain ALL
//       significant digits :
//
//       20 digit number x 20 digit number =  40 digits
//       40 digit number x 40 digit number =  80 digits
//       80 digit number x 80 digit number = 160 digits
//       etc.
//
//       So after numerous iterations, the number of significant digits
//       was growing without bound. The easy way to fix the problem is
//       to simply *round* your result after a multiply or some other
//       complex operation. For example:
//
//       #define MAX_DIGITS 256
//
//       p1 = (p0 * sin(b1) * exp(1 + u1)) / sqrt(1 + b1);
//       p1 = p1.round(MAX_DIGITS);
//
//       If you 'round' as shown here, your program will likely be
//       nearly as fast as a straight 'C' version.
