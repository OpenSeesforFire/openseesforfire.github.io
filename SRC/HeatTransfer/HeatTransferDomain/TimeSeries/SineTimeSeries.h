/* ****************************************************************** **
**    OpenSees - Open System for Earthquake Engineering Simulation    **
**          Pacific Earthquake Engineering Research Center            **
**                                                                    **
**                                                                    **
** (C) Copyright 2001, The Regents of the University of California    **
** All Rights Reserved.                                               **
**                                                                    **
** Commercial use of this program without express permission of the   **
** University of California, Berkeley, is strictly prohibited.  See   **
** file 'COPYRIGHT'  in main directory for information on usage and   **
** redistribution,  and for a DISCLAIMER OF ALL WARRANTIES.           **
**                                                                    **
** Developed by:                                                      **
**   Frank McKenna (fmckenna@ce.berkeley.edu)                         **
**   Gregory L. Fenves (fenves@ce.berkeley.edu)                       **
**   Filip C. Filippou (filippou@ce.berkeley.edu)                     **
**                                                                    **
** Fire & Heat Transfer modules developed by:                         **
**   Yaqiang Jiang (y.jiang@ed.ac.uk)                                 **
**   Asif Usmani (asif.usmani@ed.ac.uk)                               **
**                                                                    **
** ****************************************************************** */

//Adapted by Yaqiang Jiang (y.jiang@ed.ac.uk)

#ifndef SineTimeSeries_h
#define SineTimeSeries_h

// Written: Andreas Schellenberg (andreas.schellenberg@gmx.net)
// Created: 02/04
// Revision: A
//
// Purpose: This file contains the class definition for PulseSeries.
// PulseSeries is a concrete class. A PulseSeries object provides
// a pulse time series. The factor is given by the pseudoTime (t),
// pulse period (T), pulse width (pw) and phase shift (phi),
// and a constant factor provided in the constructor,
// the duration by tStart and tFinal;
//
// What: "@(#) PulseSeries.h, revA"

#include <TimeSeries.h>

class SineTimeSeries : public TimeSeries
{
public:
    // constructors
    SineTimeSeries(int tag, double tStart, double tFinish, double period,
		           double peakAmplitude = 1.0);

    SineTimeSeries();
    
    // destructor
    ~SineTimeSeries();

    TimeSeries* getCopy(void);
    
    // method to get load factor
    double getFactor(double time);
    double getDuration () {return tFinish-tStart;}
    double getPeakFactor () {return peakAmplitude;}
    double getTimeIncr (double pseudoTime) {return tFinish-tStart;}
    
    // methods for output    
    int sendSelf(int commitTag, Channel &theChannel);
    int recvSelf(int commitTag, Channel &theChannel, 
		         FEM_ObjectBroker &theBroker);

    void Print(OPS_Stream &s, int flag =0);    
    
protected:
	
private:
    double tStart;    // start time of time series (sec)
    double tFinish;   // end time of time series (sec)
    double period;    // period of pulse series (sec)
    double peakAmplitude;   // amplitude of pulse series
};

#endif