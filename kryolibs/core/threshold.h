/*****************************************************************************************
                            threshold.h  -  description
                             -------------------
This file is part of the KryoMol project.
For more information, see <http://kryomol.sourceforge.io/>

This program is free software; you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation version 2 of the License.
******************************************************************************************/

#ifndef THRESHOLD_H
#define THRESHOLD_H

/**
A little structure to store threshold values
*/
//use this structure to store Threshold values fron QM computations
struct Threshold
{
  Threshold() {}
  Threshold(double mf, double rmsf, double maxd, double rmsd) : maxforce(mf), rmsforce(rmsf), maxdisplacement(maxd), rmsdisplacement(rmsd) {}
  double maxforce;
  double rmsforce;
  double maxdisplacement;
  double rmsdisplacement;
};

#endif
