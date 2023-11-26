/*****************************************************************************************
                            chemicalshift.h  -  description
                             -------------------
This file is part of the KryoMol project.
For more information, see <http://kryomol.sourceforge.io/>

This program is free software; you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation version 2 of the License.
******************************************************************************************/

#ifndef CHEMICALSHIFT_H
#define CHEMICALSHIFT_H

#include "couplingconstant.h"

namespace kryomol
{
class KRYOMOLCORE_API ChemicalShift {
public:
    /** Build an undefined value coupling constant  between atoms i and j*/
    ChemicalShift(size_t i);
    /** Build a coupling constant between atoms i and j with experimental value ev*/
    ChemicalShift(size_t i,double ev);
    /** Build a coupling constant between atoms i and j with experimental value ev and experimental standard deviation error err*/
    ChemicalShift(size_t i, double ev,double err);
    /** @return the conformationally averaged computed value*/
    double AverageValue() const { return _averagevalue; }
    /** @return the experimental value*/
    double ExperimentalValue() const { return _experimentalvalue; }
    /** sets the experimental value*/
    void SetExperimentalValue(double ev);
    /** true if a experimental value has been defined */
    bool HasExperimentalValue() const { return _hasexperimentalvalue; }
    /** set the computed conformationally averaged averaged value*/
    void SetAverageValue(double v) { _averagevalue=v; }
    /** @return a const vector with computed values for all conformers*/
    const std::vector<double>& ComputedValues() const { return _computedvalues; }
    /** @return a const vector with computed values for all conformers*/
    std::vector<double>& ComputedValues() { return _computedvalues; }
    /** @return index of first coupled nucleus*/
    const size_t& I() const { return _first; }
    /** @return index of first coupled nucleus*/
    size_t& I() { return _first; }
    /** Standard deviation error*/
    const kryomol::Error& StandardDeviation() const { return _standarddeviation; }
    /** set the standard deviation error*/
    void SetStandardDeviation(double err) { _standarddeviation=err; }
    /** return reference to vector with standard deviations for each frame*/
    const std::vector<kryomol::Error>& ComputedStandardDeviations() const { return _computederrors; }
    /**  return const reference to vector with standard deviations for each frame*/
    std::vector<kryomol::Error>& ComputedStandardDeviations() { return _computederrors; }
    /** get a distribution of values*/
    const std::vector< D1Array<double> >& Distribution() const { return _distribution; }
    /** Get a distribution of values*/
    std::vector < D1Array<double> >&  Distribution() { return _distribution; }

private:
    size_t _first;
    std::vector<double> _computedvalues;
    std::vector<kryomol::Error> _computederrors;
    double _averagevalue;
    double _experimentalvalue;
    bool _hasexperimentalvalue;
    kryomol::Error _standarddeviation;
    std::vector< D1Array<double> > _distribution;
};
/** @return true if indexes I and J are the same*/
KRYOMOLCORE_API bool operator==(const ChemicalShift &, const  ChemicalShift &);


}

#endif
