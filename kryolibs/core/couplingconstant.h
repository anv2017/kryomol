/*****************************************************************************************
                            couplingconstant.h  -  description
                             -------------------
This file is part of the KryoMol project.
For more information, see <http://kryomol.sourceforge.io/>

This program is free software; you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation version 2 of the License.
******************************************************************************************/

#ifndef COUPLINGCONSTANT_H
#define COUPLINGCONSTANT_H

#include "coreexport.h"
#include "mathtools.h"

#include <vector>

/** */
namespace kryomol
{
  
  /**
  @brief A general class for storing experimental errors
   
  A class for storing different scalars with existence test*/
  class KRYOMOLCORE_API Error
  {
  public:
    /** Build error with initial value 0*/
    Error() :  m_value(0) {}
    /** */
    ~Error() { delete m_value; }  
    double Value() const { return *m_value; }

    /** Build error object with value v*/
    explicit Error(double v)
    {
      m_value =new double(v);
    }

    /** copy constructor*/
    Error(const Error& e)
    {
        if (  e )
        {
          m_value = new double(e.Value());

        }
        else m_value=0;
        //   if (!m_value && !e)  //Do nothing     
    }

    /** assignment operator*/
    Error& operator = (const Error& e)
    {
      if ( &e != this)
      {
        if ( m_value && e )
        {
          *m_value=e.Value();

        }
        if ( m_value && !e)
        {
          delete m_value;
          m_value=0;
        }
        if ( !m_value && e)
        {
          m_value = new double(e.Value());

        }
        //   if (!m_value && !e)  //Do nothing

      }
      return *this;
    }

     /** asignment operator*/
    Error& operator = (double v)
    {
      {
        if ( m_value ) *m_value=v;
        else
          m_value = new double(v);
      }
     return *this;
    }

  private:
    typedef void (Error::*bool_type) () const;
    void helperfunction() const {}
  private:
    //static unity m_unity;
    double* m_value;
public:
      /** @return true if error has been defined*/
    operator bool_type () const
    {
       if ( m_value != 0 ) return &Error::helperfunction;
       else return 0;
    }

};





/** @brief representation of a coupling constant between nuclei

This class allows management of experimental and computed coupling constant values between nuclei and in general pairwise properties*/
class  KRYOMOLCORE_API CouplingConstant
  {
    public:
    /** Build an undefined value coupling constant  between atoms i and j*/
    CouplingConstant(size_t i, size_t j);
     /** Build a coupling constant between atoms i and j with experimental value ev*/
    CouplingConstant(size_t i, size_t j, double ev);
    /** Build a coupling constant between atoms i and j with experimental value ev and experimental standard deviation error err*/
    CouplingConstant(size_t i, size_t j, double ev,double err);
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
    /** @return index of second coupled nucleus*/
    const size_t& J() const { return _second; }
    /** @return index of first coupled nucleus*/
    size_t& I() { return _first; }
    /** @return index of second coupled nucleus*/
    size_t& J() { return _second; }
    /** Standard deviation error*/
    const kryomol::Error& StandardDeviation() const { return _standarddeviation; }
    /** set the standard deviation error*/
    void SetStandardDeviation(double err) { _standarddeviation=err; }
    const std::vector<kryomol::Error>& ComputedStandardDeviations() const { return _computederrors; }
    std::vector<kryomol::Error>& ComputedStandardDeviations() { return _computederrors; }
    /** Get the conformer averaged standard deviation*/
    const kryomol::Error& AveragedStandardDeviation() const { return _averagedstandarddeviation; }
    /** Set the conformationally averaged standard deviation*/
    void SetAveragedStandardDeviation(double err) { _averagedstandarddeviation=err; }
    /** Get a distribution of values for each conformer*/
    const std::vector< kryomol::D1Array<double> >& Distribution() const { return _distribution; }
    /** Get a distribution of values for each conformer*/
    std::vector< kryomol::D1Array<double> >& Distribution() { return _distribution; }
    /** get a conformationally averaged distribution of values*/
    const kryomol::D1Array<double>& AveragedDistribution() const { return _averageddistribution; }
    /** get a conformationally averaged distribution of values*/
    kryomol::D1Array<double>& AveragedDistribution() { return _averageddistribution; }


private:
    size_t _first;
    size_t _second;
    std::vector<double> _computedvalues;
    double _averagevalue;
    double _experimentalvalue;
    bool _hasexperimentalvalue;
    kryomol::Error _standarddeviation;
    kryomol::Error _averagedstandarddeviation;
    std::vector< kryomol::Error> _computederrors;
    std::vector< D1Array<double> > _distribution;
    D1Array<double> _averageddistribution;
  };   
    /** @return true if indexes I and J are the same*/
    KRYOMOLCORE_API bool operator==(const CouplingConstant &, const  CouplingConstant &);
    /** print frame coordinates on stream s*/
    KRYOMOLCORE_API std::ostream& operator << ( std::ostream& s,const CouplingConstant& c );
}

#endif
