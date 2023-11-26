/*****************************************************************************************
                            energy.h  -  description
                             -------------------
This file is part of the KryoMol project.
For more information, see <http://kryomol.sourceforge.io/>

This program is free software; you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation version 2 of the License.
******************************************************************************************/

#ifndef KRYOMOLENERGY_H
#define KRYOMOLENERGY_H

#include "coreexport.h"

/** */
namespace kryomol
{

  /**
  @brief A class for storing energy values
   
  A class for storing different kind of energy values*/
  class KRYOMOLCORE_API Energy
  {
  public:
    /** energy unities*/
    enum unity { KCAL, KJ, HARTREES };
    /** Build energy with initial value 0*/
    Energy() :  m_value(0) {}
    /** */
    ~Energy() { delete m_value; }  
    double Value() const { return *m_value; }
    /** @return the magnitud of the energy value in unities @see unity*/
//     operator double () const { return *m_value; }

    /** Build energy object with value v*/
    explicit Energy(double v)
    {
      m_value =new double(v);
    }
    /** copy constructor*/
    Energy(const Energy& e)
    {
        if (  e )
        {
          m_value = new double(e.Value());

        }
        else m_value=0;
        //   if (!m_value && !e)  //Do nothing      
    }
    /** assignment operator*/
    Energy& operator = (const Energy& e)
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
    Energy& operator = (double v)
    {
      {
        if ( m_value ) *m_value=v;
        else
          m_value = new double(v);
      }
     return *this;
    }
    /** return energy unities setup globally*/
    static unity Unity();
    /** Setup globally energy unities*/
    static void SetUnity(unity u);
  private:
    typedef void (Energy::*bool_type) () const;
    void helperfunction() const {}
  private:
    //static unity m_unity;
    double* m_value;
    static unity m_unity;
public:
      /** @return true if energy has been defined*/
    operator bool_type () const
    {
       if ( m_value != 0 ) return &Energy::helperfunction;
       else return 0;
    }

public:
  /** addition of energies*/
  friend Energy operator + (const Energy& a, const Energy& b)
  {
    if ( !a || !b )
    {
       return Energy();
    }
    else
    {
       return Energy(*a.m_value + *b.m_value );
    }
  } 

  };
  
}


#endif
