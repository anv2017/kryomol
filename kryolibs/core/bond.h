/*****************************************************************************************
                            bond.h  -  description
                             -------------------
This file is part of the KryoMol project.
For more information, see <http://kryomol.sourceforge.io/>

This program is free software; you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation version 2 of the License.
******************************************************************************************/
#ifndef BOND_H
#define BOND_H

#include <cstddef>
#include <vector>

#include "coreexport.h"

namespace kryomol
{
/** @brief representation of a chemical bond

This structure represents a chemical bond as the indexes of the atoms being connected and the bond order between  them, each bond is stored internally as a coupled of indexes i and j being i < j. Therefore a the 1-4 bond can be built as Bond(1,4) or Bond(4,1). For bond comparation algorithms only the indexes and not the order matters.
*/
class KRYOMOLCORE_API Bond
{
  public:
    /** The order of the bond*/
    enum order { SINGLE=1, AROMATIC, DOUBLE, TRIPLE };
    /** Construct a bond between atoms @param i and @param h with order @param e*/
    Bond ( size_t i,  size_t j, order e= SINGLE );
    /** @return the order or the bond*/
    order& Order() { return m_order; }
    /** @return the order of the bond*/
    const order& Order() const { return m_order; }
    /** @return index of atom i*/
    const size_t& I() const { return m_i; }
    /** @return index of atom j*/
    const size_t& J() const { return m_j; }
    /**  @return true if indexes and only indexes of atoms are the same*/
    friend bool operator == ( const Bond& a,const Bond& b ) { return ( ( a.m_i == b.m_i ) && ( a.m_j == b.m_j ) ); }
  private:
    size_t m_i;
    size_t m_j;
    order m_order;
};

}

#endif

