/*****************************************************************************************
                            pdbtools.h  -  description
                             -------------------
This file is part of the KryoMol project.
For more information, see <http://kryomol.sourceforge.io/>

This program is free software; you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation version 2 of the License.
******************************************************************************************/

#ifndef PDBRESIDUE_H
#define PDBRESIDUE_H

#include "coreexport.h"
#include <string>

namespace kryomol
{
  class KRYOMOLCORE_API PDBResidue
  {
    public :
      PDBResidue(const std::string& name="",const std::string& index="");
      const std::string& Name() const { return _name; }
      const std::string& Index() const { return _index; }
      bool Visible() const { return _visible; }
      void SetVisible(bool b) { _visible=b; }
      friend bool operator == ( const kryomol::PDBResidue& a,const kryomol::PDBResidue& b );
    private:
      std::string _name;
      std::string _index;
      bool _visible;
  };
}

#endif
