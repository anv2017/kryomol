/*****************************************************************************************
                            handler.h  -  description
                             -------------------
This file is part of the KryoMol project.
For more information, see <http://kryomol.sourceforge.io/>

This program is free software; you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation version 2 of the License.
******************************************************************************************/

#ifndef HANDLER_H
#define HANDLER_H

#include <vector>
#include "coordinate.h"
#include "export.h"

namespace kryomol
{
  class MoleculeHandlerPrivate;

  /** \brief A class for translation and trackball rotation of molecules

  This class use the trackball quaterntion methodology
  A MoleculeHandler object is bulit for every molecule loaded in
  order to keep the same scale for all conformers/frames, however, rotation matrixes
  are handled indepently for each of them*/
  class KRYOMOL_API MoleculeHandler
  {
    public:
      MoleculeHandler();
      MoleculeHandler& operator= ( const MoleculeHandler& );
      MoleculeHandler ( const MoleculeHandler& );
      ~MoleculeHandler();
      float& Scale();
      const float& Scale() const;
      const D2Array<float>& Matrix ( size_t frame ) const;
      D2Array<float>& Matrix ( size_t frame );
      void Rotate ( float rotangle,const Coordinate& rotvector,size_t frame );
      void ApplyTransformation ( size_t frame );
      void SetRotationCenter ( const Coordinate& c,size_t frame );
      const Coordinate& RotationCenter ( size_t frame ) const;
      void SetNFrames ( size_t nframes );
      size_t NFrames() const;
      bool IsBlocked() const;
      void SetBlocked(bool b);
      void Reset();
    private:
    void RotateFrame( float rotangle,const Coordinate& rotvector,size_t frame );
    private:
      MoleculeHandlerPrivate* _d;
  };

}
#endif
