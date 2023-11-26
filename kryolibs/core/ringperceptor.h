/*****************************************************************************************
                            ringperceptor.h  -  description
                             -------------------
This file is part of the KryoMol project.
For more information, see <http://kryomol.sourceforge.io/>

This program is free software; you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation version 2 of the License.
******************************************************************************************/

#ifndef RINGPERCEPTOR_H

#include <vector>
#include "molecule.h"

namespace kryomol
{

  /** @brief Implementation of a ring perception algorithm

  Implementation of the ring perception algorithm of Hansser et al
  J. Chem. Inf. Comput. Sci. 1996, 36, 1146-1152
  */
  class RingPerceptor
  {
      class Edge
      {
        public:
          Edge();
          std::vector<size_t>& Path() { return _path; }
          const std::vector<size_t>& Path() const { return _path;}
          //Edge operator ^ (Edge& b);
          Edge operator+ ( Edge& b );
        private:
          std::vector< size_t > _path;
      };
    public:
      RingPerceptor ( const kryomol::Molecule* molecule );
      RingPerceptor ( const kryomol::Frame* frame );
      ~RingPerceptor();
      /** @return a vector with the ring paths*
          @param size if size==0 return all rings */
      std::vector< std::vector<size_t> >Rings ( size_t size=0 );
    private:
      void Convert();
      void Remove();
      size_t CommonElements ( const Edge& edg1, const Edge& edg2 );
    private:
      const Molecule* m_molecule;
      const Frame* m_frame;
      std::vector< size_t > m_vertex;
      std::vector < Edge > m_edges;
      std::vector<Edge> m_rings;
  };
}
#endif
