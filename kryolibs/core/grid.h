/*****************************************************************************************
                            grid.h  -  description
                             -------------------
This file is part of the KryoMol project.
For more information, see <http://kryomol.sourceforge.io/>

This program is free software; you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation version 2 of the License.
******************************************************************************************/

#ifndef GRID_H
#define GRID_H

#include "coreexport.h"
#include "vector"

namespace kryomol
{

class KRYOMOLCORE_API Grid
{
    public:
        Grid() {}
        Grid(float x, float y, float z, float step);
        Grid(float x, float y, float z, float step, float l);
        ~Grid() {}
        int Nx() {return m_nx;}
        int Ny() {return m_ny;}
        int Nz() {return m_nz;}
        int Nl() {return m_nl;}
        float X() {return m_x;}
        float Y() {return m_y;}
        float Z() {return m_z;}
        float L() {return m_l;}
        float Step() {return m_step;}
        void SetStep(float step);
        void CalculateSubgrid(const std::vector<float>& alpha);

    private:
        int m_nx;
        int m_ny;
        int m_nz;
        int m_nl;
        float m_x;
        float m_y;
        float m_z;
        float m_l;
        float m_ox;
        float m_oy;
        float m_oz;
        float m_ol;
        float m_step;
};

}

#endif // GRID_H
