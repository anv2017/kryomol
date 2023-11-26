/*****************************************************************************************
                            grid.cpp  -  description
                             -------------------
This file is part of the KryoMol project.
For more information, see <http://kryomol.sourceforge.io/>

This program is free software; you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation version 2 of the License.
******************************************************************************************/

#include <QDebug>
#include <math.h>
#include "grid.h"

using namespace kryomol;

Grid::Grid(float x, float y, float z, float step) : m_x(x), m_y(y), m_z(z), m_ox(x), m_oy(y), m_oz(z), m_step(step)
{
    // Calculate and assure the number of points of the grid are multiple of 4 for using SSE
    m_nx = floor(x/step)+1;
    while (m_nx%4 != 0)
    {
        x += step;
        m_x += step;
        m_nx = floor(x/step)+1;
    }
    m_ny = floor(y/step)+1;
    while (m_ny%4 != 0)
    {
        y += step;
        m_y += step;
        m_ny = floor(y/step)+1;
    }
    m_nz = floor(z/step)+1;
    while (m_nz%4 != 0)
    {
        z += step;
        m_z += step;
        m_nz = floor(z/step)+1;
    }
}

Grid::Grid(float x, float y, float z, float step, float l) : m_x(x), m_y(y), m_z(z), m_l(l), m_ox(x), m_oy(y), m_oz(z), m_ol(l), m_step(step)
{
    // Calculate and assure the number of points of the grid are multiple of 4 for using SSE
    m_nx = floor(x/step)+1;
    while (m_nx%4 != 0)
    {
        x += step;
        m_x += step;
        m_nx = floor(x/step)+1;
    }
    m_ny = floor(y/step)+1;
    while (m_ny%4 != 0)
    {
        y += step;
        m_y += step;
        m_ny = floor(y/step)+1;
    }
    m_nz = floor(z/step)+1;
    while (m_nz%4 != 0)
    {
        z += step;
        m_z += step;
        m_nz = floor(z/step)+1;
    }
    m_nl = floor(l/step)+1;
    while (m_nl%4 != 0)
    {
        l += step;
        m_l += step;
        m_nl = floor(l/step)+1;
    }
}

void Grid::SetStep(float step)
{
    m_step = step;

    m_x = m_ox;
    m_y = m_oy;
    m_z = m_oz;
    m_l = m_ol;

    float x = m_ox;
    float y = m_oy;
    float z = m_oz;
    float l = m_ol;

    // Calculate the new values of the grid and assure the number of points of the grid are multiple of 4 for using SSE
    m_nx = floor(x/step)+1;
    while (m_nx%4 != 0)
    {
        x += step;
        m_x += step;
        m_nx = floor(x/step)+1;
    }
    m_ny = floor(y/step)+1;
    while (m_ny%4 != 0)
    {
        y += step;
        m_y += step;
        m_ny = floor(y/step)+1;
    }
    m_nz = floor(z/step)+1;
    while (m_nz%4 != 0)
    {
        z += step;
        m_z += step;
        m_nz = floor(z/step)+1;
    }
    m_nl = floor(l/step)+1;
    while (m_nl%4 != 0)
    {
        l += step;
        m_l += step;
        m_nl = floor(l/step)+1;
    }
}

void Grid::CalculateSubgrid(const std::vector<float>& valpha)
{
    //For a confidence of 98% we use a subgrid's length = 2.32635/sqrt(2*alpha) <- (gaussian distribution)
    float alpha = valpha.at(0);
    for (size_t i=0;i<valpha.size();++i)
        if (alpha>valpha.at(i))
            alpha = valpha.at(i);
    m_l =  2*2.32635/(sqrt(2.0*alpha));

    //If the value of subgrid is bigger than the maximum side of the grid, limite the subgrid to the maximum value of the grid
    float v_max = std::max(m_x,m_y);
    v_max = std::max(v_max,m_z);
    if (m_l>v_max)
        m_l = v_max;

    // Assure the number of points of the subgrid are multiple of 4 for using SSE
    m_nl = floor(m_l/m_step)+1;
    while (m_nl%4 != 0)
    {
        m_l += m_step;
        m_nl = floor(m_l/m_step)+1;
    }
}
