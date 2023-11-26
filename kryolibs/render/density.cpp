/*****************************************************************************************
                            density.cpp  -  description
                             -------------------
This file is part of the KryoMol project.
For more information, see <http://kryomol.sourceforge.io/>

This program is free software; you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation version 2 of the License.
******************************************************************************************/

#include "density.h"
#include "renderdensity.h"

using namespace kryomol;

Density::Density(int nx, int ny, int nz, float dx, float dy, float dz, const Coordinate &origin) : m_nx(nx), m_ny(ny), m_nz(nz), m_dx(dx), m_dy(dy), m_dz(dz), m_origin(origin)
{
    m_isovalue = 0.0032;
    m_resolution = 0;
}

Density::Density(int nx, int ny, int nz, int nl, float dx, float dy, float dz, float dl, const Coordinate &origin) : m_nx(nx), m_ny(ny), m_nz(nz), m_nl(nl), m_dx(dx), m_dy(dy), m_dz(dz), m_dl(dl), m_origin(origin)
{
    m_isovalue = 0.0032;
    m_resolution = 0;
}

void Density::RenderDensityData()
{
    PositiveRenderDensityData();
    NegativeRenderDensityData();
}

void Density::PositiveRenderDensityData()
{
    m_positiverenderdensityvector.erase(m_positiverenderdensityvector.begin(),m_positiverenderdensityvector.end());

    int x=0;
    while (x<m_nx-1)
    {
        GLfloat cX = x*m_dx + m_origin.x();

        int y=0;
        while (y<m_ny-1)
        {
            GLfloat cY = y*m_dy+m_origin.y();

            int z=0;
            while (z<m_nz-1)
            {
                GLfloat cZ = z*m_dz+m_origin.z();

                //Find which vertices of the cube are inside of the surface and which are outside
                GLfloat cubevalues[8] = {m_densitymatrix(x,y,z), m_densitymatrix(x+1,y,z), m_densitymatrix(x+1,y,z+1), m_densitymatrix(x,y,z+1), m_densitymatrix(x,y+1,z), m_densitymatrix(x+1,y+1,z), m_densitymatrix(x+1,y+1,z+1), m_densitymatrix(x,y+1,z+1)};

                GLint flagsindex = 0;
                for (GLint vertex=0; vertex<8; vertex++)
                {
                    if (cubevalues[vertex] <= m_isovalue)
                    {
                        flagsindex |= 1<<vertex;
                    }
                }
                //Find which edges of the cube are intersected by the surface
                GLint edgeflags = cubeEdgeFlags[flagsindex];

                //If the cube is entirely inside or outside of the surface, then there will be no intersections
                if (edgeflags != 0)
                {
                    std::vector<RenderDensity::GLcoordinate> v_edgevertex;
                    std::vector<RenderDensity::GLcoordinate> v_edgenorm;

                    //Find the point of intersection of the surface with each edge
                    //Then find the normal to the surface at those points
                    for(GLint edge = 0; edge < 12; edge++)
                    {
                        RenderDensity::GLcoordinate edgevertex;
                        RenderDensity::GLcoordinate edgenorm;

                        //if there is an intersection on this edge
                        if(edgeflags & (1<<edge))
                        {
                            GLfloat offset = GetPosition(cubevalues[cubeEdgeConnection[edge][0] ],cubevalues[cubeEdgeConnection[edge][1] ], m_isovalue);

                            edgevertex.X = cX + (cubeVertexPosition[ cubeEdgeConnection[edge][0] ][0]  +  offset * cubeEdgeDirection[edge][0])*m_dx;
                            edgevertex.Y = cY + (cubeVertexPosition[ cubeEdgeConnection[edge][0] ][1]  +  offset * cubeEdgeDirection[edge][1])*m_dy;
                            edgevertex.Z = cZ + (cubeVertexPosition[ cubeEdgeConnection[edge][0] ][2]  +  offset * cubeEdgeDirection[edge][2])*m_dz;

                            edgenorm = GetNormal(x, y, z, false);
                        }
                        v_edgevertex.push_back(edgevertex);
                        v_edgenorm.push_back(edgenorm);
                    }
                    RenderDensity renderdensity(flagsindex,v_edgevertex,v_edgenorm);
                    m_positiverenderdensityvector.push_back(renderdensity);
                }
                ++z;//z = z++;
            }
            ++y;//y = y++;
        }
        ++x;//x = x++;
    }

}


void Density::NegativeRenderDensityData()
{
    m_negativerenderdensityvector.erase(m_negativerenderdensityvector.begin(),m_negativerenderdensityvector.end());

    int x=0;
    while (x<m_nx-1)
    {
        GLfloat cX = x*m_dx + m_origin.x();

        int y=0;
        while (y<m_ny-1)
        {
            GLfloat cY = y*m_dy + m_origin.y();

            int z=0;
            while (z<m_nz-1)
            {
                GLfloat cZ = z*m_dz + m_origin.z();

                //Find which vertices of the cube are inside of the surface and which are outside
                GLfloat cubevalues[8] = {-m_densitymatrix(x,y,z), -m_densitymatrix(x+1,y,z), -m_densitymatrix(x+1,y,z+1), -m_densitymatrix(x,y,z+1), -m_densitymatrix(x,y+1,z), -m_densitymatrix(x+1,y+1,z), -m_densitymatrix(x+1,y+1,z+1), -m_densitymatrix(x,y+1,z+1)};

                GLint flagsindex = 0;
                for (GLint vertex=0; vertex<8; vertex++)
                {
                    if (cubevalues[vertex] <= m_isovalue)
                    {
                        flagsindex |= 1<<vertex;
                    }
                }
                //Find which edges of the cube are intersected by the surface
                GLint edgeflags = cubeEdgeFlags[flagsindex];

                //If the cube is entirely inside or outside of the surface, then there will be no intersections
                if (edgeflags != 0)
                {
                    std::vector<RenderDensity::GLcoordinate> v_edgevertex;
                    std::vector<RenderDensity::GLcoordinate> v_edgenorm;

                    //Find the point of intersection of the surface with each edge
                    //Then find the normal to the surface at those points
                    for(GLint edge = 0; edge < 12; edge++)
                    {
                        RenderDensity::GLcoordinate edgevertex;
                        RenderDensity::GLcoordinate edgenorm;

                        //if there is an intersection on this edge
                        if(edgeflags & (1<<edge))
                        {
                            GLfloat offset = GetPosition(cubevalues[cubeEdgeConnection[edge][0] ],cubevalues[cubeEdgeConnection[edge][1] ], m_isovalue);

                            edgevertex.X = cX + (cubeVertexPosition[ cubeEdgeConnection[edge][0] ][0]  +  offset * cubeEdgeDirection[edge][0])*m_dx;
                            edgevertex.Y = cY + (cubeVertexPosition[ cubeEdgeConnection[edge][0] ][1]  +  offset * cubeEdgeDirection[edge][1])*m_dy;
                            edgevertex.Z = cZ + (cubeVertexPosition[ cubeEdgeConnection[edge][0] ][2]  +  offset * cubeEdgeDirection[edge][2])*m_dz;

                            edgenorm = GetNormal(x, y, z,true);
                        }
                        v_edgevertex.push_back(edgevertex);
                        v_edgenorm.push_back(edgenorm);
                    }
                    RenderDensity renderdensity(flagsindex,v_edgevertex,v_edgenorm);
                    m_negativerenderdensityvector.push_back(renderdensity);
                }
                ++z;//z = z++;
            }
            ++y;//y = y++;
        }
        ++x;//x = x++;
    }

}



//Find the approximate point of intersection of the surface between two points with the values v1 and v2
GLfloat Density::GetPosition(const GLfloat &v1, const GLfloat &v2, const GLfloat &isovalue)
{
    GLdouble delta = v2 - v1;
    if(delta == 0.0)
        return 0.5;
    return (isovalue - v1)/delta;
}

//Find the gradient of the scalar field at a point
RenderDensity::GLcoordinate Density::GetNormal(const GLfloat &cX, const GLfloat &cY, const GLfloat &cZ, const bool negative)
{
    RenderDensity::GLcoordinate normal;

    GLfloat X,Y,Z;
    if (m_resolution<0)
    {
        if (cX < (abs(m_resolution)+1))
            X=cX;
        else
            X=cX-(abs(m_resolution)+1);
        if (cY < (abs(m_resolution)+1))
            Y=cY;
        else
            Y=cY-(abs(m_resolution)+1);
        if (cZ < (abs(m_resolution)+1))
            Z=cZ;
        else
            Z=cZ-(abs(m_resolution)+1);

        normal.X = m_densitymatrix(X, cY, cZ) - m_densitymatrix(cX+(abs(m_resolution)+1), cY, cZ);
        normal.Y = m_densitymatrix(cX, Y, cZ) - m_densitymatrix(cX, cY+(abs(m_resolution)+1), cZ);
        normal.Z = m_densitymatrix(cX, cY, Z) - m_densitymatrix(cX, cY, cZ+(abs(m_resolution)+1));
    }
    else
    {
        if (cX < 1)
            X=cX;
        else
            X=cX-1;
        if (cY < 1)
            Y=cY;
        else
            Y=cY-1;
        if (cZ < 1)
            Z=cZ;
        else
            Z=cZ-1;

        normal.X = m_densitymatrix(X, cY, cZ) - m_densitymatrix(cX+1, cY, cZ);
        normal.Y = m_densitymatrix(cX, Y, cZ) - m_densitymatrix(cX, cY+1, cZ);
        normal.Z = m_densitymatrix(cX, cY, Z) - m_densitymatrix(cX, cY, cZ+1);

    }

    if (negative)
    {
        normal.X = -normal.X;
        normal.Y = -normal.Y;
        normal.Z = -normal.Z;
    }

    normal = NormalizeVector(normal);

    return normal;
}


RenderDensity::GLcoordinate Density::NormalizeVector(RenderDensity::GLcoordinate &v)
{
        RenderDensity::GLcoordinate vnormal;
        GLfloat length;

        length = sqrt((v.X*v.X)+(v.Y*v.Y)+(v.Z*v.Z));

        if(length == 0.0)
        {
                vnormal.X = v.X;
                vnormal.Y = v.Y;
                vnormal.Z = v.Z;
        }
        else
        {
                vnormal.X = v.X/length;
                vnormal.Y = v.Y/length;
                vnormal.Z = v.Z/length;
        }

        return vnormal;
}

bool Density::ExistsDensityData()
{
    bool b = false;
    if (!(m_densitymatrix.Empty()))
        b = true;

    return b;
}

