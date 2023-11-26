/*****************************************************************************************
                            density.h  -  description
                             -------------------
This file is part of the KryoMol project.
For more information, see <http://kryomol.sourceforge.io/>

This program is free software; you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation version 2 of the License.
******************************************************************************************/
#ifndef DENSITY_H
#define DENSITY_H

#include <QtOpenGL/QGLWidget>
#include <vector>

#include "coreexport.h"
#include "coordinate.h"
#include "mathtools.h"
#include "orbitalarray.h"
#include "renderdensity.h"
#include "renderexport.h"
#include "coordinate.h"
//#include "frame.h"

namespace kryomol
{
/** @brief representation of an electronic density

This structure represents an electronic density of a molecule.
*/

static const GLint cubeEdgeFlags[256] =
{
        0x000, 0x109, 0x203, 0x30a, 0x406, 0x50f, 0x605, 0x70c, 0x80c, 0x905, 0xa0f, 0xb06, 0xc0a, 0xd03, 0xe09, 0xf00,
        0x190, 0x099, 0x393, 0x29a, 0x596, 0x49f, 0x795, 0x69c, 0x99c, 0x895, 0xb9f, 0xa96, 0xd9a, 0xc93, 0xf99, 0xe90,
        0x230, 0x339, 0x033, 0x13a, 0x636, 0x73f, 0x435, 0x53c, 0xa3c, 0xb35, 0x83f, 0x936, 0xe3a, 0xf33, 0xc39, 0xd30,
        0x3a0, 0x2a9, 0x1a3, 0x0aa, 0x7a6, 0x6af, 0x5a5, 0x4ac, 0xbac, 0xaa5, 0x9af, 0x8a6, 0xfaa, 0xea3, 0xda9, 0xca0,
        0x460, 0x569, 0x663, 0x76a, 0x066, 0x16f, 0x265, 0x36c, 0xc6c, 0xd65, 0xe6f, 0xf66, 0x86a, 0x963, 0xa69, 0xb60,
        0x5f0, 0x4f9, 0x7f3, 0x6fa, 0x1f6, 0x0ff, 0x3f5, 0x2fc, 0xdfc, 0xcf5, 0xfff, 0xef6, 0x9fa, 0x8f3, 0xbf9, 0xaf0,
        0x650, 0x759, 0x453, 0x55a, 0x256, 0x35f, 0x055, 0x15c, 0xe5c, 0xf55, 0xc5f, 0xd56, 0xa5a, 0xb53, 0x859, 0x950,
        0x7c0, 0x6c9, 0x5c3, 0x4ca, 0x3c6, 0x2cf, 0x1c5, 0x0cc, 0xfcc, 0xec5, 0xdcf, 0xcc6, 0xbca, 0xac3, 0x9c9, 0x8c0,
        0x8c0, 0x9c9, 0xac3, 0xbca, 0xcc6, 0xdcf, 0xec5, 0xfcc, 0x0cc, 0x1c5, 0x2cf, 0x3c6, 0x4ca, 0x5c3, 0x6c9, 0x7c0,
        0x950, 0x859, 0xb53, 0xa5a, 0xd56, 0xc5f, 0xf55, 0xe5c, 0x15c, 0x055, 0x35f, 0x256, 0x55a, 0x453, 0x759, 0x650,
        0xaf0, 0xbf9, 0x8f3, 0x9fa, 0xef6, 0xfff, 0xcf5, 0xdfc, 0x2fc, 0x3f5, 0x0ff, 0x1f6, 0x6fa, 0x7f3, 0x4f9, 0x5f0,
        0xb60, 0xa69, 0x963, 0x86a, 0xf66, 0xe6f, 0xd65, 0xc6c, 0x36c, 0x265, 0x16f, 0x066, 0x76a, 0x663, 0x569, 0x460,
        0xca0, 0xda9, 0xea3, 0xfaa, 0x8a6, 0x9af, 0xaa5, 0xbac, 0x4ac, 0x5a5, 0x6af, 0x7a6, 0x0aa, 0x1a3, 0x2a9, 0x3a0,
        0xd30, 0xc39, 0xf33, 0xe3a, 0x936, 0x83f, 0xb35, 0xa3c, 0x53c, 0x435, 0x73f, 0x636, 0x13a, 0x033, 0x339, 0x230,
        0xe90, 0xf99, 0xc93, 0xd9a, 0xa96, 0xb9f, 0x895, 0x99c, 0x69c, 0x795, 0x49f, 0x596, 0x29a, 0x393, 0x099, 0x190,
        0xf00, 0xe09, 0xd03, 0xc0a, 0xb06, 0xa0f, 0x905, 0x80c, 0x70c, 0x605, 0x50f, 0x406, 0x30a, 0x203, 0x109, 0x000
};

//edgeConnection lists the index of the endpoint vertices for each of the 12 edges of the cube
static const GLint cubeEdgeConnection[12][2] =
{
        {0,1}, {1,2}, {2,3}, {3,0},
        {4,5}, {5,6}, {6,7}, {7,4},
        {0,4}, {1,5}, {2,6}, {3,7}
};

//edgeDirection lists the direction vector (vertex1-vertex0) for each edge in the cube
static const GLfloat cubeEdgeDirection[12][3] =
{
        {1.0, 0.0, 0.0},{0.0, 0.0, 1.0},{-1.0, 0.0, 0.0},{0.0, 0.0, -1.0},
        {1.0, 0.0, 0.0},{0.0, 0.0, 1.0},{-1.0, 0.0, 0.0},{0.0, 0.0, -1.0},
        {0.0, 1.0, 0.0},{0.0, 1.0, 0.0},{ 0.0, 1.0, 0.0},{0.0, 1.0, 0.0}
};

//vertexPosition lists the positions, relative to vertex0, of each of the 8 vertices of a cube
static const GLfloat cubeVertexPosition[8][3] =
{
        {0.0, 0.0, 0.0},{1.0, 0.0, 0.0},{1.0, 0.0, 1.0},{0.0, 0.0, 1.0},
        {0.0, 1.0, 0.0},{1.0, 1.0, 0.0},{1.0, 1.0, 1.0},{0.0, 1.0, 1.0}
};

class KRYOMOLRENDER_API Density
{
  public:

    enum unity {ANGSTROM, BOHR};
    Density() { m_isovalue=0.0032; m_resolution=0;}
    Density(int nx, int ny, int nz, float dx, float dy, float dz, const Coordinate& origin=Coordinate(0,0,0));
    Density(int nx, int ny, int nz, int nl, float dx, float dy, float dz, float dl, const Coordinate& origin=Coordinate(0,0,0));
    ~Density() {}

    int Nx() const { return m_nx; }
    int Ny() const { return m_ny; }
    int Nz() const { return m_nz; }
    int Nl() const { return m_nl; }
    float Dx() const {return m_dx; }
    float Dy() const {return m_dy; }
    float Dz() const {return m_dz; }
    float Dl() const {return m_dl; }
    float Dl() {return m_dl;}
    const Coordinate& Origin() const { return m_origin; }
    float Isovalue() const { return m_isovalue; }
    int Resolution() const {return m_resolution;}
    const OrbitalArray& DensityMatrix() const { return m_densitymatrix; }
    bool ExistsDensityData();
    const std::vector<RenderDensity>& PositiveRenderDensityVector() const {return m_positiverenderdensityvector;}
    //std::vector<RenderDensity>& PositiveRenderDensityVector() {return m_positiverenderdensityvector;}
    const std::vector<RenderDensity>& NegativeRenderDensityVector() const {return m_negativerenderdensityvector;}
    //std::vector<RenderDensity> NegativeRenderDensityVector() {return m_negativerenderdensityvector;}
    void SetIsovalue(float isovalue) {m_isovalue = isovalue;}
    void SetResolution(int resolution) {m_resolution = resolution;}
    void SetOrigin(const Coordinate& origin) {m_origin = origin;}
    void SetDensityMatrix(const OrbitalArray& density) {m_densitymatrix = density;}
    void RenderDensityData();
    void PositiveRenderDensityData();
    void NegativeRenderDensityData();
    GLfloat GetPosition(const GLfloat &v1, const GLfloat &v2, const GLfloat &isovalue);
    RenderDensity::GLcoordinate GetNormal(const GLfloat &cX, const GLfloat &cY, const GLfloat &cZ, const bool b);
    RenderDensity::GLcoordinate NormalizeVector(RenderDensity::GLcoordinate &vector);

  private:
    int m_nx;
    int m_ny;
    int m_nz;
    int m_nl;
    float m_dx;
    float m_dy;
    float m_dz;
    float m_dl;
    float m_isovalue;
    int m_resolution;
    Coordinate m_origin;
    OrbitalArray m_densitymatrix;
    std::vector<RenderDensity> m_positiverenderdensityvector;
    std::vector<RenderDensity> m_negativerenderdensityvector;

};

}

#endif
