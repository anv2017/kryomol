/*****************************************************************************************
                            renderdensity.h  -  description
                             -------------------
This file is part of the KryoMol project.
For more information, see <http://kryomol.sourceforge.io/>

This program is free software; you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation version 2 of the License.
******************************************************************************************/
#ifndef RENDERDENSITY_H
#define RENDERDENSITY_H

#include <QtOpenGL/QGLWidget>
#include <vector>

#include "mathtools.h"
#include "coreexport.h"

namespace kryomol
{
/** @brief data for drawing an electronic density

This structure represents an electronic density of a molecule.
*/

class KRYOMOLCORE_API RenderDensity
{
  public:

    struct GLcoordinate
    {
            GLfloat X;
            GLfloat Y;
            GLfloat Z;
    };

    RenderDensity() { }
    RenderDensity(GLint flagindex, const std::vector<GLcoordinate>& edgevertex, const std::vector<GLcoordinate>& edgenorm);

    GLint FlagIndex() const { return m_flagindex;}
    const std::vector<GLcoordinate>& EdgeVertex() const {return m_edgevertex;}
    std::vector<GLcoordinate>& EdgeVertex() {return m_edgevertex;}
    const std::vector<GLcoordinate>& EdgeNorm() const {return m_edgenorm;}
    std::vector<GLcoordinate>& EdgeNorm() {return m_edgenorm;}
    const std::vector<GLcoordinate>& Color() const {return m_color;}
    std::vector<GLcoordinate>& Color() {return m_color;}

  private:
    GLint m_flagindex;
    std::vector<GLcoordinate> m_edgevertex;
    std::vector<GLcoordinate> m_edgenorm;
    std::vector<GLcoordinate> m_color;

};

}

#endif
