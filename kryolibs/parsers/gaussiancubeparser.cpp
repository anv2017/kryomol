/*****************************************************************************************
                            gaussiancubeparser.cpp  -  description
                             -------------------
This file is part of the KryoMol project.
For more information, see <http://kryomol.sourceforge.io/>

This program is free software; you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation version 2 of the License.
******************************************************************************************/


#include <iostream>
#include <sstream>
#include "gaussiancubeparser.h"
#include "molecule.h"
#include "stringtools.h"
#include "orbitalarray.h"

using namespace kryomol;

GaussianCubeParser::GaussianCubeParser(const char* file) : Parser(file)
{}

GaussianCubeParser::GaussianCubeParser(std::istream* stream) : Parser(stream)
{}


GaussianCubeParser::~GaussianCubeParser()
{}


bool GaussianCubeParser::ParseFile(std::streampos pos)
{
  std::string line;
  m_file->clear();
  m_file->seekg(pos,std::ios::beg);
  std::getline ( *m_file,line );
  std::getline ( *m_file,line );

  std::cout << line << std::endl;

  Molecules()->push_back(Molecule());
  Molecule& molecule=Molecules()->back();

  int n_atoms;
  bool bohr;
  Coordinate origin;
  std::getline ( *m_file,line );
  StringTokenizer token1(line," \t\r");
  if (token1.size() == 4)
  {
      n_atoms = abs(atof(token1.at(0)));
      if (atof(token1.at(0))<0)
          bohr = true;

      origin.x() = atof(token1.at(1));
      origin.y() = atof(token1.at(2));
      origin.z() = atof(token1.at(3));

      if (bohr)
      {
          origin.x() = 0.529177249*origin.x();
          origin.y() = 0.529177249*origin.y();
          origin.z() = 0.529177249*origin.z();
      }
  }
  else return false;

  size_t nx;
  float dx;
  std::getline ( *m_file,line );
  StringTokenizer token2(line," \t\r");
  if (token2.size()== 4)
  {
      nx = abs(atof(token2.at(0)));      
      dx = atof(token2.at(1));
      if (bohr)
          dx = 0.529177249*dx;
  }
  else return false;

  size_t ny;
  float dy;
  std::getline ( *m_file,line );
  StringTokenizer token3(line," \t\r");
  if (token3.size()== 4)
  {
      ny = abs(atof(token3.at(0)));
      dy = atof(token3.at(2));
      if (bohr)
          dy = 0.529177249*dy;
  }
  else return false;

  size_t nz;
  float dz;
  std::getline ( *m_file,line );
  StringTokenizer token4(line," \t\r");
  if (token4.size()== 4)
  {
      nz = abs(atof(token4.at(0)));
      dz = atof(token4.at(3));
      if (bohr)
          dz = 0.529177249*dz;
  }
  else return false;

  ElectronicDensity density(nx,ny,nz,dx,dy,dz,origin);

  molecule.Frames().push_back(Frame(&molecule));
  Frame& frame= molecule.Frames().back();
  int k=0;
  while ((std::getline ( *m_file,line )) && (k < n_atoms))
  {
      StringTokenizer token(line," \t\r");
      if (token.size()== 5)
      {
          Atom atom(atoi(token.at(0).c_str()));
          molecule.Atoms().push_back(atom);

          Coordinate c;
          c.x()=atof(token.at(2));
          c.y()=atof(token.at(3));
          c.z()=atof(token.at(4));

          if (bohr)
          {
              c.x() = 0.529177249*c.x();
              c.y() = 0.529177249*c.y();
              c.z() = 0.529177249*c.z();
          }

          frame.XYZ().push_back(c);

          k++;

      }
      else return false;
  }

  size_t x=0,y=0,z=0;
  OrbitalArray matrix(nx,ny,nz,0);
  while (std::getline ( *m_file,line ))
  {
      StringTokenizer token(line," \t\r");
      for (size_t i=0; i<token.size(); i++)
      {
          matrix(x,y,z)=atof(token.at(i));
          if (z == nz-1)
          {
              if (y == ny-1)
              {
                  z = 0;
                  y = 0;
                  x++;
              }
              else
              {
                  z = 0;
                  y++;
              }
          }
          else
              z++;
      }
  }
  density.SetDensity(matrix);
  frame.SetElectronicDensityData(density);


#ifdef __GNUC__
#warning supressed move to centroid
#endif

  return true;

}

