/*****************************************************************************************
                            irspectrum.cpp  -  description
                             -------------------
This file is part of the KryoMol project.
For more information, see <http://kryomol.sourceforge.io/>

This program is free software; you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation version 2 of the License.
******************************************************************************************/

#include <QApplication>
#include <QClipboard>

#include <algorithm>
#include <fstream>
#include <iostream>
#include "irspectrum.h"
#include "qplotspectrum.h"

#ifndef M_PI
# define M_PI       3.14159265358979323846  // matches value in gcc v2 math.h
#endif

const std::string jcampversion="4.24";

IRSpectrum::IRSpectrum() : m_spectrumtype(QPlotSpectrum::IR)
{
  m_title= "gaussian spectrum";
  m_linewidth = 10;
  m_npoints = 8192;
  m_shift = 0.0;

}

IRSpectrum::~IRSpectrum()
{}

void IRSpectrum::SetFrequencies(const std::vector< std::vector<Frequency> >& v,double scale/*=1.0*/)
{
  m_frequencysets=v;
  size_t nsets=v.size();
  m_sinusoidsets.resize(nsets);


  for(size_t i=0;i<nsets;++i)
  {
      auto& fset=m_frequencysets[i];
      m_sinusoidsets[i].resize(fset.size());
      for(auto& f: fset)
      {
          f.x*=scale;
      }
  }

  m_sinusoidsets.clear();
  m_sinusoidsets.resize(v.size());

  for(size_t i=0;i<nsets;++i)
  {
      const auto& fset=m_frequencysets[i];
      auto& sset=m_sinusoidsets[i];
      for(size_t j=0;j<fset.size();++j)
      {
          auto& s=sset[j];
          float amplitude=std::numeric_limits<float>::max();
          switch( m_spectrumtype)
          {
          case QPlotSpectrum::VCD:
              amplitude=fset[j].z;
              break;
          case QPlotSpectrum::IR:
              amplitude=fset[j].y;
              break;
          default:
              throw QString("Not handled QPlotSpectrum type");
          }
          s.m_zamplitude=s.m_flamplitude=amplitude;
          s.m_fldecay=m_linewidth;
          s.m_flphase=0;
          s.m_flt1=0;
          s.m_multiplicity=1;
          s.InitData(m_npoints);
      }
  }

}


void IRSpectrum::SetType( QPlotSpectrum::SpectrumType type )
{
  m_spectrumtype=type;
  SetFrequencies(m_frequencysets);
}


void IRSpectrum::CalculateSpectrum()
{
    m_data.clear();
    m_data.resize(m_frequencysets.size());

    for(const auto& ss : m_sinusoidsets)
    {
        for(auto& d : m_data)
        {
            d.grow(m_npoints);
            float max=m_max;
            float min=m_min;
            float width= (max -min);

            float fdelta_v=width/(d.size()-1);
            for(size_t k=0;k<m_data.size();++k)
            {
                float m=min+k*fdelta_v;
                d[k]=GetIntensityAt(m,ss);

            }
        }
    }


}

float IRSpectrum::GetIntensityAt(float lambda,const std::vector<Sinusoid>& sv)
{

  float sumsignals=0;
  //http://en.wikipedia.org/wiki/Cauchy_distribution
  for( const auto& s: sv )
  {
    float gamma=s.m_fldecay/2;
    float sigma2=gamma*gamma;
    float diff=lambda-s.m_flfrequency;
    //std::cout << "Amplitude=" << it->m_flamplitude << std::endl;
    float x=s.m_flamplitude*gamma/( diff*diff+sigma2);

    sumsignals+=x;
  }
  sumsignals/=M_PI;
  return sumsignals;
}

bool IRSpectrum::WriteJCampDX()
{
  std::ofstream output(m_file.c_str());

  if(!output)
    return false;

  output <<  "##TITLE= " << m_title.c_str() << std::endl;
  output <<  "##JCAMP-DX= " << jcampversion.c_str() << std::endl;

  if (m_spectrumtype == QPlotSpectrum::IR)
      output << "##DATA TYPE= " << "IR SPECTRUM" << std::endl;
  if (m_spectrumtype == QPlotSpectrum::VCD)
      output << "##DATA TYPE= " << "VCD SPECTRUM" << std::endl;

  output << "##NBLOCKS= " << m_frequencysets.size() << std::endl;

  //loop over the number of conformations
  for(size_t i=0;i<m_data.size();++i)
  {

      const fidarray& d=m_data.at(i);
      output << "BLOCK_ID= " << (i+1) << std::endl;
      output << "##XUNITS= 1/CM" << std::endl;
      qWarning() << "take care here with the units, revise" << Qt::endl;
      output << "##YUNITS= ABSORBANCE" << std::endl;
      output << "##XFACTOR= " << 1.0f << std::endl;
      output << "##YFACTOR= " << 0.01f << std::endl;  //0.01 f should be appropiate
      output << "##FIRSTX= " <<  m_min << std::endl;
      output << "##LASTX= " << m_max << std::endl;
      output << "##NPOINTS= " << m_npoints << std::endl;

      //Put now the data
      output << "##XYDATA= (X++(Y..Y))" << std::endl;
      float increment=(m_max-m_min)/(m_npoints-1);
      for (size_t  i=0; i< d.size();i++)
      {
          output << m_min+i*increment << "   " << d[m_npoints-1-i].real() << std::endl;

      }

      output << "##END= " << std::endl;
  }

  return true;
}


void IRSpectrum::CopyData()
{
    QClipboard* c=QApplication::clipboard();
    if ( c == NULL )
    {
        std::cerr << "Could not open clipboard";
    }

    std::stringstream s;

    float width= (m_max - m_min);

    //Get the curve data size
    size_t n=m_expdata.size();
    float fdelta_v=width/(n-1);
    s << "Exp" << "\t" << "Total" << "\t";
    for(size_t i=0;i<n;++i)
    {
        s << "Conf." << (i+1) << "\t";
    }
    for(size_t k=0;k<n;++k)
    {
        float m=m_min+k*fdelta_v;
        s << m << "\t";
        for(const auto& d : m_data )
        {
            s << d[k].real() <<  "\t";
        }
        s << std::endl;;

    }

    c->setText ( s.str().c_str() );

    return;
}


void IRSpectrum::SetFile(const char* file)
{
  m_file=file;
}


void IRSpectrum::SetLineWidth(float lw)
{   
    for(auto& sset : m_sinusoidsets )
    {
        for(auto& s: sset)
        {
            s.m_fldecay=lw;
        }
    }
}


void IRSpectrum::SetShift(float shift)
{
    std::vector<Sinusoid>::iterator it;
    std::vector<Frequency>::iterator vit;

    for(size_t idx=0;idx<m_sinusoidsets.size();++idx )
    {
        std::vector<Sinusoid>& set = m_sinusoidsets[idx];
        for (size_t jdx=0;jdx<set.size();++jdx)
        {
            set[jdx].m_flfrequency=m_frequencysets[idx][jdx].x+shift;
        }
    }

}


void IRSpectrum::SetAutoLimits()
{
  //Get the limits of the spectrum
  ;
  float min=1e10f;
  float max=-1e10f;

  for(auto& fset : m_frequencysets )
  {
        for(auto& f : fset )
        {
            if(f.x >max)
                max=f.x;
            if(f.x < min)
                min=f.x;
        }
  }

  min*=0.9f;
  max*=1.1f;
  if (min < 100 ) min=100; //No very low frequencies by default
  m_min=min;
  m_max=max;
}
