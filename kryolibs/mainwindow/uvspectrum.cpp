/*****************************************************************************************
                            uvspectrum.cpp  -  description
                             -------------------
This file is part of the KryoMol project.
For more information, see <http://kryomol.sourceforge.io/>

This program is free software; you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation version 2 of the License.
******************************************************************************************/

#include <QClipboard>
#include <QApplication>

#include <algorithm>
#include <fstream>
#include <iostream>

#include "uvspectrum.h"
#include "irspectrum.h"
#include "qplotspectrum.h"

#ifndef M_PI
# define M_PI       3.14159265358979323846  // matches value in gcc v2 math.h
#endif

const std::string jcampversion="4.24";
const bool useboltzmannweighting=false;

UVSpectrum::UVSpectrum() : m_spectrumtype(QPlotSpectrum::UV)
{
    m_title= "gaussian spectrum";
    m_linewidth = 0.3;
    m_npoints = 512;
    m_shift = 0.0;
    m_formalism=length;
    m_benantiomer=false;
    m_boltzw=useboltzmannweighting;
    m_populations=nullptr;
}

UVSpectrum::~UVSpectrum()
{}

void UVSpectrum::SetLines(std::vector< std::vector<Spectralline> >& v,double scale/*=1.0*/)
{
    m_linesets=v;
    m_sinusoidsets.resize(m_linesets.size());
    for(size_t idx=0;idx<m_linesets.size();++idx)
    {
        size_t nlines=m_linesets[idx].size();
        m_sinusoidsets[idx].resize(nlines);
    }

    m_weights=std::vector<float>(m_linesets.size(),1.0f/m_linesets.size());

    for (auto& lineset : v )
    {
        for(auto& l : lineset)
        {
            l.x*=scale;
        }

    }

    for(size_t idx=0;idx<m_linesets.size();++idx)
    {
        std::vector<Sinusoid>& sv=m_sinusoidsets[idx];
        std::vector<Spectralline>& sl=m_linesets[idx];
        for(size_t jdx=0;jdx<sv.size();++jdx)
        {
            Sinusoid& svv=sv[jdx];
            const Spectralline& sll=sl[jdx];
            if ( m_spectrumtype == QPlotSpectrum::ECD )
            {
                int enantiomerize=1;
                if ( m_benantiomer ) enantiomerize=-1;

                double y;
                if ( m_formalism == length)
                    y=sll.RotatoryStrengthLength()*enantiomerize; //proportional to rotatory strength
                else
                    y=sll.RotatoryStrengthVelocity()*enantiomerize;

                svv.m_zamplitude=svv.m_flamplitude=y;

            } else if ( m_spectrumtype == QPlotSpectrum::UV )
            {
                svv.m_zamplitude=svv.m_flamplitude=sll.y0;
            }

            svv.m_fldecay=m_linewidth;
            svv.m_flfrequency=sll.x+m_shift;
            svv.m_flphase=0.0f;
            svv.m_flt1=0.0f;
            svv.m_J=0.0f;
            svv.m_multiplicity=1;
            svv.InitData(m_npoints);
        }
    }

}


void UVSpectrum::SetType( QPlotSpectrum::SpectrumType type )
{
    m_spectrumtype=type;
    SetLines(m_linesets);
}


void UVSpectrum::CalculateSpectrum()
{
    m_data.clear();
    m_data.resize(m_linesets.size());
    m_totaldata.clear();
    m_totaldata.grow(m_npoints);
    for(size_t i=0;i<m_npoints;++i)
    {
        m_totaldata[i]=std::complex<float>(0,0);
    }
    for(size_t idx=0;idx<m_sinusoidsets.size();++idx)
    {
        auto& d=m_data[idx];
        d.grow(m_npoints);
        float max=m_max;
        float min=m_min;
        float width=(max-min);
        float fdelta_v=width/(d.size()-1);
        for(size_t k=0;k<d.size();++k)
        {
            float m=min+k*fdelta_v;
            d[k]=GetIntensityAt(m,m_sinusoidsets[idx]);

            if ( this->BoltzmannWeighting() )
            {
                double p=(*m_populations)[idx];
                d[k]*=p;
            } else
            {
                d[k]*=m_weights[idx];
            }
            m_totaldata[k]+=(d[k]);
        }    
    }
}

//We are using the following equation
float UVSpectrum::GetIntensityAt(float lambda,const std::vector<Sinusoid>& sdd)
{
    float eval=1240/lambda;
    float prefactor=1/(22.97*sqrt(2*M_PI));
    float sumsignals=0;
    for( const auto& s : sdd)
    {
        //Get the transition uv value in ev
        float deltae=1240/s.m_flfrequency;
        float sigma=s.m_fldecay/2.3548; //FWHM
        float expfactor=(eval-deltae)/(sigma);
        float cd=(1/sigma)*deltae*s.m_flamplitude*exp(-( expfactor*expfactor )/2 );

        sumsignals+=cd;
    }
    sumsignals*=prefactor;
    return sumsignals;
}

bool UVSpectrum::WriteJCampDX()
{
    std::ofstream output(m_file.c_str());

    if(!output)
    {
        std::cout << "cannot open file:" << m_file.c_str() << std::endl;
        return false;
    }

    output <<  "##TITLE= " << m_title.c_str() << std::endl;
    output <<  "##JCAMP-DX= " << jcampversion.c_str() << std::endl;
    output << "##DATA TYPE= " << "LINK" << std::endl;


    output << "##ORIGIN= " << "KryoMol" << std::endl;

    //include total data
    output << "##BLOCKS= " << m_linesets.size()+1 << std::endl;

    //loop over the number of conformations
    for(size_t i=0;i<m_data.size()+1;++i)
    {
        const fidarray* pd=nullptr;
        if (i < m_data.size() )
        {
            pd=&m_data.at(i);
        }
        else pd=&m_totaldata;

        const fidarray& d=*pd;
        output <<  "##TITLE= " << m_title.c_str() << std::endl;
        output <<  "##JCAMP-DX= " << jcampversion.c_str() << std::endl;
        if (m_spectrumtype == QPlotSpectrum::UV)
            output << "##DATA TYPE= " << "UV SPECTRUM" << std::endl;

        if (m_spectrumtype == QPlotSpectrum::ECD)
            output << "##DATA TYPE= " << "ECD SPECTRUM" << std::endl;
        output << "##BLOCK_ID= " << (i+1) << std::endl;
        output << "##XUNITS= NM" << std::endl;
        qWarning() << "take care here with the units, revise" << Qt::endl;
        output << "##YUNITS= ABSORBANCE" << std::endl;
        output << "##XFACTOR= " << 1.0f << std::endl;
        output << "##YFACTOR= " << 1.0f << std::endl;  //0.01 f should be appropiate
        output << "##FIRSTX= " <<  m_min << std::endl;
        output << "##LASTX= " << m_max << std::endl;
        output << "##NPOINTS= " << m_npoints << std::endl;

        //Put now the data
        output << "##XYDATA= (X++(Y..Y))" << std::endl;
        float increment=(m_max-m_min)/(m_npoints-1);
        for (size_t  i=0; i< d.size();i++)
        {
            output << m_min+i*increment << "   " << std::uppercase << d[i].real() << std::nouppercase<< std::endl;

        }

        output << "##END= " << std::endl;
    }

    return true;
}

void UVSpectrum::CopyData()
{
    QClipboard* c=QApplication::clipboard();
    if ( c == nullptr )
    {
        std::cerr << "Could not open clipboard";
        return;
    }

    std::stringstream s;

    float width= (m_max - m_min);

    const std::vector<fidarray>* datasets=this->GetData();
    const fidarray* tdata=this->GetTotalData();
    std::string sep="\t";

    s << "Lambda(nm)" << sep;
    for(size_t i=0;i<datasets->size();++i)
    {
        s << "#" << (i+1) << sep;

    }
    if ( tdata != nullptr )
    {
        s << "Total" << std::endl;
    }

    size_t npoints=datasets->front().size();
    double delta=width/(npoints-1);
    for(size_t k=0;k<npoints;++k)
    {
        //write x coordinates
        double x=m_min+k*delta;
        s << x << sep;
        for(const auto& d : *datasets)
        {
            s << d[k].real() << sep;
        }

        if ( tdata != nullptr )
        {
            s << (*tdata)[k].real();
        }
        s << std::endl;

    }

    c->setText ( s.str().c_str() );

    return;
}


void UVSpectrum::SetFile(const char* file)
{
    m_file=file;
}


void UVSpectrum::SetLineWidth(float sigma)
{
    //set the global linewidth and each line
    m_linewidth=sigma;
    for(auto& sd : m_sinusoidsets )
    {
        for(auto& s : sd)
        {
            s.m_fldecay=sigma;
        }
    }
}


void UVSpectrum::SetShift(float shift)
{
    m_shift=shift;
    RecalculateX();

}

void UVSpectrum::SubstractSolventShift(bool b)
{
    m_bsubstractsolventshift=b;
    RecalculateX();
}

void UVSpectrum::RecalculateX()
{

    for(size_t i=0;i<m_sinusoidsets.size();++i)
    {
        std::vector<Sinusoid>& sd=m_sinusoidsets[i];
        const std::vector<Spectralline>& ld=m_linesets[i];
        for(size_t j=0;j<sd.size();++j)
        {
            sd[j].m_flfrequency=ld[j].x+m_shift;
            if ( m_bsubstractsolventshift )
            {
                sd[j].m_flfrequency-=ld[j].SolventShift();
            }
        }
    }

}

void UVSpectrum::SetAutoLimits()
{
    //Get the limits of the spectrum
    if ( m_linesets.empty() ) return;
    std::vector<Spectralline>::iterator fit;
    float min=std::numeric_limits<float>::max();
    float max=std::numeric_limits<float>::min();
    const std::vector<Spectralline>& lines=m_linesets.front();

    for(const auto& l : lines )
    {
        if ( l.x > max ) max=l.x;
        if ( l.x < min ) min=l.x;
    }

    min*=0.9;
    max*=1.1;

    //I do not why I hard set this to 100 and 600, check if is better working now
    m_min=min;
    m_max=max;
}

void UVSpectrum::SetEnantiomer(bool b)
{
    m_benantiomer=b;
    RecalculateX();
}

void UVSpectrum::SetFormalism(formalism f)
{
    m_formalism=f;
    RecalculateX();
}
