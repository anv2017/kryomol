/*****************************************************************************************
                            ringperceptor.cpp  -  description
                             -------------------
This file is part of the KryoMol project.
For more information, see <http://kryomol.sourceforge.io/>

This program is free software; you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation version 2 of the License.
******************************************************************************************/

#include "ringperceptor.h"
#include "molecule.h"

using namespace kryomol;


RingPerceptor::Edge::Edge()
{}

RingPerceptor::Edge RingPerceptor::Edge::operator+ ( RingPerceptor::Edge& b )
{

    Edge& a= ( *this );
    Edge newedge;
    newedge.Path().reserve ( this->Path().size() +b.Path().size() );

    if ( a.Path().back() == b.Path().back() ) // (a x ) (b x)
    {
        newedge.Path() =a.Path();
        std::vector<size_t>::reverse_iterator it;
        for ( it=b.Path().rbegin() +1;it!=b.Path().rend();++it )
        {
            newedge.Path().push_back ( *it );
        }
        return newedge;
    }

    if ( a.Path().back() == b.Path().front() ) //  (a x ) (x b)
    {
        newedge.Path() =a.Path();
        std::vector<size_t>::const_iterator it;
        for ( it=b.Path().begin() +1;it!=b.Path().end();++it )
        {
            newedge.Path().push_back ( *it );
        }
        return newedge;
    }


    if ( a.Path().front() == b.Path().front() )  // (x a) (x b)
    {
        std::vector<size_t>::reverse_iterator it;
        for ( it=a.Path().rbegin();it!=a.Path().rend();++it )
        {
            std::cout << ( *it ) << std::endl;
            newedge.Path().push_back ( *it );
        }
        std::vector<size_t>::iterator jjt;
        for ( jjt=b.Path().begin() +1;jjt!=b.Path().end();++jjt )
        {
            newedge.Path().push_back ( *jjt );
        }
        return newedge;
    }
    if ( a.Path().front() == b.Path().back() ) // ( x a ) (b x)
    {
        std::vector<size_t>::reverse_iterator it;
        for ( it=a.Path().rbegin();it!=a.Path().rend();++it )
        {
            newedge.Path().push_back ( *it );
        }
        for ( it=b.Path().rbegin() +1;it!=b.Path().rend();++it )
        {
            newedge.Path().push_back ( *it );
        }
        return newedge;
    }

    return newedge;
}

RingPerceptor::RingPerceptor ( const Molecule* molecule ) : m_molecule ( molecule ), m_frame ( NULL )
{}

RingPerceptor::RingPerceptor ( const Frame* frame ) : m_molecule ( NULL ), m_frame ( frame )
{}

RingPerceptor::~RingPerceptor()
{}

std::vector< std::vector<size_t> > RingPerceptor::Rings ( size_t size/*=0*/ )
{
    Convert();
    while ( !m_vertex.empty() )
        Remove();

    //Output Rings
    std::vector< std::vector<size_t> > rings;
    std::vector<Edge>::iterator it;
    for ( it=m_rings.begin();it!=m_rings.end();++it )
    {
        std::vector<size_t>::iterator st;
        if ( size == 0 || size == it->Path().size() -1 )
        {
            rings.push_back ( it->Path() );
            rings.back().pop_back();
        }
    }
    return rings;
}

void RingPerceptor::Convert()
{
    m_vertex.resize ( m_molecule->Atoms().size() );
    for ( size_t i=0;i< m_molecule->Atoms().size();++i )
    {
        m_vertex[i]=i;
    }
    std::vector<Bond>::const_iterator bt;
    const std::vector<Bond>* bonds;
    if ( m_molecule )
        bonds=&m_molecule->Bonds();
    else bonds=&m_frame->Bonds();

    for ( bt=bonds->begin();bt!=bonds->end();++bt )
    {
        m_edges.push_back ( Edge() );
        m_edges.back().Path().push_back ( bt->I() );
        m_edges.back().Path().push_back ( bt->J() );

    }
}



void RingPerceptor::Remove()
{
    size_t vertex=m_vertex.back();

    std::vector<Edge>::iterator it;
    std::vector<Edge> newedges;
    for ( it=m_edges.begin();it!=m_edges.end();++it )
    {
        if ( it->Path().front() == vertex || it->Path().back() ==vertex )
        {
            std::vector<Edge>::iterator it1;
            //bool bfound=false;
            for ( it1=it+1;it1!=m_edges.end();++it1 )
            {

                if ( it1->Path().front() == vertex || it1->Path().back() == vertex )
                {
                    if ( ( it->Path().front() != it->Path().back() ) && ( it1->Path().front() != it1->Path().back() ) )
                    {
                        if ( CommonElements ( *it,*it1 ) == 1 )
                        {
                            newedges.push_back ( ( *it ) + ( *it1 ) );

                        }
                        if ( CommonElements ( *it,*it1 ) ==  2 )
                        {
                            if ( (it->Path().back() == it1->Path().back() && it->Path().front() == it1->Path().front()) ||
                                 (it->Path().back() == it1->Path().front() && it->Path().front() == it1->Path().back()) )
                            {
                                newedges.push_back ( ( *it ) + ( *it1 ) );

                            }
                        }
                    }
                }
                //   if ( bfound ) break;
            }
        }

    }   //append new edges

    bool found;
    do
    {
#ifdef __GNUC__
#warning optimize this loop
#endif
        found=false;
        std::vector<Edge>::iterator it;
        for ( it=m_edges.begin();it!=m_edges.end();++it )
        {
            if ( it->Path().front() ==vertex || it->Path().back() == vertex ) break;
        }
        if ( it!=m_edges.end() )
        {
            if ( it->Path().front() == it->Path().back() ) m_rings.push_back ( *it );
            found=true;
            m_edges.erase ( it );
        }
    }
    while ( found );

    for ( it=newedges.begin();it!=newedges.end();++it )
    {
        m_edges.push_back ( *it );
    }

    m_vertex.pop_back();

}

size_t RingPerceptor::CommonElements ( const Edge& edg1, const Edge& edg2 )
{
    size_t counter=0;
    std::vector<size_t>::const_iterator it;
    std::vector<size_t>::const_iterator it1;
    for ( it=edg1.Path().begin();it!=edg1.Path().end();++it )
    {
        for ( it1=edg2.Path().begin();it1!=edg2.Path().end();++it1 )
        {
#ifdef __GNUC__
#warning optimization
#endif
            if ( ( *it ) == ( *it1 ) ) counter++;
        }
    }
    return counter;
}
