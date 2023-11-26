/*****************************************************************************************
                            pcmodelparser.cpp  -  description
                             -------------------
This file is part of the KryoMol project.
For more information, see <http://kryomol.sourceforge.io/>

This program is free software; you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation version 2 of the License.
******************************************************************************************/


#include <iostream>
#include <sstream>
#include "pcmodelparser.h"
#include "molecule.h"
#include "stringtools.h"

using namespace kryomol;

PCModelParser::PCModelParser(const char* file) : Parser(file)
{}

PCModelParser::PCModelParser(std::istream* stream) : Parser(stream)
{}


PCModelParser::~PCModelParser()
{}


bool PCModelParser::ParseFile(std::streampos pos)
{
    std::cout << "parsing" << std::endl;
    std::string line;
    m_file->clear();
    m_file->seekg(pos,std::ios::beg);
    bool find_mol=false;
    bool find_type=false;
    while ( (std::getline ( *m_file,line )) && (find_mol == false))
    {
        if ( line.find ( "{PCM " ) != std::string::npos )
        {
            find_mol=true;
        }
    }
    while ( (std::getline ( *m_file,line )) && (find_mol == true) && (find_type == false))
    {

        if ( line.find ( "mmenergy " ) != std::string::npos )
        {
            StringTokenizer token(line, " \t\r");
            if (token.size()>=2)
            {
                std::string forcefield=kryomol::toupper(token.at(1));
                std::cout << "forcefield is " << forcefield << std::endl;
                TypeTable type=UNKNOWN;
                if (forcefield == "AMBER")
                    type = AMBER;
                if (forcefield == "MM2TEST")
                    type = MM2TEST;
                if (forcefield == "MM3")
                    type = MM3;
                if (forcefield == "MMFF94")
                    type = MMFF94;
                if (forcefield == "MMXCONST")
                    type = MMXCONST;
                if (forcefield == "OPLSAA")
                    type = OPLSAA;
                if ( type == UNKNOWN )
                    throw kryomol::Exception("Invalid force field in pcmodel file");
                BuildAtomTable(type);
                find_type=true;
            }
            else throw kryomol::Exception("Could not determine force field in pcmodel file");
        }
    }
    if (find_mol && find_type)
    {
        m_file->clear();
        m_file->seekg(pos,std::ios::beg);
        Molecules()->push_back(Molecule());
        GetMolecule();

        m_file->clear();
        m_file->seekg(pos,std::ios::beg);
        while ( std::getline ( *m_file,line ) )
        {
            if ( line.find ( "{PCM " ) != std::string::npos )
            {
                Molecule& molecule=Molecules()->back();
                molecule.Frames().push_back(Frame(&molecule));
                StringTokenizer tok(line," \t\r");
                molecule.Frames().back().PotentialEnergy()=std::atof(tok.at(2).c_str());
                GetConformer();
            }
        }
    }
    return true;

}


void PCModelParser::GetMolecule()
{
    std::string line;
    Molecule& molecule=Molecules()->back();
    while ((std::getline(*m_file,line)) && (line.find ("mmenergy ") == std::string::npos))
    {
        if (line.find ("AT ") != std::string::npos)
        {
            StringTokenizer token(line," \t\r");
            if (token.size() >= 3)
            {
                Atom atom(ZFromString(token.at(2)));
                molecule.Atoms().push_back(atom);
            }
        }
    }
}

bool PCModelParser::GetConformer()
{
    std::string line;
    Molecule& molecule=Molecules()->back();
    Frame& frame= molecule.Frames().back();
    while ((std::getline(*m_file,line)) && (line.find ("dipole_moment") == std::string::npos))
    {
        if (line.find ("AT ") != std::string::npos)
        {
            StringTokenizer token(line," \t\r");
            if ( (token.size() >= 6) &&  kryomol::isnum(token.at(3)) && kryomol::isnum(token.at(4)) && kryomol::isnum(token.at(5)) )
            {
                Coordinate c;
                c.x()=kryomol::atof(token.at(3));
                c.y()=kryomol::atof(token.at(4));
                c.z()=kryomol::atof(token.at(5));

                frame.XYZ().push_back(c);
            }
        }
    }

    return true;
}

void PCModelParser::BuildAtomTable(TypeTable type)
{
    m_atomtable.clear();

    switch (type)
    {
    case AMBER:
        m_atomtable.insert(AtomPair("1",6));
        m_atomtable.insert(AtomPair("2",6));
        m_atomtable.insert(AtomPair("3",6));
        m_atomtable.insert(AtomPair("4",6));
        m_atomtable.insert(AtomPair("5",6));
        m_atomtable.insert(AtomPair("6",6));
        m_atomtable.insert(AtomPair("7",6));
        m_atomtable.insert(AtomPair("8",6));
        m_atomtable.insert(AtomPair("9",6));
        m_atomtable.insert(AtomPair("10",6));
        m_atomtable.insert(AtomPair("11",6));
        m_atomtable.insert(AtomPair("12",6));
        m_atomtable.insert(AtomPair("13",6));
        m_atomtable.insert(AtomPair("14",7));
        m_atomtable.insert(AtomPair("15",7));
        m_atomtable.insert(AtomPair("16",7));
        m_atomtable.insert(AtomPair("17",7));
        m_atomtable.insert(AtomPair("18",7));
        m_atomtable.insert(AtomPair("19",7));
        m_atomtable.insert(AtomPair("20",7));
        m_atomtable.insert(AtomPair("21",8));
        m_atomtable.insert(AtomPair("22",8));
        m_atomtable.insert(AtomPair("23",8));
        m_atomtable.insert(AtomPair("24",8));
        m_atomtable.insert(AtomPair("25",8));
        m_atomtable.insert(AtomPair("26",16));
        m_atomtable.insert(AtomPair("27",16));
        m_atomtable.insert(AtomPair("28",15));
        m_atomtable.insert(AtomPair("29",1));
        m_atomtable.insert(AtomPair("30",1));
        m_atomtable.insert(AtomPair("31",1));
        m_atomtable.insert(AtomPair("32",1));
        m_atomtable.insert(AtomPair("33",1));
        m_atomtable.insert(AtomPair("34",1));
        m_atomtable.insert(AtomPair("35",1));
        m_atomtable.insert(AtomPair("36",1));
        m_atomtable.insert(AtomPair("37",1));
        m_atomtable.insert(AtomPair("38",1));
        m_atomtable.insert(AtomPair("39",1));
        m_atomtable.insert(AtomPair("40",1));
        m_atomtable.insert(AtomPair("41",8));
        m_atomtable.insert(AtomPair("42",1));
        m_atomtable.insert(AtomPair("43",1));
        m_atomtable.insert(AtomPair("50",18));

        break;
    case MM2TEST:
        m_atomtable.insert(AtomPair("1",6));
        m_atomtable.insert(AtomPair("2",6));
        m_atomtable.insert(AtomPair("3",6));
        m_atomtable.insert(AtomPair("4",6));
        m_atomtable.insert(AtomPair("5",1));
        m_atomtable.insert(AtomPair("6",8));
        m_atomtable.insert(AtomPair("7",8));
        m_atomtable.insert(AtomPair("8",7));
        m_atomtable.insert(AtomPair("9",7));
        m_atomtable.insert(AtomPair("10",7));
        m_atomtable.insert(AtomPair("11",9));
        m_atomtable.insert(AtomPair("12",17));
        m_atomtable.insert(AtomPair("13",35));
        m_atomtable.insert(AtomPair("14",53));
        m_atomtable.insert(AtomPair("15",16));
        m_atomtable.insert(AtomPair("16",16));
        m_atomtable.insert(AtomPair("17",16));
        m_atomtable.insert(AtomPair("18",16));
        m_atomtable.insert(AtomPair("19",14));
        m_atomtable.insert(AtomPair("21",1));
        m_atomtable.insert(AtomPair("22",6));
        m_atomtable.insert(AtomPair("23",1));
        m_atomtable.insert(AtomPair("24",1));
        m_atomtable.insert(AtomPair("25",15));
        m_atomtable.insert(AtomPair("26",5));
        m_atomtable.insert(AtomPair("27",5));
        m_atomtable.insert(AtomPair("28",1));
        m_atomtable.insert(AtomPair("29",6));
        m_atomtable.insert(AtomPair("30",6));
        m_atomtable.insert(AtomPair("31",32));
        m_atomtable.insert(AtomPair("32",50));
        m_atomtable.insert(AtomPair("33",82));
        m_atomtable.insert(AtomPair("34",34));
        m_atomtable.insert(AtomPair("35",52));
        m_atomtable.insert(AtomPair("36",1));
        m_atomtable.insert(AtomPair("37",7));
        m_atomtable.insert(AtomPair("38",16));
        m_atomtable.insert(AtomPair("39",34));
        m_atomtable.insert(AtomPair("40",6));
        m_atomtable.insert(AtomPair("41",7));
        m_atomtable.insert(AtomPair("42",8));
        m_atomtable.insert(AtomPair("43",5));
        m_atomtable.insert(AtomPair("44",13));
        m_atomtable.insert(AtomPair("45",1));
        m_atomtable.insert(AtomPair("46",8));
        m_atomtable.insert(AtomPair("47",15));
        m_atomtable.insert(AtomPair("48",6));
        m_atomtable.insert(AtomPair("49",6));
        m_atomtable.insert(AtomPair("50",6));
        m_atomtable.insert(AtomPair("51",6));
        m_atomtable.insert(AtomPair("52",6));
        m_atomtable.insert(AtomPair("53",8));
        m_atomtable.insert(AtomPair("54",53));
        m_atomtable.insert(AtomPair("55",7));
        m_atomtable.insert(AtomPair("56",6));
        m_atomtable.insert(AtomPair("57",6));
        m_atomtable.insert(AtomPair("58",13));
        m_atomtable.insert(AtomPair("60",1));
        m_atomtable.insert(AtomPair("61",6));
        m_atomtable.insert(AtomPair("62",6));
        m_atomtable.insert(AtomPair("63",6));
        m_atomtable.insert(AtomPair("64",8));
        m_atomtable.insert(AtomPair("65",7));
        m_atomtable.insert(AtomPair("66",8));
        m_atomtable.insert(AtomPair("67",15));
        m_atomtable.insert(AtomPair("68",7));
        m_atomtable.insert(AtomPair("70",1));
        m_atomtable.insert(AtomPair("71",6));
        m_atomtable.insert(AtomPair("72",7));
        m_atomtable.insert(AtomPair("73",9));
        m_atomtable.insert(AtomPair("74",17));
        m_atomtable.insert(AtomPair("75",35));
        m_atomtable.insert(AtomPair("76",53));
        m_atomtable.insert(AtomPair("77",7));

        break;
    case MM3:
        m_atomtable.insert(AtomPair("1",6));
        m_atomtable.insert(AtomPair("2",6));
        m_atomtable.insert(AtomPair("3",6));
        m_atomtable.insert(AtomPair("4",6));
        m_atomtable.insert(AtomPair("5",1));
        m_atomtable.insert(AtomPair("6",8));
        m_atomtable.insert(AtomPair("7",8));
        m_atomtable.insert(AtomPair("8",7));
        m_atomtable.insert(AtomPair("9",7));
        m_atomtable.insert(AtomPair("10",7));
        m_atomtable.insert(AtomPair("11",9));
        m_atomtable.insert(AtomPair("12",17));
        m_atomtable.insert(AtomPair("13",35));
        m_atomtable.insert(AtomPair("14",53));
        m_atomtable.insert(AtomPair("15",16));
        m_atomtable.insert(AtomPair("16",16));
        m_atomtable.insert(AtomPair("17",16));
        m_atomtable.insert(AtomPair("18",16));
        m_atomtable.insert(AtomPair("19",14));
        m_atomtable.insert(AtomPair("21",1));
        m_atomtable.insert(AtomPair("22",6));
        m_atomtable.insert(AtomPair("23",1));
        m_atomtable.insert(AtomPair("24",1));
        m_atomtable.insert(AtomPair("25",15));
        m_atomtable.insert(AtomPair("26",5));
        m_atomtable.insert(AtomPair("27",5));
        m_atomtable.insert(AtomPair("28",1));
        m_atomtable.insert(AtomPair("29",6));
        m_atomtable.insert(AtomPair("30",6));
        m_atomtable.insert(AtomPair("31",32));
        m_atomtable.insert(AtomPair("32",50));
        m_atomtable.insert(AtomPair("33",82));
        m_atomtable.insert(AtomPair("34",34));
        m_atomtable.insert(AtomPair("35",52));
        m_atomtable.insert(AtomPair("36",1));
        m_atomtable.insert(AtomPair("37",7));
        m_atomtable.insert(AtomPair("38",6));
        m_atomtable.insert(AtomPair("39",7));
        m_atomtable.insert(AtomPair("40",7));
        m_atomtable.insert(AtomPair("41",8));
        m_atomtable.insert(AtomPair("42",16));
        m_atomtable.insert(AtomPair("43",7));
        m_atomtable.insert(AtomPair("44",1));
        m_atomtable.insert(AtomPair("45",7));
        m_atomtable.insert(AtomPair("46",7));
        m_atomtable.insert(AtomPair("47",8));
        m_atomtable.insert(AtomPair("48",1));
        m_atomtable.insert(AtomPair("49",8));
        m_atomtable.insert(AtomPair("50",6));
        m_atomtable.insert(AtomPair("51",2));
        m_atomtable.insert(AtomPair("52",10));
        m_atomtable.insert(AtomPair("53",18));
        m_atomtable.insert(AtomPair("54",36));
        m_atomtable.insert(AtomPair("55",54));
        m_atomtable.insert(AtomPair("56",6));
        m_atomtable.insert(AtomPair("57",6));
        m_atomtable.insert(AtomPair("58",6));
        m_atomtable.insert(AtomPair("59",12));
        m_atomtable.insert(AtomPair("60",15));
        m_atomtable.insert(AtomPair("61",26));
        m_atomtable.insert(AtomPair("62",26));
        m_atomtable.insert(AtomPair("63",28));
        m_atomtable.insert(AtomPair("64",28));
        m_atomtable.insert(AtomPair("65",27));
        m_atomtable.insert(AtomPair("66",27));
        m_atomtable.insert(AtomPair("67",6));
        m_atomtable.insert(AtomPair("68",6));
        m_atomtable.insert(AtomPair("69",8));
        m_atomtable.insert(AtomPair("70",8));
        m_atomtable.insert(AtomPair("71",6));
        m_atomtable.insert(AtomPair("72",7));
        m_atomtable.insert(AtomPair("73",1));
        m_atomtable.insert(AtomPair("74",16));
        m_atomtable.insert(AtomPair("75",8));
        m_atomtable.insert(AtomPair("76",8));
        m_atomtable.insert(AtomPair("77",8));
        m_atomtable.insert(AtomPair("78",8));
        m_atomtable.insert(AtomPair("79",8));
        m_atomtable.insert(AtomPair("80",8));
        m_atomtable.insert(AtomPair("81",8));
        m_atomtable.insert(AtomPair("82",8));
        m_atomtable.insert(AtomPair("83",8));
        m_atomtable.insert(AtomPair("84",8));
        m_atomtable.insert(AtomPair("85",8));
        m_atomtable.insert(AtomPair("86",8));
        m_atomtable.insert(AtomPair("87",8));
        m_atomtable.insert(AtomPair("88",8));
        m_atomtable.insert(AtomPair("89",8));
        m_atomtable.insert(AtomPair("90",8));
        m_atomtable.insert(AtomPair("91",8));
        m_atomtable.insert(AtomPair("92",8));
        m_atomtable.insert(AtomPair("93",8));
        m_atomtable.insert(AtomPair("94",8));
        m_atomtable.insert(AtomPair("95",8));
        m_atomtable.insert(AtomPair("96",8));
        m_atomtable.insert(AtomPair("97",8));
        m_atomtable.insert(AtomPair("98",8));
        m_atomtable.insert(AtomPair("99",8));
        m_atomtable.insert(AtomPair("100",8));
        m_atomtable.insert(AtomPair("101",8));
        m_atomtable.insert(AtomPair("102",8));
        m_atomtable.insert(AtomPair("103",8));
        m_atomtable.insert(AtomPair("104",16));
        m_atomtable.insert(AtomPair("105",16));
        m_atomtable.insert(AtomPair("106",6));
        m_atomtable.insert(AtomPair("107",7));
        m_atomtable.insert(AtomPair("108",7));
        m_atomtable.insert(AtomPair("109",7));
        m_atomtable.insert(AtomPair("110",7));
        m_atomtable.insert(AtomPair("111",7));
        m_atomtable.insert(AtomPair("113",6));
        m_atomtable.insert(AtomPair("114",6));
        m_atomtable.insert(AtomPair("115",8));
        m_atomtable.insert(AtomPair("116",8));
        m_atomtable.insert(AtomPair("117",8));
        m_atomtable.insert(AtomPair("118",8));
        m_atomtable.insert(AtomPair("119",8));
        m_atomtable.insert(AtomPair("120",8));
        m_atomtable.insert(AtomPair("121",8));
        m_atomtable.insert(AtomPair("124",1));
        m_atomtable.insert(AtomPair("125",20));
        m_atomtable.insert(AtomPair("126",38));
        m_atomtable.insert(AtomPair("127",56));
        m_atomtable.insert(AtomPair("128",57));
        m_atomtable.insert(AtomPair("129",58));
        m_atomtable.insert(AtomPair("130",59));
        m_atomtable.insert(AtomPair("131",60));
        m_atomtable.insert(AtomPair("132",61));
        m_atomtable.insert(AtomPair("133",62));
        m_atomtable.insert(AtomPair("134",63));
        m_atomtable.insert(AtomPair("135",64));
        m_atomtable.insert(AtomPair("136",65));
        m_atomtable.insert(AtomPair("137",66));
        m_atomtable.insert(AtomPair("138",67));
        m_atomtable.insert(AtomPair("139",68));
        m_atomtable.insert(AtomPair("140",69));
        m_atomtable.insert(AtomPair("141",70));
        m_atomtable.insert(AtomPair("142",71));
        m_atomtable.insert(AtomPair("143",7));
        m_atomtable.insert(AtomPair("144",7));
        m_atomtable.insert(AtomPair("145",8));
        m_atomtable.insert(AtomPair("146",7));
        m_atomtable.insert(AtomPair("148",8));
        m_atomtable.insert(AtomPair("149",8));
        m_atomtable.insert(AtomPair("150",7));
        m_atomtable.insert(AtomPair("151",7));
        m_atomtable.insert(AtomPair("153",15));
        m_atomtable.insert(AtomPair("154",16));
        m_atomtable.insert(AtomPair("155",7));
        m_atomtable.insert(AtomPair("156",6));
        m_atomtable.insert(AtomPair("157",6));
        m_atomtable.insert(AtomPair("158",6));
        m_atomtable.insert(AtomPair("159",8));
        m_atomtable.insert(AtomPair("160",6));
        m_atomtable.insert(AtomPair("161",6));
        m_atomtable.insert(AtomPair("162",6));

        break;
    case MMFF94:
        m_atomtable.insert(AtomPair("1",6));
        m_atomtable.insert(AtomPair("2",6));
        m_atomtable.insert(AtomPair("3",6));
        m_atomtable.insert(AtomPair("4",6));
        m_atomtable.insert(AtomPair("5",1));
        m_atomtable.insert(AtomPair("6",8));
        m_atomtable.insert(AtomPair("7",8));
        m_atomtable.insert(AtomPair("8",7));
        m_atomtable.insert(AtomPair("9",7));
        m_atomtable.insert(AtomPair("10",7));
        m_atomtable.insert(AtomPair("11",9));
        m_atomtable.insert(AtomPair("12",17));
        m_atomtable.insert(AtomPair("13",35));
        m_atomtable.insert(AtomPair("14",53));
        m_atomtable.insert(AtomPair("15",16));
        m_atomtable.insert(AtomPair("16",16));
        m_atomtable.insert(AtomPair("17",16));
        m_atomtable.insert(AtomPair("18",16));
        m_atomtable.insert(AtomPair("19",14));
        m_atomtable.insert(AtomPair("20",6));
        m_atomtable.insert(AtomPair("21",1));
        m_atomtable.insert(AtomPair("22",6));
        m_atomtable.insert(AtomPair("23",1));
        m_atomtable.insert(AtomPair("24",1));
        m_atomtable.insert(AtomPair("25",15));
        m_atomtable.insert(AtomPair("26",15));
        m_atomtable.insert(AtomPair("27",1));
        m_atomtable.insert(AtomPair("28",1));
        m_atomtable.insert(AtomPair("29",1));
        m_atomtable.insert(AtomPair("30",6));
        m_atomtable.insert(AtomPair("31",1));
        m_atomtable.insert(AtomPair("32",8));
        m_atomtable.insert(AtomPair("33",1));
        m_atomtable.insert(AtomPair("34",7));
        m_atomtable.insert(AtomPair("35",8));
        m_atomtable.insert(AtomPair("36",1));
        m_atomtable.insert(AtomPair("37",6));
        m_atomtable.insert(AtomPair("38",7));
        m_atomtable.insert(AtomPair("39",7));
        m_atomtable.insert(AtomPair("40",7));
        m_atomtable.insert(AtomPair("41",6));
        m_atomtable.insert(AtomPair("42",7));
        m_atomtable.insert(AtomPair("43",7));
        m_atomtable.insert(AtomPair("44",16));
        m_atomtable.insert(AtomPair("45",7));
        m_atomtable.insert(AtomPair("46",7));
        m_atomtable.insert(AtomPair("47",8));
        m_atomtable.insert(AtomPair("48",8));
        m_atomtable.insert(AtomPair("49",8));
        m_atomtable.insert(AtomPair("50",1));
        m_atomtable.insert(AtomPair("51",8));
        m_atomtable.insert(AtomPair("52",1));
        m_atomtable.insert(AtomPair("53",7));
        m_atomtable.insert(AtomPair("54",7));
        m_atomtable.insert(AtomPair("55",7));
        m_atomtable.insert(AtomPair("56",7));
        m_atomtable.insert(AtomPair("57",6));
        m_atomtable.insert(AtomPair("58",7));
        m_atomtable.insert(AtomPair("59",8));
        m_atomtable.insert(AtomPair("60",6));
        m_atomtable.insert(AtomPair("61",7));
        m_atomtable.insert(AtomPair("62",7));
        m_atomtable.insert(AtomPair("63",6));
        m_atomtable.insert(AtomPair("64",6));
        m_atomtable.insert(AtomPair("65",7));
        m_atomtable.insert(AtomPair("66",7));
        m_atomtable.insert(AtomPair("67",7));
        m_atomtable.insert(AtomPair("68",7));
        m_atomtable.insert(AtomPair("69",7));
        m_atomtable.insert(AtomPair("70",8));
        m_atomtable.insert(AtomPair("71",1));
        m_atomtable.insert(AtomPair("72",16));
        m_atomtable.insert(AtomPair("73",16));
        m_atomtable.insert(AtomPair("74",16));
        m_atomtable.insert(AtomPair("75",15));
        m_atomtable.insert(AtomPair("76",7));
        m_atomtable.insert(AtomPair("77",17));
        m_atomtable.insert(AtomPair("78",6));
        m_atomtable.insert(AtomPair("79",7));
        m_atomtable.insert(AtomPair("80",6));
        m_atomtable.insert(AtomPair("81",7));
        m_atomtable.insert(AtomPair("82",7));
        m_atomtable.insert(AtomPair("87",26));
        m_atomtable.insert(AtomPair("88",26));
        m_atomtable.insert(AtomPair("89",9));
        m_atomtable.insert(AtomPair("90",17));
        m_atomtable.insert(AtomPair("91",35));
        m_atomtable.insert(AtomPair("92",3));
        m_atomtable.insert(AtomPair("93",11));
        m_atomtable.insert(AtomPair("94",19));
        m_atomtable.insert(AtomPair("95",30));
        m_atomtable.insert(AtomPair("96",20));
        m_atomtable.insert(AtomPair("97",29));
        m_atomtable.insert(AtomPair("98",29));
        m_atomtable.insert(AtomPair("99",12));
        m_atomtable.insert(AtomPair("100",8));
        m_atomtable.insert(AtomPair("101",8));

        break;
    case MMXCONST:
        m_atomtable.insert(AtomPair("1",6));
        m_atomtable.insert(AtomPair("2",6));
        m_atomtable.insert(AtomPair("3",6));
        m_atomtable.insert(AtomPair("4",6));
        m_atomtable.insert(AtomPair("5",1));
        m_atomtable.insert(AtomPair("6",8));
        m_atomtable.insert(AtomPair("7",8));
        m_atomtable.insert(AtomPair("8",7));
        m_atomtable.insert(AtomPair("9",7));
        m_atomtable.insert(AtomPair("10",7));
        m_atomtable.insert(AtomPair("11",9));
        m_atomtable.insert(AtomPair("12",17));
        m_atomtable.insert(AtomPair("13",35));
        m_atomtable.insert(AtomPair("14",53));
        m_atomtable.insert(AtomPair("15",16));
        m_atomtable.insert(AtomPair("16",16));
        m_atomtable.insert(AtomPair("17",16));
        m_atomtable.insert(AtomPair("18",16));
        m_atomtable.insert(AtomPair("19",14));
        m_atomtable.insert(AtomPair("20",0));
        m_atomtable.insert(AtomPair("21",1));
        m_atomtable.insert(AtomPair("22",6));
        m_atomtable.insert(AtomPair("23",1));
        m_atomtable.insert(AtomPair("24",1));
        m_atomtable.insert(AtomPair("25",15));
        m_atomtable.insert(AtomPair("26",5));
        m_atomtable.insert(AtomPair("27",5));
        m_atomtable.insert(AtomPair("28",1));
        m_atomtable.insert(AtomPair("29",6));
        m_atomtable.insert(AtomPair("30",6));
        m_atomtable.insert(AtomPair("31",32));
        m_atomtable.insert(AtomPair("32",50));
        m_atomtable.insert(AtomPair("33",82));
        m_atomtable.insert(AtomPair("34",34));
        m_atomtable.insert(AtomPair("35",52));
        m_atomtable.insert(AtomPair("36",1));
        m_atomtable.insert(AtomPair("37",7));
        m_atomtable.insert(AtomPair("38",16));
        m_atomtable.insert(AtomPair("39",34));
        m_atomtable.insert(AtomPair("40",6));
        m_atomtable.insert(AtomPair("41",7));
        m_atomtable.insert(AtomPair("42",8));
        m_atomtable.insert(AtomPair("43",5));
        m_atomtable.insert(AtomPair("44",13));
        m_atomtable.insert(AtomPair("45",1));
        m_atomtable.insert(AtomPair("46",8));
        m_atomtable.insert(AtomPair("47",15));
        m_atomtable.insert(AtomPair("48",6));
        m_atomtable.insert(AtomPair("49",6));
        m_atomtable.insert(AtomPair("50",6));
        m_atomtable.insert(AtomPair("51",6));
        m_atomtable.insert(AtomPair("52",6));
        m_atomtable.insert(AtomPair("53",8));
        m_atomtable.insert(AtomPair("54",53));
        m_atomtable.insert(AtomPair("55",7));
        m_atomtable.insert(AtomPair("56",6));
        m_atomtable.insert(AtomPair("57",6));
        m_atomtable.insert(AtomPair("58",13));
        m_atomtable.insert(AtomPair("60",1));
        m_atomtable.insert(AtomPair("61",6));
        m_atomtable.insert(AtomPair("62",6));
        m_atomtable.insert(AtomPair("63",6));
        m_atomtable.insert(AtomPair("64",8));
        m_atomtable.insert(AtomPair("65",7));
        m_atomtable.insert(AtomPair("66",8));
        m_atomtable.insert(AtomPair("67",15));
        m_atomtable.insert(AtomPair("68",7));
        m_atomtable.insert(AtomPair("70",1));
        m_atomtable.insert(AtomPair("71",6));
        m_atomtable.insert(AtomPair("72",7));
        m_atomtable.insert(AtomPair("73",9));
        m_atomtable.insert(AtomPair("74",17));
        m_atomtable.insert(AtomPair("75",35));
        m_atomtable.insert(AtomPair("76",53));
        m_atomtable.insert(AtomPair("77",7));
        m_atomtable.insert(AtomPair("78",18));

        break;
    case OPLSAA:
        m_atomtable.insert(AtomPair("1",6));
        m_atomtable.insert(AtomPair("2",6));
        m_atomtable.insert(AtomPair("3",6));
        m_atomtable.insert(AtomPair("4",6));
        m_atomtable.insert(AtomPair("5",6));
        m_atomtable.insert(AtomPair("6",6));
        m_atomtable.insert(AtomPair("7",6));
        m_atomtable.insert(AtomPair("8",6));
        m_atomtable.insert(AtomPair("9",6));
        m_atomtable.insert(AtomPair("10",6));
        m_atomtable.insert(AtomPair("11",6));
        m_atomtable.insert(AtomPair("12",6));
        m_atomtable.insert(AtomPair("13",6));
        m_atomtable.insert(AtomPair("14",6));
        m_atomtable.insert(AtomPair("15",6));
        m_atomtable.insert(AtomPair("16",6));
        m_atomtable.insert(AtomPair("17",6));
        m_atomtable.insert(AtomPair("18",6));
        m_atomtable.insert(AtomPair("19",6));
        m_atomtable.insert(AtomPair("20",6));
        m_atomtable.insert(AtomPair("21",6));
        m_atomtable.insert(AtomPair("22",6));
        m_atomtable.insert(AtomPair("23",6));
        m_atomtable.insert(AtomPair("24",6));
        m_atomtable.insert(AtomPair("25",6));
        m_atomtable.insert(AtomPair("26",6));
        m_atomtable.insert(AtomPair("27",6));
        m_atomtable.insert(AtomPair("28",6));
        m_atomtable.insert(AtomPair("29",1));
        m_atomtable.insert(AtomPair("30",1));
        m_atomtable.insert(AtomPair("31",1));
        m_atomtable.insert(AtomPair("32",1));
        m_atomtable.insert(AtomPair("33",1));
        m_atomtable.insert(AtomPair("34",1));
        m_atomtable.insert(AtomPair("35",1));
        m_atomtable.insert(AtomPair("36",1));
        m_atomtable.insert(AtomPair("37",1));
        m_atomtable.insert(AtomPair("38",1));
        m_atomtable.insert(AtomPair("39",1));
        m_atomtable.insert(AtomPair("40",1));
        m_atomtable.insert(AtomPair("41",1));
        m_atomtable.insert(AtomPair("42",1));
        m_atomtable.insert(AtomPair("43",1));
        m_atomtable.insert(AtomPair("44",8));
        m_atomtable.insert(AtomPair("45",8));
        m_atomtable.insert(AtomPair("46",8));
        m_atomtable.insert(AtomPair("47",8));
        m_atomtable.insert(AtomPair("48",8));
        m_atomtable.insert(AtomPair("49",8));
        m_atomtable.insert(AtomPair("50",8));
        m_atomtable.insert(AtomPair("51",8));
        m_atomtable.insert(AtomPair("52",8));
        m_atomtable.insert(AtomPair("53",7));
        m_atomtable.insert(AtomPair("54",7));
        m_atomtable.insert(AtomPair("55",7));
        m_atomtable.insert(AtomPair("56",7));
        m_atomtable.insert(AtomPair("57",7));
        m_atomtable.insert(AtomPair("58",7));
        m_atomtable.insert(AtomPair("59",7));
        m_atomtable.insert(AtomPair("60",7));
        m_atomtable.insert(AtomPair("61",7));
        m_atomtable.insert(AtomPair("62",9));
        m_atomtable.insert(AtomPair("63",9));
        m_atomtable.insert(AtomPair("64",9));
        m_atomtable.insert(AtomPair("65",9));
        m_atomtable.insert(AtomPair("66",9));
        m_atomtable.insert(AtomPair("67",9));
        m_atomtable.insert(AtomPair("68",16));
        m_atomtable.insert(AtomPair("69",16));
        m_atomtable.insert(AtomPair("70",3));
        m_atomtable.insert(AtomPair("71",11));
        m_atomtable.insert(AtomPair("72",19));
        m_atomtable.insert(AtomPair("73",37));
        m_atomtable.insert(AtomPair("74",55));
        m_atomtable.insert(AtomPair("75",12));
        m_atomtable.insert(AtomPair("76",20));
        m_atomtable.insert(AtomPair("77",38));
        m_atomtable.insert(AtomPair("78",56));
        m_atomtable.insert(AtomPair("79",9));
        m_atomtable.insert(AtomPair("80",17));
        m_atomtable.insert(AtomPair("81",35));
        m_atomtable.insert(AtomPair("82",2));
        m_atomtable.insert(AtomPair("83",10));
        m_atomtable.insert(AtomPair("84",18));
        m_atomtable.insert(AtomPair("85",36));
        m_atomtable.insert(AtomPair("86",54));
        m_atomtable.insert(AtomPair("87",7));

        break;
    default:
        break;
    }
}

const PCModelParser::AtomTable PCModelParser::GetAtomTable()
{
    return m_atomtable;
}

int PCModelParser::ZFromString(const std::string symbol)
{
    const AtomTable table=GetAtomTable();
    std::map<const std::string,const int>::const_iterator it;
    for ( it=table.begin();it!=table.end();++it )
    {
        if ( it->first == symbol )
        {
            return it->second;
        }
    }
    return 0;
}


