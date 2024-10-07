/*****************************************************************************************
                            density.h  -  description
                             -------------------
This file is part of the KryoMol project.
For more information, see <http://kryomol.sourceforge.io/>

This program is free software; you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation version 2 of the License.
******************************************************************************************/

#ifndef PHYSICALCONSTANTS_H
#define PHYSICALCONSTANTS_H

class PC
{
public:
    PC();

    static constexpr double Pi = 3.1415926535897932384626433832795;
    static constexpr double AngstromToBohr = 1.889726132885643067222711130708;
    static constexpr double StandDevUniform = 0.28867513459481288225457439025098; //Standard Deviation of an uniform distribution
    static constexpr double Tan30 = 0.57735026918962576450914878050196;
    static constexpr double Sqrt15Div15 = 0.25819888974716112567861769331883;
    static constexpr double ConstDY0 = 0.31539156525252000603089369029571; //0.25*sqrt(5/PC::Pi)
    static constexpr double ConstDY1 = 1.0925484305920790705433857058027; //(0.5)*sqrt(15/(PC::Pi))
    static constexpr double ConstDY2 = 1.0925484305920790705433857058027; //(0.5)*sqrt(15/(PC::Pi))
    static constexpr double ConstDY3 = 0.54627421529603953527169285290134; //0.25*sqrt(15/(PC::Pi))
    static constexpr double ConstDY4 = 1.0925484305920790705433857058027; //0.5*sqrt(15/(PC::Pi))
    static constexpr double ConstFY0 = 0.373176332590115391414395913199; //0.25*sqrt(7/PC::Pi)
    static constexpr double ConstFY1 = 0.45704579946446573615802069691665; //(0.25)*sqrt(21/(2*PC::Pi))
    static constexpr double ConstFY2 = 0.45704579946446573615802069691665; //(0.25)*sqrt(21/(2*PC::Pi))
    static constexpr double ConstFY3 = 1.445305721320277027694690077199; //0.25*sqrt(105/(PC::Pi))
    static constexpr double ConstFY4 = 2.890611442640554055389380154398; //0.5*sqrt(105/(PC::Pi))
    static constexpr double ConstFY5 = 0.5900435899266435103456102775415; //(0.25)*sqrt(35/(2*PC::Pi)
    static constexpr double ConstFY6 = 0.5900435899266435103456102775415; //(0.25)*sqrt(35/(2*PC::Pi))
    static constexpr double h = 6.62607015e-34;
    static constexpr double jtoev = 6.241509074e18;
    static constexpr double c = 299792458;


};

#endif // PHYSICALCONSTANTS_H
