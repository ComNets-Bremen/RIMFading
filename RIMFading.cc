/******************************************************************************
 * RIMFading - A RIMFading implementation for the INET Framework of the OMNeT++
 * Simulator.
 *
 * Copyright (C) 2016, Sustainable Communication Networks, University of Bremen, Germany
 *
 * This program is free software; you can redistribute it and/or modify it
 * under the terms of the GNU General Public License as published by the
 * Free Software Foundation; version 3 of the License.
 *
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY
 * or FITNESS FOR A PARTICULAR PURPOSE.
 * See the GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License along
 * with this program; if not, see <http://www.gnu.org/licenses/>
 *
 *
 ******************************************************************************/

/**
 *  The C++ implementation file of the RIMFading Radio Propagation Model for the INET Framework
 * in OMNeT++.
 *
 * @authors : Behruz Khalilov (behruz@uni-bremen.de), Anas bin Muslim (anas1@uni-bremen.de)
 *
 */
#include <iostream>
#include <fstream>
#include "inet/physicallayer/pathloss/RIMFading.h"
#include "inet/common/ModuleAccess.h"
#include "inet/physicallayer/contract/packetlevel/IRadioMedium.h"

using namespace std;

namespace inet {

namespace physicallayer {

Define_Module(RIMFading);

RIMFading::RIMFading() :
            a(1.5),
            b(1),
            model(2),
            DOI()

{
}

void RIMFading::initialize(int stage)
{
    FreeSpacePathLoss::initialize(stage);
    if (stage == INITSTAGE_LOCAL) {
        a=par("a");
        b=par("b");
        model=par("model");
        DOI=par("DOI");

    }

    ofstream myfile;

    //Irregularity varies for each degree => ki varies
    double ki = 1;
    double ran = weibull(a,b);
    int sign = 0;

    // this loop generates Ki and checks condition |K0-k359|<=DOI; if so, result is recorded
    do{
        for(int i=0; i<360; i++){
            sign = intuniform(0, 1);
            if(sign==0){sign = sign -1;}
            ran = weibull(a,b);
            if(i==0)
                ki=1;
            else{
                ki += sign*(ran * DOI);
            }
            DifferenceInPathLoss[i] = ki;
        }
    } while(abs(DifferenceInPathLoss[0]-DifferenceInPathLoss[359])>DOI);
    for(int i=0; i<360; i++){
        //write the result of coefficient of irregularity in the file
        myfile.open("DifferenceInPathLoss.txt",std::ios::app);
        myfile << DifferenceInPathLoss[i] <<"\n";
        myfile.close();
    }
}

std::ostream& RIMFading::printToStream(std::ostream& stream, int level) const
{
    stream << "RIMFading";
    if (level >= PRINT_LEVEL_TRACE)
        stream << ", alpha = " << alpha
        << ", system loss = " << systemLoss
        << ", a = " << a
        << ", b = " << b
        << ", Model = " << model
        << ", DOI = " << DOI;
    return stream;
}

double RIMFading::computeAngles2D(const ITransmission *transmission, const IArrival *arrival) const
{
    const Coord transmissionStartPosition = transmission->getStartPosition();
    const Coord arrivalStartPosition = arrival->getStartPosition();

    double tx_x = transmissionStartPosition.x;
    double tx_y = transmissionStartPosition.y;

    double rx_x = arrivalStartPosition.x;
    double rx_y = arrivalStartPosition.y;

    double dx = rx_x - tx_x;
    double dy = rx_y - tx_y;
    double angleRxTx = 0;

    if(dx == 0) //Both points along a vertical line
    {
        if (dy > 0)
            angleRxTx = 0;
        else if (dy < 0)
            angleRxTx = 180;
        else
            EV << "Nodes at the same point; No angle \n";
    }
    else if (dx > 0) // Rx to the right of Tx
    {
        if (dy == 0) //Both points along a horizontal line
            angleRxTx = 90;
        else if (dy < 0) //Rx above Tx
            angleRxTx = 90 - abs(180/M_PI*atan(dy/dx));
        else // dy > 0, Rx below Tx
            angleRxTx = 90 + (180/M_PI*atan(dy/dx));
    }
    else // dx < 0; Rx to the left of Tx
    {
        if (dy == 0) //Both points along a horizontal line
            angleRxTx = 270;
        else if (dy < 0) //Rx above Tx
            angleRxTx = 270 + (180/M_PI*atan(dy/dx));
        else // dy > 0, Rx below Tx
            angleRxTx = 270 - abs(180/M_PI*atan(dy/dx));
    }

    return angleRxTx;
}


double RIMFading::computeAngles3D(const ITransmission *transmission, const IArrival *arrival, double *angle) const
{
    const Coord transmissionStartPosition = transmission->getStartPosition();
    const Coord arrivalStartPosition = arrival->getStartPosition();

    double x1 = transmissionStartPosition.x;
    double y1 = transmissionStartPosition.y;
    double z1 = transmissionStartPosition.z;

    double x2 = arrivalStartPosition.x;
    double y2 = arrivalStartPosition.y;
    double z2 = arrivalStartPosition.z;

    //If Nodes are present in 2D then use 2D, no need for using 3D Model
    if(z1==z2){
        *angle = 0;
        return computeAngles2D(transmission,arrival);
    }

    //Shift the origin to transmitter's location & calculate new coordinates for receiver
    x2=x2-x1;
    y2=y2-y1;
    z2=z2-z1;

    //Calculations for angles
    double phi, theta;
    theta = computeAngles2D(transmission,arrival);
    phi = (atan(sqrt(pow(x2,2) + pow(y2,2))/z2)) * (180/M_PI);
    *angle = phi;

    return theta;
}


/*
 * Path Loss calculation of RIM Propagation Model
 */

double RIMFading::RIMPathLossCalculation(double freeSpacePathLoss, int iter) const
{
    ofstream myfile;

    EV<<"Path Loss : "<<freeSpacePathLoss*DifferenceInPathLoss[iter-1]<<"\n";

    //write the result of a path-loss in the file

    myfile.open("pathloss.txt",std::ios::app);
    myfile << "Path Loss : " << freeSpacePathLoss*DifferenceInPathLoss[iter-1] <<"\n";
    myfile.close();

    return freeSpacePathLoss*DifferenceInPathLoss[iter-1];

}

double RIMFading::computePathLoss(const ITransmission *transmission, const IArrival *arrival) const {
    auto radioMedium = transmission->getTransmitter()->getMedium();
    auto narrowbandSignalAnalogModel = check_and_cast<const INarrowbandSignal *>(transmission->getAnalogModel());
    auto transmitterPosition = transmission->getStartPosition();
    auto recepiverPosition = arrival->getStartPosition();
    mps propagationSpeed = radioMedium->getPropagation()->getPropagationSpeed();
    Hz carrierFrequency = narrowbandSignalAnalogModel->getCarrierFrequency();
    m distance = m(recepiverPosition.distance(transmitterPosition));
    m waveLength = propagationSpeed / carrierFrequency;


    double theta = 0, loopIterations = 0, phi=0;
    double *angle = &phi;

    double freeSpacePathLoss = computeFreeSpacePathLoss(waveLength, distance, alpha, systemLoss);

    if(model==2){
        theta = computeAngles2D(transmission,arrival);
        phi=0;
    }
    else if(model==3){
        theta = computeAngles3D(transmission,arrival,angle);
    }
    else{
        EV<<"Please give dimension model either in 2D or 3D"<<"\n";
    }
    loopIterations = theta + phi;
    //convert double to int -- since angle can be any double value it should be rounded to integer value
    int resangle= nearbyint(loopIterations);

   double resultPL = RIMPathLossCalculation(freeSpacePathLoss, resangle);

    return resultPL;
}

} // namespace physicallayer

} // namespace inet
