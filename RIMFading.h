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
* The C++ include file of the RIMFading Radio Propagation Model for the INET Framework
* in OMNeT++.
*
* @authors : Behruz Khalilov (behruz@uni-bremen.de), Anas bin Muslim (anas1@uni-bremen.de)
*
*/
#ifndef __INET_RIMFading_H
#define __INET_RIMFADING_H

#include "inet/common/INETDefs.h"
#include "math.h"
#include "inet/physicallayer/pathloss/FreeSpacePathLoss.h"
#include "inet/environment/contract/IPhysicalEnvironment.h"
namespace inet {

namespace physicallayer {

using namespace inet::physicalenvironment;
/**
 * This class implements the Radio Irregularity Path Loss model
 */
class INET_API RIMFading : public FreeSpacePathLoss
{
protected:
    double a,b,model;
    double DifferenceInPathLoss[360];


protected:
    virtual void initialize(int stage) override;


private:
    double DOI; // Degree of Irregularity

 //   double DifferenceInPathLoss[360];// Array to storing the difference in path loss "Ki" on start

public:
    RIMFading();
    virtual std::ostream& printToStream(std::ostream& stream, int level) const override;

    /**
     * returns the loss factor for the provided transmission and arrival.
     */
    virtual double computePathLoss(const ITransmission *transmission, const IArrival *arrival) const override;

    /**
     * returns the loss factor as a function of propagation speed, carrier, frequency and distance.
     */

    virtual double RIMPathLossCalculation(double, int) const;
    virtual double computeAngles2D(const ITransmission *transmission, const IArrival *arrival) const;
    virtual double computeAngles3D(const ITransmission *transmission, const IArrival *arrival, double *angle) const;

};

} // namespace physicallayer

} // namespace inet

#endif // ifndef __INET_RIMFADING_H
