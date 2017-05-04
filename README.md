# RIMFading

RIMFading is an implementation of the RIMFading Radio Propagation Model [[1](https://www.cs.virginia.edu/wsn/docs/papers/tosn06-models.pdf)] in 2D and 3D in the OMNeT++ Simulator [[2](https://omnetpp.org)]. This implementation work has been done at `Sustainable Communication Networks Group, University of Bremen, Germany`.

Introduction
============

RIMFading is a Radio Propagation Model that characaterizes a radio wave propagation as a function of frequency, distance and other parameters. The model predicts the behavior of propagation with a given constrains and under certain circumstances. The parameter Ki represents the Degree of Irregularity. This parameter elaborates the change of radio range with respect to the direction of propagation.

This code is an implementation of the RIMFading Model in OMNeT++ Simulator for the INET framework [[3](https://inet.omnetpp.org)].


Radio Propagation Models in OMNeT++
==========================
Radio Propagaiton Models in OMNet++ can be seen in INET framewwork under the following path: src/inet/physicallayer/pathloss


To simplify the implementation of RIMFading, a set of base implementations are provided in OMNeT++ that implements some of the abstract methods of pathloss. These base implementations focus on implementing the basic functionality required for pathloss calculation that is used for further calculaction.

This RIMFading implementation extends the FreeSpacePathLoss class.


Compiling the RIM Fading Model
=================================

Place the following 3 files in the following folder of the INET framework and rebuild: src/inet/physicallayer/pathloss

- `RIMFading.ned`
- `RIMFading.h`
- `RIMFading.cc`


Using the RIMFading Model
=============================

There are two approaches for using the RIMFading model:
- Implement a new interface
- Use existing interfaces.

The RIMfading model has a number of configurable parameters that are defined in the 'RIMFading.ned'. Each of these parameters has a default value and if they are required to be changed, use the `omnetpp.ini` file to set these changed values. The following list provides the parameters specific to the RIMFading model.

- 'a' & 'b' - Distribution parameters of a Weibull distribution. These parameters  are used for random number distribution                 that produces floating-point(default a = 1.5 , b = 1); 
- 'model'   - The RIMfading model is implemented in 2D and 3D. This parameter is used to specify the type of the model                     (default is 2D);
- 'k'       - For the fading phenomena the Rician distribution is used. This parameter is used to calculate path loss for                   Rician fading (default is 8dB).
- 'DOI'     - Unlike other propagation models, RIMfading is close to real scenarious. One of the parameter that makes the                   model more realistic is degree of irregularity(DOI). This parameter shows how irregular the radio range is                   (default is 0.006).


Support
=======

If you have any questions or comments regarding this code, please write to,

- Behruz Khalilov (behruz@uni-bremen.de)
- Anas bin Muslim (anas1@uni-bremen.de),
- Asanga Udugama (adu@comnets.uni-bremen.de) or 
- Anna Foerster (anna.foerster@comnets.uni-bremen.de)

