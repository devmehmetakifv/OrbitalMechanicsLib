#ifndef spaceScience_H
#define spaceScience_H
#pragma once
#include <string>
#include <iostream>
using namespace std;

// spaceScience.cpp : Defines the functions for the static library.
//
class OrbitalFormulas {
public:
    static long double SOI(long double semiMajorAxis, long double smallerObjectMass, long double biggerObjectMass);
    static void OrbitalShape(double eccentricity);
    static long double GravitationalParameter(long double firstObjectMass, long double secondObjectMass);
    static long double OrbitalPeriod(long double gravitationalParameter, long double semiMajorAxis);
    static long double OrbitalEnergy(long double gravitationalParameter, long double semiMajorAxis);
    class OrbitalVelocity {
    public:
        static long double CircularOrbitVelocity(long double planetMass, long double altitudeFromCenter);
        static long double EllipticOrbitVelocity(long double gravitationalParameter, long double distanceBetweenTwoOrbitingBodies, long double semiMajorAxis);
        static long double HyperbolicOrbitVelocity(long double gravitationalParameter, long double distanceBetweenTwoOrbitingBodies, long double semiMajorAxis);
        static long double ParabolicOrbitVelocity(long double planetMass, long double altitudeFromCenter);
        static long double SemiLatusCircular(long double semiMajorAxis);
        static long double SemiLatusElliptical(long double semiMajorAxis, long double eccentricity);
        static long double SemiLatusParabolic(long double periapsisDistance);
        static long double SemiLatusHyperbolic(long double semiMajorAxis, long double eccentricity);
    };
    static long double PeriapsisDistance(std::string orbitShape, long double eccentricity, long double semiMajorAxis);
    class DistanceFromBody {
    public:
        static long double CircularOrbit(long double semiMajorAxis);
        static long double EllipticalOrbit(long double specificAngularMomentum, long double gravitationalParameter, long double eccentricity, long double trueAnomaly);
        static long double ParabolicOrbit(long double periapsisDistance, long double trueAnomaly);
        static long double HyperbolicOrbit(long double specificAngularMomentum, long double gravitationalParameter, long double eccentricity, long double trueAnomaly);
    };
    class ArealVelocity {
    public:
        long double Circular(long double planetMass, long double semiMajorAxis);
        long double Elliptical(long double planetMass, long double semiMajorAxis, long double eccentricity);
        long double Parabolic(long double periapsisDistance, long double planetMass);
        long double Hyperbolic(long double planetMass, long double semiMajorAxis, long double eccentricity);
    };
    class PeriapsisVelocity {
    public:
        long double Circular(long double planetMass, long double semiMajorAxis);
        long double Elliptical(long double planetMass, long double semiMajorAxis, long double eccentricity);
        long double Parabolic(long double planetMass, long double periapsisDistance);
        long double Hyperbolic(long double planetMass, long double semiMajorAxis, long double eccentricity);
    };
};
class MechanicFormulas {
public:
    static double DeltaV(double specificImpulse, double standardGravity, long double initialTotalMass, long double finalTotalMass);
    static long double GravitationalForce(long double smallerObjectMass, long double biggerObjectMass, long double distance);
    static double GravitationalAcceleration(long double planetMass, long double planetRadius);
    static long double PlanetaryMass(double gravitationalAcceleration, long double planetRadius);
    static long double PlanetaryDensity(long double planetaryMass, long double planetRadius);
    static long double PlanetArea(long double planetRadius);
    static long double MetersToAU(long double meter);
    static long double AUToMeters(double au);
    static long double metersToLightYear(long double meters);
    static long double lightYearToMeters(long double lightYears);
};

#endif spaceScience_H