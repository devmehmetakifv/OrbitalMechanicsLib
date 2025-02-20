#pragma once

// spaceScience.cpp : Defines the functions for the static library.
//
#include <iostream>
#include <string>
#include <math.h>
#include "spaceScience.h"
#include "pch.h"
#define GRAV_CONSTANT 0.0000000000667
#define AU 149600000000
#define PI 3.14159265359
constexpr auto LY = 9460528400000000;
using std::string;

long double OrbitalFormulas::SOI(long double semiMajorAxis, long double smallerObjectMass, long double biggerObjectMass) {
    long double result;
    result = semiMajorAxis * pow((smallerObjectMass / biggerObjectMass), (2 / 5));
    return result;
}
void OrbitalFormulas::OrbitalShape(double eccentricity) {
    if (eccentricity == 0) {
        std::cout << "Circular Orbit" << std::endl;
        std::cout << "Semi-major axis is equal to radius of the orbit." << std::endl;
        std::cout << "Energy value is smaller than 0." << std::endl;
    }
    else if (eccentricity < 1 && eccentricity > 1) {
        std::cout << "Ellipse Orbit" << std::endl;
        std::cout << "Semi-major axis is bigger than 0." << std::endl;
        std::cout << "Energy value is smaller than 0." << std::endl;
    }
    else if (eccentricity == 1) {
        std::cout << "Parabolic Orbit" << std::endl;
        std::cout << "Semi-major axis goes to infinity." << std::endl;
        std::cout << "Energy is equal to 0." << std::endl;
    }
    else if (eccentricity > 1) {
        std::cout << "Hyperbolic Orbit" << std::endl;
        std::cout << "Semi-major axis is smaller than 0." << std::endl;
        std::cout << "Energy value is bigger than 0." << std::endl;
    }
}
long double OrbitalFormulas::GravitationalParameter(long double firstObjectMass, long double secondObjectMass) {
    long double result;
    result = GRAV_CONSTANT * (firstObjectMass + secondObjectMass);
    return result;
}
long double OrbitalFormulas::OrbitalPeriod(long double gravitationalParameter, long double semiMajorAxis) {
    long double result;
    result = 2 * PI * sqrt(pow(semiMajorAxis, 3) / gravitationalParameter);
    return result;
}
long double OrbitalFormulas::OrbitalEnergy(long double gravitationalParameter, long double semiMajorAxis) {
    long double result;
    result = -1 * (gravitationalParameter / (2 * semiMajorAxis));
    return result;
}

long double OrbitalFormulas::OrbitalVelocity::CircularOrbitVelocity(long double planetMass, long double altitudeFromCenter) {
    long double result;
    result = sqrt((GRAV_CONSTANT * planetMass) / altitudeFromCenter);
    return result;
}
long double OrbitalFormulas::OrbitalVelocity::EllipticOrbitVelocity(long double gravitationalParameter, long double distanceBetweenTwoOrbitingBodies, long double semiMajorAxis) {
    long double result;
    result = sqrt(gravitationalParameter * ((2 / distanceBetweenTwoOrbitingBodies) - (1 / semiMajorAxis)));
    return result;
}
long double OrbitalFormulas::OrbitalVelocity::HyperbolicOrbitVelocity(long double gravitationalParameter, long double distanceBetweenTwoOrbitingBodies, long double semiMajorAxis) {
    long double result;
    result = sqrt(gravitationalParameter * ((2 / distanceBetweenTwoOrbitingBodies) + abs(1 / semiMajorAxis)));
    return result;
}
long double OrbitalFormulas::OrbitalVelocity::ParabolicOrbitVelocity(long double planetMass, long double altitudeFromCenter) {
    long double result;
    result = sqrt((2 * GRAV_CONSTANT * planetMass) / altitudeFromCenter);
    return result;
}
long double OrbitalFormulas::OrbitalVelocity::SemiLatusCircular(long double semiMajorAxis) {
    long double result;
    result = semiMajorAxis;
    return result;
}
long double OrbitalFormulas::OrbitalVelocity::SemiLatusElliptical(long double semiMajorAxis, long double eccentricity) {
    long double result;
    result = semiMajorAxis * (1 - pow(eccentricity, 2));
    return result;
}
long double OrbitalFormulas::OrbitalVelocity::SemiLatusParabolic(long double periapsisDistance) {
    long double result;
    result = 2 * periapsisDistance;
    return result;
}
long double OrbitalFormulas::OrbitalVelocity::SemiLatusHyperbolic(long double semiMajorAxis, long double eccentricity) {
    long double result;
    result = semiMajorAxis * (1 - pow(eccentricity, 2));
    return result;
}

long double OrbitalFormulas::PeriapsisDistance(string orbitShape, long double eccentricity, long double semiMajorAxis) {
    if (orbitShape == "Circular") {
        long double result;
        result = semiMajorAxis;
        return result;
    }
    else if (orbitShape == "Elliptic") {
        long double result;
        result = semiMajorAxis * (1 - eccentricity);
        return result;
    }
    else if (orbitShape == "Parabolic") {
        long double result;
        result = NULL;
        return result;
    }
    else if (orbitShape == "Hyperbolic") {
        long double result;
        result = semiMajorAxis * (1 - eccentricity);
        return result;
    }
    else if (orbitShape != "Circular" || "Elliptic" || "Parabolic" || "Hyperbolic") {
        std::cout << "You've entered invalid type of orbit shape. Try followings: 'Circular', 'Elliptic', 'Parabolic', 'Hyperbolic'" << std::endl;
        return NULL;
    }
}
long double OrbitalFormulas::DistanceFromBody::CircularOrbit(long double semiMajorAxis) {
    long double result;
    result = semiMajorAxis;
    return result;
}
long double OrbitalFormulas::DistanceFromBody::EllipticalOrbit(long double specificAngularMomentum, long double gravitationalParameter, long double eccentricity, long double trueAnomaly) {
    long double result;
    long double p;
    p = pow(specificAngularMomentum, 2) / gravitationalParameter;
    result = p * (1 / (1 + eccentricity * cos(trueAnomaly)));
    return result;
}
long double OrbitalFormulas::DistanceFromBody::ParabolicOrbit(long double periapsisDistance, long double trueAnomaly) {
    long double result;
    result = 2 * periapsisDistance * (1 / 1 + cos(trueAnomaly));
    return result;
}
long double OrbitalFormulas::DistanceFromBody::HyperbolicOrbit(long double specificAngularMomentum, long double gravitationalParameter, long double eccentricity, long double trueAnomaly) {
    long double result;
    long double p;
    p = pow(specificAngularMomentum, 2) / gravitationalParameter;
    result = p * (1 / 1 + eccentricity * cos(trueAnomaly));
    return result;
}
long double OrbitalFormulas::ArealVelocity::Circular(long double planetMass, long double semiMajorAxis) {
    long double result;
    result = sqrt(GRAV_CONSTANT * planetMass * semiMajorAxis);
    return result;
}
long double OrbitalFormulas::ArealVelocity::Elliptical(long double planetMass, long double semiMajorAxis, long double eccentricity) {
    long double result;
    result = sqrt(GRAV_CONSTANT * planetMass * semiMajorAxis * ((1 + eccentricity) / (1 - eccentricity)));
    return result;
}
long double OrbitalFormulas::ArealVelocity::Parabolic(long double periapsisDistance, long double planetMass) {
    long double result;
    result = sqrt((GRAV_CONSTANT * planetMass * periapsisDistance) / 2);
    return result;
}
long double OrbitalFormulas::ArealVelocity::Hyperbolic(long double planetMass, long double semiMajorAxis, long double eccentricity) {
    long double result;
    result = sqrt(-GRAV_CONSTANT * planetMass * semiMajorAxis * ((1 + eccentricity) / (eccentricity - 1)));
    return result;
}
long double OrbitalFormulas::PeriapsisVelocity::Circular(long double planetMass, long double semiMajorAxis) {
    long double result;
    result = sqrt((GRAV_CONSTANT * planetMass) / semiMajorAxis);
    return result;
}
long double OrbitalFormulas::PeriapsisVelocity::Elliptical(long double planetMass, long double semiMajorAxis, long double eccentricity) {
    long double result;
    result = sqrt(((GRAV_CONSTANT * planetMass) / semiMajorAxis) * ((1 + eccentricity) / (1 - eccentricity)));
    return result;
}
long double OrbitalFormulas::PeriapsisVelocity::Parabolic(long double planetMass, long double periapsisDistance) {
    long double result;
    result = sqrt((2 * GRAV_CONSTANT * planetMass) / periapsisDistance);
    return result;
}
long double OrbitalFormulas::PeriapsisVelocity::Hyperbolic(long double planetMass, long double semiMajorAxis, long double eccentricity) {
    long double result;
    result = sqrt(((-GRAV_CONSTANT * planetMass) / semiMajorAxis) * ((1 + eccentricity) / (eccentricity - 1)));
    return result;
}
double MechanicFormulas::DeltaV(double specificImpulse, double standardGravity, long double initialTotalMass, long double finalTotalMass) {
    double result;
    result = specificImpulse * standardGravity * log(initialTotalMass / finalTotalMass);
    return result;
}
long double MechanicFormulas::GravitationalForce(long double smallerObjectMass, long double biggerObjectMass, long double distance) {
    long double result;
    result = (GRAV_CONSTANT * biggerObjectMass * smallerObjectMass) / pow(distance, 2);
    return result;
}
double MechanicFormulas::GravitationalAcceleration(long double planetMass, long double planetRadius) {
    double result;
    result = (GRAV_CONSTANT * planetMass) / pow(planetRadius, 2);
    return result;
}
long double MechanicFormulas::PlanetaryMass(double gravitationalAcceleration, long double planetRadius) {
    long double result;
    result = (gravitationalAcceleration * pow(planetRadius, 2)) / GRAV_CONSTANT;
    return result;
}
long double MechanicFormulas::PlanetaryDensity(long double planetaryMass, long double planetRadius) {
    long double result;
    long double volume;
    volume = (4 * PI * pow(planetRadius, 3)) / 3;
    result = (planetaryMass / volume) / 1000;
    return result;
}
long double MechanicFormulas::PlanetArea(long double planetRadius) {
    long double result;
    result = 4 * PI * pow(planetRadius, 2);
    return result;
}
long double MechanicFormulas::MetersToAU(long double meter) {
    long double result;
    result = meter / AU;
    return result;
}
long double MechanicFormulas::AUToMeters(double au) {
    long double result;
    result = au * AU;
    return result;
}
long double MechanicFormulas::metersToLightYear(long double meters) {
    long double result;
    result = meters / LY;
    return result;
}
long double MechanicFormulas::lightYearToMeters(long double lightYears) {
    long double result;
    result = lightYears * LY;
    return result;
}
