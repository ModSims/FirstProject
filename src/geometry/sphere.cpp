#include "geometry.h"
#include <iostream>
#include <cmath>

Sphere::Sphere(float radius) : _radius(radius) {
}

Sphere::~Sphere () {
    std::cout << "Destructing sphere with radius " << _radius << std::endl;
}

void Sphere::draw() {
    std::cout << "Sphere with radius " << _radius << std::endl;
}

float Sphere::getVolume() { 
    return 4.0f/3.0f * 3.14159f * pow(_radius, 3);
}