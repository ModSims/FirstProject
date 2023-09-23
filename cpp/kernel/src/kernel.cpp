#include "geometry.h"

#include <iostream>
#include <iomanip>
#include <vector>
#include <memory>
#include <algorithm>
#include <numeric>

void printGeo(Geometry* g) {
    std::cout << "This Geometry has Volume: \t";
    std::cout << g->getVolume() << "\n";
}

int main(int argc, char** argv) {
    std::cout << std::fixed;
    std::cout << std::setprecision(10);

    std::vector<std::unique_ptr<Geometry>> geometryList;

    geometryList.push_back(std::make_unique<Sphere>(24.0f));
    geometryList.push_back(std::make_unique<Sphere>(36.0f));

    std::for_each(
        geometryList.begin(), 
        geometryList.end(), 
        [](const std::unique_ptr<Geometry> &g){
            printGeo(g.get());
        }
    );

    std::cout << std::endl;

    return 0;
}
