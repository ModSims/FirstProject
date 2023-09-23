#pragma once

class Geometry {
public:
    Geometry() = default;
    virtual ~Geometry() = default;
    virtual void draw() = 0;
    virtual float getVolume() = 0;
};

class Sphere : public Geometry {
private:
    float _radius;
public:
    Sphere(float radius = 1.0f);
    ~Sphere();
    void draw();
    float getVolume();
};