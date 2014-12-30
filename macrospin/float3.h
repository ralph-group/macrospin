
#include <math.h>

struct float3 {
    float x, y, z;
    
    float3(float x, float y, float z)
        : x(x), y(y), z(z)
    {}

    float3()
        : x(0.0), y(0.0), z(0.0)
    {}

    float3 operator+(const float3 o) {
        return float3(x+o.x, y+o.y, z+o.z);
    }

    float3 operator-(const float3 o) {
        return float3(x-o.x, y-o.y, z-o.z);
    }

    float3 operator*(const float3 o) {
        return float3(x*o.x, y*o.y, z*o.z);
    }

    float3 operator*(const float o) {
        return float3(x*o, y*o, z*o);
    }

    float3 operator/(const float o) {
        return float3(x/o, y/o, z/o);
    }

    float3 cross(const float3& o) {
        return float3(
            y*o.z - z*o.y,
            z*o.x - x*o.z,
            x*o.y - y*o.x
        );
    }

    float dot(const float3& o) {
        return x*o.x + y*o.y + z*o.z;
    }

    float mag() {
        return sqrt(x*x + y*y + z*z);
    }

    void normalize() {
        float mag = this->mag();
        x = x/mag;
        y = y/mag;
        z = z/mag;
    }

};