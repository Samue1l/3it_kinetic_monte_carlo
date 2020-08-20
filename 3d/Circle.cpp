#include "Circle.h"
#include "libs.h"
Circle::Circle(int xc, int yc, float r) : _xc(xc), _yc(yc), _radius(r)
{
}

bool Circle::is_in(int x, int y) const {
    bool r_value = pow(x-_xc,2) + pow(y-_yc,2) <= pow(_radius,2);
    return r_value;
}
