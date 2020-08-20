#ifndef CIRCLE_H
#define CIRCLE_H


class Circle
{
    public:
        Circle(int xc, int yc, float r);
        bool is_in(int x, int y) const;

    private:
        const float _radius;
        const int _xc, _yc;

};

#endif // CIRCLE_H
