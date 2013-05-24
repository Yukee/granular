#ifndef C_POS
#define C_POS

#include <iostream>
class Position
{
public :
    Position();
    Position(unsigned int dim);
    Position(int x, int y);
    Position & operator+=(const Position & addPos);
    int operator[](const int i) const;
    int & operator[](const int i);
    friend bool operator==(const Position & p1, const Position & p2);
    friend Position operator+(const Position & p1, const Position & p2);
    friend Position operator*(const int & a, const Position & p);
    friend bool operator>=(const Position & p, const int & m);
    friend bool operator<=(const Position & p, const int & M);

    void set_position(int x, int y);
    void aff_position();
    int get_position(int i) const;
    inline unsigned int get_dimension() const
    {
        return m_dim;
    }

private :
    int m_x;
    int m_y;

    unsigned int m_dim;//dimension of space, here 2 but to be extended
};

#endif
