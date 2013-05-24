#include "Position.h"
#include <exception>
#include <stdexcept>

using namespace std;

Position::Position() : m_x(0), m_y(0), m_dim(2)
{
}

Position::Position(unsigned int dim) : m_x(0), m_y(0), m_dim(dim)
{
}

Position::Position(int x, int y) : m_x(x), m_y(y), m_dim(2)
{
}

Position & Position::operator+=(const Position & addPos)
{
    m_x+=addPos.m_x;
    m_y+=addPos.m_y;
    return *this;
}

int Position::operator[](const int i) const
{
    if( i!=0 && i!=1 ) throw invalid_argument("In Position::operator[] argument must be either 0 or 1");
    if(i==0) return m_x;
    if(i==1) return m_y;
}

int & Position::operator[](const int i)
{
    if( i!=0 && i!=1 ) throw invalid_argument("In Position::operator[] argument must be either 0 or 1");
    if(i==0) return m_x;
    if(i==1) return m_y;
}

bool operator==(const Position & p1, const Position & p2)
{
    return p1.m_x == p2.m_x && p1.m_y == p2.m_y;
}

Position operator+(const Position & p1, const Position & p2)
{
    return Position(p1.m_x + p2.m_x, p1.m_y + p2.m_y);
}

Position operator*(const int & a, const Position & p)
{
    return Position(a*p.m_x, a*p.m_y);
}

bool operator>=(const Position & p, const int & m)
{
    return min(p.m_x, p.m_y) >= m;
}

bool operator<=(const Position & p, const int & M)
{
    return max(p.m_x, p.m_y) <= M;
}

void Position::set_position(int x, int y)
{
    m_x=x;
    m_y=y;
}

void Position::aff_position()
{
    cout << m_x << " " << m_y << endl;
}
int Position::get_position(int i) const
{
    if(i==0) return m_x;
    if(i==1) return m_y;
    else return 0;
}
