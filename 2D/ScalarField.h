#ifndef SCALARFIELD_H
#define SCALARFIELD_H

#include "Vector.h"
#include <cstring>
#include <exception>
#include <stdexcept>
#include <iostream>

class ScalarField
{

public:
  ScalarField(); // by default is constructed a scalard field on a 2D space, equal to 0 at the origin and not defined elsewhere.
    ScalarField(const ScalarField &); // overloading the copy constructor because the default one let both m_data pointers point the same memory
    ScalarField(Vector<int> range); // creates a scalar field on a range.size()D space, defined on a rectangle of lengths range[i]. Note that the data is NOT initialized.
    virtual ~ScalarField();
    
    void resize_field(Vector<int> range);
    inline unsigned int get_space_dimension() const
    {
        return m_r_len;
    }
    inline Vector<int> get_range() const
    {
        return m_r;
    }
    inline int get_size() const
    {
        return m_data_len;
    }
    ScalarField & operator=(const ScalarField &);
    bool operator==(const ScalarField &);
    friend std::ostream & operator<<(std::ostream &, const ScalarField &);
    friend ScalarField operator+(const double &, const ScalarField &);
    friend ScalarField operator+(const ScalarField &, const ScalarField &);
    friend ScalarField operator-(const double &, const ScalarField &);
    friend ScalarField operator-(const ScalarField &, const ScalarField &);
    friend ScalarField operator*(const double &, const ScalarField &);
    friend ScalarField operator*(const ScalarField &, const ScalarField &);
    friend ScalarField operator/(const ScalarField &, const ScalarField &);
    virtual double & operator()(const Vector<int> & component);
    double & operator[](const int &);
    ScalarField max_field(const ScalarField) const;
    ScalarField module() const;
    double get_max() const;
    void write_in_file(std::ostream & output, const Vector<double> deltaX, const Vector<double> lowerLeftCorner);
    Vector<int> get_pos(int) const;//retrieves the Position of the ith element of m_data

protected:
    Vector<int> m_r;//has as much components as the number of dimensions of the space, each component is the number of points in a direction of the space
    unsigned int m_r_len;//lenght of the m_r array; corresponds to the dimension of the space where the field lives
    double *m_data;
    unsigned int m_data_len;//lenght of the m_data array; corresponds to the number of points of the field times the dimension of the space
};

#endif // SCALARFIELD_H
