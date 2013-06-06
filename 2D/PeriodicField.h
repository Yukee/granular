#ifndef PFIELD_H
#define PFIELD_H

#include "ScalarField.h"
#include <stdlib.h> // using the abs function in modulo

class PeriodicField : public ScalarField
{
 public:
  virtual double & operator()(const Vector<int> &);
 PeriodicField(Vector<int> range) : ScalarField(range) {}

  /**************************************/

 PeriodicField() : ScalarField() {} // by default is constructed a scalard field on a 2D space, equal to 0 at the origin and not defined elsewhere.
 PeriodicField(const PeriodicField & u) : ScalarField(u) {} // overloading the copy constructor because the default one let both m_data pointers point the same memory
    virtual ~PeriodicField();
    
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
    PeriodicField & operator=(const PeriodicField &);
    bool operator==(const PeriodicField &);
    friend std::ostream & operator<<(std::ostream &, const PeriodicField &);
    friend PeriodicField operator+(const double &, const PeriodicField &);
    friend PeriodicField operator+(const PeriodicField &, const PeriodicField &);
    friend PeriodicField operator-(const double &, const PeriodicField &);
    friend PeriodicField operator-(const PeriodicField &, const PeriodicField &);
    friend PeriodicField operator*(const double &, const PeriodicField &);
    friend PeriodicField operator*(const PeriodicField &, const PeriodicField &);
    friend PeriodicField operator/(const PeriodicField &, const PeriodicField &);
    double & operator[](const int &);
    PeriodicField max_field(const PeriodicField) const;
    PeriodicField module() const;
    double get_max() const;
    void write_in_file(std::ostream & output, const Vector<double> deltaX, const Vector<double> lowerLeftCorner);
    Vector<int> get_pos(int) const;//retrieves the Position of the ith element of m_data

 private:
  int modulo(int, int);
};

#endif
