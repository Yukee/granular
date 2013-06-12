#ifndef NFIELD_H
#define NFIELD_H

#include "ScalarField.h"

class NullField : public ScalarField
{
 public :
  virtual double & operator()(const Vector<int> &);
    NullField(Vector<int> range) : ScalarField(range) {}
    
    NullField() : ScalarField() {}
    NullField(const NullField & u) : ScalarField(u) {}
    virtual ~NullField();
    
    NullField & operator=(const double &);
    
 /*************************/
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
    NullField & operator=(const NullField &);
    bool operator==(const NullField &);
    friend std::ostream & operator<<(std::ostream &, const NullField &);
    friend NullField operator+(const double &, const NullField &);
    friend NullField operator+(const NullField &, const NullField &);
    friend NullField operator-(const double &, const NullField &);
    friend NullField operator-(const NullField &, const NullField &);
    friend NullField operator*(const double &, const NullField &);
    friend NullField operator*(const NullField &, const NullField &);
    friend NullField operator/(const NullField &, const NullField &);
    double & operator[](const int &);
    NullField max_field(const NullField) const;
    NullField module() const;
    double get_max() const;
    void write_in_file(std::ostream & output, const Vector<double> deltaX, const Vector<double> lowerLeftCorner);
    void write_in_file_matrixform(std::ostream & output);
    Vector<int> get_pos(int) const;//retrieves the Position of the ith element of m_data

 private :
  double zero;
};

#endif
