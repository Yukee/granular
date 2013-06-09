#ifndef PRFIELD_H
#define PRFIELD_H

#include "ScalarField.h"

class PrescribedField : public ScalarField
{
 public:
  virtual double & operator()(Vector<int> );

 PrescribedField() {
	 PrescribedField(Vector<int> (2,1));}
	 
 PrescribedField(Vector<int> range): ScalarField(range) {
	 m_bounds.resize(2*m_r_len);
	 Vector<int> range_surf = m_r;
	 Vector< Vector<int> > b = range_surf.get_base_vectors(1,0); for(unsigned int i=0;i<m_r_len;i++) range_surf = range_surf + b[i];
	 
	 for(unsigned int d=0;d<m_r_len;d++)
	 {
		 m_bounds[2*d].resize_field(range_surf.drop(d));
		 m_bounds[2*d+1].resize_field(range_surf.drop(d));
	 }
 }
	 
 PrescribedField(const PrescribedField & u): ScalarField(u) {
	 m_bounds = u.m_bounds;}
	 
 virtual ~PrescribedField();
 
 // sets surface orthogonal to the d axis of the cube (i is the direction, positive if i=1, negative if i=-1)
 void set_bound(const int d, const int i, const ScalarField & u); 
    
  /**************************************/

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
    PrescribedField & operator=(const PrescribedField &);
    bool operator==(const PrescribedField &);
    friend std::ostream & operator<<(std::ostream &, const PrescribedField &);
    friend PrescribedField operator+(const double &, const PrescribedField &);
    friend PrescribedField operator+(const PrescribedField &, const PrescribedField &);
    friend PrescribedField operator-(const double &, const PrescribedField &);
    friend PrescribedField operator-(const PrescribedField &, const PrescribedField &);
    friend PrescribedField operator*(const double &, const PrescribedField &);
    friend PrescribedField operator*(const PrescribedField &, const PrescribedField &);
    friend PrescribedField operator/(const PrescribedField &, const PrescribedField &);
    double & operator[](const int &);
    PrescribedField max_field(const PrescribedField) const;
    PrescribedField module() const;
    double get_max() const;
    void write_in_file(std::ostream & output, const Vector<double> deltaX, const Vector<double> lowerLeftCorner);
    void write_in_file_matrixform(std::ostream & output);
    Vector<int> get_pos(int) const;//retrieves the Position of the ith element of m_data  
    
 private:
	Vector<ScalarField> m_bounds; // value of the field on each face of the cube
};

#endif
