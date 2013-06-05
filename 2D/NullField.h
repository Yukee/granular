#ifndef NFIELD_H
#define NFIELD_H

#include "ScalarField.h"

class NullField : public ScalarField
{
 public :
  virtual double & operator()(const Vector<int> &);
    NullField(Vector<int> range) : ScalarField(range) {}

 private :
  double zero;
};

#endif
