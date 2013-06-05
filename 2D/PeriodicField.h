#ifndef PFIELD_H
#define PFIELD_H

#include "ScalarField.h"
#include <stdlib.h> // using the abs function in modulo

class PeriodicField : public ScalarField
{
 public:
  virtual double & operator()(const Vector<int> &);
 PeriodicField(Vector<int> range) : ScalarField(range) {}

 private:
  int modulo(int, int);
};

#endif
