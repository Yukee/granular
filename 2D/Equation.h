#ifndef EQUATION_H
#define	EQUATION_H

#include "Vector.h"
#include "ScalarField.h"

//first argument : u(t), second : x(t)
typedef Vector<ScalarField> VectorField;
typedef VectorField (*convectionFlux)(ScalarField);
typedef VectorField (*convectionFluxJacobian)(ScalarField);

class Equation
{
protected:
    convectionFlux m_f;
    convectionFluxJacobian m_df;

public:
    Equation();
	Equation(convectionFlux f, convectionFluxJacobian df) : m_f(f), m_df(df){}
	virtual VectorField get_convectionFlux(ScalarField un) = 0;
	virtual VectorField get_convectionFluxJacobian(ScalarField un) = 0;
	virtual void set_speed_field(VectorField, ScalarField) = 0;
	virtual VectorField get_speed_field() = 0;
};


#endif


