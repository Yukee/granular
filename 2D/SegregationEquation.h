#ifndef SEGEQ_H
#define SEGEQ_H
#include "Equation.h"
#include "ScalarField.h"

/* In this special case 3 (fu, fv, fw) of the 4 fluxes (fpsi, fu, fv, fw) are null.
So u = u0, v = v0, w = w0 and there is no need to solve the whole system of 4 equations.
We may just solve the psi equation. In that case we need to give the program the known in advance
values of the u, v, w fields, and this is the purpose of this class */

typedef Vector<ScalarField> (*vectorFunction)(Vector<ScalarField>);
typedef ScalarField (*scalarFunction)(Vector<ScalarField>);

class SegregationEquation : public Equation
{
public:
  SegregationEquation(vectorFunction speed_fun, scalarFunction dv, double segregationRate, unsigned int dim); //we need to know the number of space dimensions to specify the number of components of the vectors
    inline void set_speed_field(Vector<ScalarField> position, ScalarField boundaries)
    {
        m_speed = boundaries*m_speed_fun(position);
        m_segregationRate = m_sr*boundaries;
    }
    inline void set_transverse_speed_gradient(Vector<ScalarField> position, ScalarField boundaries)
    {
      m_dv = boundaries*m_transverse_speed_gradient(position);
    }
    inline Vector<ScalarField> get_speed_field()
    {
        return m_speed;
    }
    inline Vector<ScalarField> get_convectionFlux(ScalarField un)
    {
        return m_f(un);
    }
    inline Vector<ScalarField> get_convectionFluxJacobian(ScalarField un)
    {
        return m_df(un);
    }
    inline ScalarField get_sourceTerm(ScalarField un)
    {
      return m_d(un);
    }

private:
    vectorFunction m_speed_fun;
    scalarFunction m_transverse_speed_gradient;
    double m_sr;
    ScalarField m_segregationRate; //the sr is null outside the boundary surface
    Vector<ScalarField> m_speed;
    ScalarField m_dv;
    Vector<ScalarField> m_flux;
    Vector<ScalarField> m_fluxJabobian;

    inline Vector<ScalarField> m_f(ScalarField u)
    {
        m_flux[0]=u*m_speed[0];
        m_flux[1]=u*(m_speed[1]+(-m_sr)*(1+(-1)*u));

        return m_flux;
    }

    inline Vector<ScalarField> m_df(ScalarField u)
    {
        m_fluxJabobian[0]=m_speed[0];
        m_fluxJabobian[1]=m_speed[1]+(-1)*m_segregationRate*(1+(-2)*u);

        return m_fluxJabobian;
    }

    inline ScalarField m_d(ScalarField u)
    {
      return (-1)*u*m_dv;
    }
};
#endif // SEGEQ_H
