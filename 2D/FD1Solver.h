#ifndef FD1SOLVER_H
#define	FD1SOLVER_H

#include "Equation.h"
#include "Vector.h"
#include "ScalarField.h"
#include "PeriodicField.h"

enum grid_position {lft, rght};
enum bound {lower, upper};

class FD1Solver
{
public:
    FD1Solver(Vector<double> deltaX, Vector<double> xInterval, Equation *eq, Flux *bs, Vector<double> lowerLeftCorner);

    VectorField get_numerical_flux_gradient(VectorField un);

    inline VectorField get_position()
    {
        return m_position;
    }

    double check_CFL(double deltaT);

    ~FD1Solver();

    inline Vector<double> get_deltaX()
    {
        return m_deltaX;
    }

    inline Vector<double> get_lowerLeftCorner()
    {
        return m_lowerLeftCorner;
    }

    inline Vector<TensorField> get_flux_jacobian(VectorField u0)
    {
        get_numerical_flux_gradient(u0);
        return m_eq->get_convectionFluxJacobian(upper_right_intermediate_un_values[0]);
    }

    inline VectorField get_initial_field(VectorField u0)
    {
        get_numerical_flux_gradient(u0);
        return m_un;
    }

private:
    int m_m; //number of solved dimensions
    int m_n; //number of space dimensions

    Vector< Vector<int> > m_b; // base vectors of the physical space

    Vector<double> m_deltaX;
    Vector<double> m_xInterval;
    Equation *m_eq;
    Flux *m_bs;
    Vector<double> m_lowerLeftCorner;
    Vector<int> m_nxSteps;

    VectorField m_position;

    VectorField m_un;

    TensorField left_convection_flux;
    TensorField right_convection_flux;
    TensorField left_diffusion_flux;
    TensorField right_diffusion_flux;

    TensorField un_derivatives;
    TensorField upper_right_intermediate_un_values;
    TensorField lower_right_intermediate_un_values;
    TensorField upper_left_intermediate_un_values;
    TensorField lower_left_intermediate_un_values;

    VectorField right_localSpeed;
    VectorField left_localSpeed;

    void compute_un_derivatives();
    void compute_intermediate_un_values();
    void compute_localSpeed();
    void compute_numerical_convection_flux();
    void compute_numerical_diffusion_flux();
    double minmod(double a, double b);
    double three_pts_derivative(Vector<int> j, int d, int i); //du/dx at point j in the x or y direction

    // caca
    ScalarField unity;
};


#endif
