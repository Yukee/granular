#ifndef FD1SOLVER_H
#define	FD1SOLVER_H

#include "Equation.h"
#include "Vector.h"
#include "ScalarField.h"

typedef ScalarField (*boundarySurface)(VectorField); //returns a field which is 0 outside the surface, 1 inside

enum grid_position {lft, rght};
enum bound {lower, upper};
enum boundaryConditions {periodic, null, prescribedWestAndSouth};

class FD1Solver
{
public:
    FD1Solver(int spaceDimension, Vector<double> deltaX, Vector<double> xInterval, Equation *eq, boundaryConditions bc, boundarySurface bs, Vector<double> lowerLeftCorner);
    ScalarField get_numerical_flux_gradient(ScalarField un);

    inline VectorField get_position()
    {
        return m_position;
    }

    inline VectorField get_resized_position()
    {
        return m_resizedPos;
    }

    double check_CFL(double deltaT);
    ~FD1Solver();
    void set_un(ScalarField un);
    inline void set_uWest(Vector<double> uw)
    {
        uWest = uw;
    }

    inline void set_uSouth(Vector<double> us)
    {
        uSouth = us;
    }

    inline ScalarField get_test_boundaries()
    {
        return m_test_boundaries;
    }

    inline Vector<double> get_deltaX()
    {
        return m_deltaX;
    }

    inline Vector<double> get_lowerLeftCorner()
    {
        return m_lowerLeftCorner;
    }

    inline VectorField get_speed_field()
    {
        return m_eq->get_speed_field();
    }

    inline VectorField get_flux_jacobian(ScalarField u0)
    {
        get_numerical_flux_gradient(u0);
        return m_eq->get_convectionFluxJacobian(upper_right_intermediate_un_values[0]);
    }

    inline ScalarField get_initial_field(ScalarField u0)
    {
        get_numerical_flux_gradient(u0);
        return m_un;
    }

private:
    int m_spaceDimension;
    Vector<double> m_deltaX;
    Vector<double> m_xInterval;
    Equation *m_eq;
    boundaryConditions m_bc;
    boundarySurface m_bs;
    Vector<double> m_lowerLeftCorner;
    ScalarField m_test_boundaries;//the evalutation of m_bs over all positions m_position
    Vector<int> m_nxSteps;
    ScalarField m_un;
    VectorField left_convection_flux;//convection_flux[i] = flux in the i direction of space
    VectorField right_convection_flux;
    VectorField upper_right_intermediate_un_values;//un values cell boundary
    VectorField lower_right_intermediate_un_values;
    VectorField upper_left_intermediate_un_values;
    VectorField lower_left_intermediate_un_values;
    VectorField m_position;
    VectorField m_resizedPos;
    VectorField right_localSpeed;
    VectorField left_localSpeed;

    void compute_intermediate_un_values();
    void compute_localSpeed();
    void compute_numerical_convection_flux();
    double minmod(double a, double b);
    double three_pts_derivative(Vector<int> j, int direction); //du/dx at point j in the x or y direction
    double intermediate_uxt_values(Vector<int> j, grid_position p, bound b, int direction); //un(+-)(j+-1/2) at point j in the x or y direction
    double un(Vector<int> j);
    void resize_pos();//compute a position vector m_resizedPos with non extended BC

    //For the prescribed boundary case
    Vector<double> uWest;
    Vector<double> uSouth;

    Vector<int> x;
    Vector<int> y;
};


#endif
