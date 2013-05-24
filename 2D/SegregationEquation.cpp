#include "SegregationEquation.h"
#include "Equation.h"

using namespace std;

SegregationEquation::SegregationEquation(vectorFunction speed_fun, double segregationRate, unsigned int dim) : m_speed_fun(speed_fun), m_sr(segregationRate)
{
        m_flux.resize(dim);
        m_fluxJabobian.resize(dim);
        m_speed.resize(dim);
}
