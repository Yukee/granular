#include "Flux.h"

Flux::Flux() :  m_space_dimensions(1), m_solved_dimensions(1)
{
  m_evaluated_flux_d.resize(m_solved_dimensions);
}

Flux::Flux(int space_dimensions, int solved_dimensions) : m_space_dimensions(space_dimensions), m_solved_dimensions(solved_dimensions)
{
  m_evaluated_flux_d.resize(m_solved_dimensions);
}

int Flux::get_space_dimensions()
{
  return m_space_dimensions;
}

int Flux::get_solved_dimensions()
{
  return m_solved_dimensions;
}
