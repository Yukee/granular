#include "PeriodicField.h"

int PeriodicField::modulo (int m, int n)
{ 
  int mod = m % n;
  if(m < 0) mod = ( n - abs ( mod ) ) % n; 
  return mod; 
}

double & PeriodicField::operator()(const Vector<int> & component)
{

#ifdef DEBUG
    if(component.size() != m_r_len) throw invalid_argument("In PeriodicField::operator(): dimension differs from the field dimension");
#endif // DEBUG

    for(unsigned int i=0; i<m_r_len; i++) component[i] = modulo(component[i], m_r[i]); // periodic

    int element = component[m_r_len-1];
    for(unsigned int i=0; i<m_r_len-1; i++)
    {
        element*=m_r[m_r_len-2-i];
        element+=component[m_r_len-2-i];
    }
    return m_data[element];
}
