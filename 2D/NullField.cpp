#include "NullField.h"

double & NullField::operator()(const Vector<int> & component)
{

#ifdef DEBUG
    if(component.size() != m_r_len) throw invalid_argument("In PeriodicField::operator(): dimension differs from the field dimension");
#endif // DEBUG

    bool within_ranges = 1;
    for(unsigned int i=0; i<m_r_len; i++) within_ranges*=( (component[i] < m_r[i]) && (component[i] >= 0) );

    if(within_ranges)
      {
	int element = component[m_r_len-1];
	for(unsigned int i=0; i<m_r_len-1; i++)
	  {
	    element*=m_r[m_r_len-2-i];
	    element+=component[m_r_len-2-i];
	  }
	return m_data[element];
      }
    
    zero = 0;
    return zero;
}
