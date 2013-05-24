#ifndef VECTOR_H
#define VECTOR_H

#include <stdexcept>

template <class T=int>
class Vector
{
public :
  Vector();
  Vector(const unsigned int);
  Vector(const unsigned int, const T);
  Vector(const Vector &);
    ~Vector();
  inline Vector<T> & operator+=(const Vector<T> & right)
  {
    if(right.N != N) throw std::invalid_argument("In Vector::operator+=");
    for(unsigned int i=0;i<N;i++) m_data[i]+=right.m_data[i];
    return *this;
  }

  inline T & operator[](const unsigned int i) const
  {
    return m_data[i];
  }

  inline T & operator[](const unsigned int i)
  {
    return m_data[i];
  }

  inline friend bool operator==(const Vector<T> & v1, const Vector<T> & v2)
  {
    bool is_eq = (v1.N == v2.N);
    if(is_eq)
      {
	for(int i=0;i<v1.N;i++) is_eq *= (v1.m_data[i] == v2.m_data[i]);
      }
    return is_eq;
  }

  inline Vector<T> & operator=(const Vector<T> & v)
  {
    if(!(&v == this))
      {
	N = v.N;
	if(m_data) delete[] m_data;
	m_data = 0;
	m_data = new T[N];
	for(unsigned int i=0;i<N;i++) m_data[i] = v.m_data[i];
      }
    return *this;
  }

  inline friend Vector<T> operator+(const Vector<T> & v1, const Vector<T> & v2)
  {
    if(v1.N != v2.N) throw std::invalid_argument("In Vector::operator+");
    Vector<T> temp;
    temp.N = v1.N;
    for(unsigned int i=0;i<temp.N;i++) temp.m_data[i] = v1.m_data[i] + v2.m_data[i];
    return temp;
  }

  inline friend Vector<T> operator-(const Vector<T> & left, const Vector<T> & right)
  {
    if(left.N != right.N) throw std::invalid_argument("In Vector::operator-");
    Vector<T> temp;
    temp.N = left.N;
    for(unsigned int i=0;i<temp.N;i++) temp.m_data[i] = left.m_data[i] - right.m_data[i];
    return temp;
  }

  inline friend Vector<T> operator*(const T & k, const Vector<T> & v)
  {
    Vector<T> temp;
    temp.N = v.N;
    for(unsigned int i=0;i<temp.N;i++) temp.m_data[i] = k*v.m_data[i];
    return temp;
  }

  unsigned int size() const;
  Vector< Vector<T> > get_base_vectors(const T, const T);
  void resize(const unsigned int n);

private :
  T *m_data;
  unsigned int N;
};

template <class T>
Vector<T>::Vector()
{
  N = 2;
  m_data = new T[N];
}

template <class T>
Vector<T>::Vector(const unsigned int n)
{
  N = n;
  m_data = new T[N];
}

template <class T>
Vector<T>::Vector(const unsigned int n, const T value)
{
  N = n;
  m_data = new T[N];
  for(unsigned int i=0;i<N;i++) m_data[i] = value;
}

template <class T>
Vector<T>::Vector(const Vector<T> & v)
{
  N = v.N;
  m_data = new T[N];
  for(unsigned int i=0;i<N;i++) m_data[i] = v.m_data[i];
}

template <class T>
Vector<T>::~Vector()
{
  if(m_data) delete[] m_data;
  m_data = 0;
}

template <class T>
unsigned int Vector<T>::size() const
{
  return N;
}

template <class T>
Vector< Vector<T> > Vector<T>::get_base_vectors(const T additiveIdentity, const T multiplicativeIdentity)
{
  Vector< Vector<T> > basis(N);
  Vector<T> tempElement(N);
  for(unsigned int d=0;d<N;d++)
    {
      for(unsigned int i=0;i<N;i++)
	{
	  if(i==d) tempElement[i] = multiplicativeIdentity;
	  else tempElement[i] = additiveIdentity;
	}
      basis[d] = tempElement;
    }
  return basis;
}

template <class T>
void Vector<T>::resize(const unsigned int n)
{
    if(n != N)
    {
        N = n;
        if(m_data) delete[] m_data;
        m_data = 0;
        m_data = new T[N];
    }
}

#endif
