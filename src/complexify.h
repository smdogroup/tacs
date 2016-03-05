/***
 *  For the complex step derivative method:
 *  f'(x) ~  Im [ f(x+ih) ] / h
 *  Define a double complex class that inherits from the
 *  library complex type and overloads appropriate operators.
 *  Mon Jan  8 22:42:20 PST 2001
 ***/

#ifndef COMPLEXIFY_H
#define COMPLEXIFY_H

#include <complex>

class cplx : public std::complex<double> {
public:
  cplx() : std::complex<double>() {};
  cplx(const double& d) : std::complex<double>(d) {};
  cplx(const double& r, const double& i) : std::complex<double>(r,i) {};
  cplx(const std::complex<double>& z) : std::complex<double>(z) {};
  cplx(const std::complex<float>& z) : std::complex<double>(z) {};
  operator double() {return this->real();}
  operator int() {return int(this->real());}
  // relational operators
  // Conversion constructor should be able to take care of the
  // operator== and != calls with double, but MIPS compiler 
  // complains of ambiguous inheritance.  This should be more
  // efficient anyway.  (A hint of what comes below.)
  friend inline bool operator==(const cplx&,const cplx&);
  friend inline bool operator==(const cplx&,const double&);
  friend inline bool operator==(const double&,const cplx&);
  friend inline bool operator!=(const cplx&,const cplx&);
  friend inline bool operator!=(const cplx&,const double&);
  friend inline bool operator!=(const double&,const cplx&);
  friend inline bool operator>(const cplx&,const cplx&);
  friend inline bool operator>(const cplx&,const double&);
  friend inline bool operator>(const double&,const cplx&);
  friend inline bool operator<(const cplx&,const cplx&);
  friend inline bool operator<(const cplx&,const double&);
  friend inline bool operator<(const double&,const cplx&);
  friend inline bool operator>=(const cplx&,const cplx&);
  friend inline bool operator>=(const cplx&,const double&);
  friend inline bool operator>=(const double&,const cplx&);
  friend inline bool operator<=(const cplx&,const cplx&);
  friend inline bool operator<=(const cplx&,const double&);
  friend inline bool operator<=(const double&,const cplx&);
  // here's the annoying thing:
  // Every function in class complex<double> that returns a
  // complex<double> causes ambiguities with function overloading
  // resolution because of the mix of types cplx and
  // complex<double> and double and int in math expressions.
  // So, although they are inherited, must redefine them
  // to return type cplx:
  // basic arithmetic
  inline cplx operator+() const;
  inline cplx operator+(const cplx&) const;
  inline cplx operator+(const double&) const;
  inline cplx operator+(const int&) const;
  inline friend cplx operator+(const double&, const cplx&);
  inline friend cplx operator+(const int&, const cplx&);
  inline cplx operator-() const;
  inline cplx operator-(const cplx&) const;
  inline cplx operator-(const double&) const;
  inline cplx operator-(const int&) const;
  inline friend cplx operator-(const double&, const cplx&);
  inline friend cplx operator-(const int&, const cplx&);
  inline cplx operator*(const cplx&) const;
  inline cplx operator*(const double&) const;
  inline cplx operator*(const int&) const;
  inline friend cplx operator*(const double&, const cplx&);
  inline friend cplx operator*(const int&, const cplx&);
  inline cplx operator/(const cplx&) const;
  inline cplx operator/(const double&) const;
  inline cplx operator/(const int&) const;
  inline friend cplx operator/(const double&, const cplx&);
  inline friend cplx operator/(const int&, const cplx&);
  // from <math.h>
  inline friend cplx sin(const cplx&);
  inline friend cplx sinh(const cplx&);
  inline friend cplx cos(const cplx&);
  inline friend cplx cosh(const cplx&);
  inline friend cplx tan(const cplx&);
  inline friend cplx tanh(const cplx&);
  inline friend cplx log10(const cplx&);
  inline friend cplx log(const cplx&);
  inline friend cplx sqrt(const cplx&);
  inline friend cplx exp(const cplx&);
  inline friend cplx pow(const cplx&, const cplx&);
  inline friend cplx pow(const cplx&, const double&);
  inline friend cplx pow(const cplx&, const int&);
  inline friend cplx pow(const double&, const cplx&);
  inline friend cplx pow(const int&, const cplx&);
  // complex versions of these are not in standard library
  // or they need to be redefined:
  // (frexp, modf, and fmod have not been dealt with)
  inline friend cplx fabs(const cplx&);
  inline friend cplx asin(const cplx&);
  inline friend cplx acos(const cplx&);
  inline friend cplx atan(const cplx&);
  inline friend cplx atan2(const cplx&, const cplx&);
  inline friend cplx ceil(const cplx&);
  inline friend cplx floor(const cplx&);
  inline friend cplx ldexp(const cplx&, const int&);
};

inline double RealPart(const cplx& c){
  return real(c);
}

inline double ImagPart(const cplx& c){
  return imag(c);
}

inline double RealPart(const double& r) {
  /***
   *  So the RealPart() statement can be used even with
   *  the double version of the code to be complexified.
   *  Most useful inside printf statements.
   ***/
  return r;
}

inline double ImagPart(const double& r) {
  return 0.;
}

inline bool operator==(const cplx& lhs, const cplx& rhs)
{
  return RealPart(lhs) == RealPart(rhs);
}

inline bool operator==(const cplx& lhs, const double& rhs)
{
  return RealPart(lhs) == rhs;
}

inline bool operator==(const double& lhs, const cplx& rhs)
{
  return lhs == RealPart(rhs);
}

inline bool operator!=(const cplx& lhs, const cplx& rhs)
{
  return RealPart(lhs) != RealPart(rhs);
}

inline bool operator!=(const cplx& lhs, const double& rhs)
{
  return RealPart(lhs) != rhs;
}

inline bool operator!=(const double& lhs, const cplx& rhs)
{
  return lhs != RealPart(rhs);
}

inline bool operator>(const cplx& lhs, const cplx& rhs)
{
  return RealPart(lhs) > RealPart(rhs);
}

inline bool operator>(const cplx& lhs, const double& rhs)
{
  return RealPart(lhs) > rhs;
}

inline bool operator>(const double& lhs, const cplx& rhs)
{
  return lhs > RealPart(rhs);
}

inline bool operator<(const cplx& lhs, const cplx& rhs)
{
  return RealPart(lhs) < RealPart(rhs);
}

inline bool operator<(const cplx& lhs, const double& rhs)
{
  return RealPart(lhs) < rhs;
}

inline bool operator<(const double& lhs, const cplx& rhs)
{
  return lhs < RealPart(rhs);
}

inline bool operator>=(const cplx& lhs, const cplx& rhs)
{
  return RealPart(lhs) >= RealPart(rhs);
}

inline bool operator>=(const cplx& lhs, const double& rhs)
{
  return RealPart(lhs) >= rhs;
}

inline bool operator>=(const double& lhs, const cplx& rhs)
{
  return lhs >= RealPart(rhs);
}

inline bool operator<=(const cplx& lhs, const cplx& rhs)
{
  return RealPart(lhs) <= RealPart(rhs);
}

inline bool operator<=(const cplx& lhs, const double& rhs)
{
  return RealPart(lhs) <= rhs;
}

inline bool operator<=(const double& lhs, const cplx& rhs)
{
  return lhs <= RealPart(rhs);
}

inline cplx cplx::operator+() const 
{
  return +std::complex<double>(*this);
}

inline cplx cplx::operator+(const cplx& z) const
{
  return std::complex<double>(*this)+std::complex<double>(z);
}

inline cplx cplx::operator+(const double& r) const
{
  return std::complex<double>(*this)+r;
}

inline cplx cplx::operator+(const int& i) const
{
  return std::complex<double>(*this)+double(i);
}

inline cplx operator+(const double& r, const cplx& z)
{
  return r+std::complex<double>(z);
}

inline cplx operator+(const int& i, const cplx& z)
{
  return double(i)+std::complex<double>(z);
}

inline cplx cplx::operator-() const
{
  return -std::complex<double>(*this);
}

inline cplx cplx::operator-(const cplx& z) const
{
  return std::complex<double>(*this)-std::complex<double>(z);
}

inline cplx cplx::operator-(const double& r) const
{
  return std::complex<double>(*this)-r;
}

inline cplx cplx::operator-(const int& i) const
{
  return std::complex<double>(*this)-double(i);
}

inline cplx operator-(const double& r, const cplx& z)
{
  return r-std::complex<double>(z);
}

inline cplx operator-(const int& i, const cplx& z)
{
  return double(i)-std::complex<double>(z);
}

inline cplx cplx::operator*(const cplx& z) const
{
  return std::complex<double>(*this)*std::complex<double>(z);
}

inline cplx cplx::operator*(const double& r) const
{
  return std::complex<double>(*this)*r;
}

inline cplx cplx::operator*(const int& i) const
{
  return std::complex<double>(*this)*double(i);
}

inline cplx operator*(const double& r, const cplx& z)
{
  return r*std::complex<double>(z);
}

inline cplx operator*(const int& i, const cplx& z)
{
  return double(i)*std::complex<double>(z);
}

inline cplx cplx::operator/(const cplx& z) const
{
  return std::complex<double>(*this)/std::complex<double>(z);
}

inline cplx cplx::operator/(const double& r) const
{
  return std::complex<double>(*this)/r;
}

inline cplx cplx::operator/(const int& i) const
{
  return std::complex<double>(*this)/double(i);
}

inline cplx operator/(const double& r, const cplx& z)
{
  return r/std::complex<double>(z);
}

inline cplx operator/(const int& i, const cplx& z)
{
  return double(i)/std::complex<double>(z);
}

inline cplx sin(const cplx& z)
{
  return sin(std::complex<double>(z));
}

inline cplx sinh(const cplx& z)
{
  return sinh(std::complex<double>(z));
}

inline cplx cos(const cplx& z)
{
  return cos(std::complex<double>(z));
}

inline cplx cosh(const cplx& z)
{
  return cosh(std::complex<double>(z));
}

#ifdef __GNUC__ // bug in gcc ?? get segv w/egcs-2.91.66 and 2.95.2
inline cplx tan(const cplx& z) 
{
  return sin(std::complex<double>(z))/cos(std::complex<double>(z));
}

inline cplx tanh(const cplx& z)
{
  return sinh(std::complex<double>(z))/cosh(std::complex<double>(z));
}

inline cplx log10(const cplx& z)
{
  return log(std::complex<double>(z))/log(10.);
}
#else
inline cplx tan(const cplx& z)
{
  return tan(std::complex<double>(z));
}

inline cplx tanh(const cplx& z)
{
  return tanh(std::complex<double>(z));
}

inline cplx log10(const cplx& z)
{
  return log10(std::complex<double>(z));
}
#endif

inline cplx log(const cplx& z)
{
  return log(std::complex<double>(z));
}

inline cplx sqrt(const cplx& z)
{
  return sqrt(std::complex<double>(z));
}

inline cplx exp(const cplx& z)
{
  return exp(std::complex<double>(z));
}

inline cplx pow(const cplx& a, const cplx& b)
{
  return pow(std::complex<double>(a),std::complex<double>(b));
}

inline cplx pow(const cplx& a, const double& b)
{
  return pow(std::complex<double>(a),b);
}

inline cplx pow(const cplx& a, const int& b)
{
  return pow(std::complex<double>(a),double(b));
}

inline cplx pow(const double& a, const cplx& b)
{
  return pow(a,std::complex<double>(b));
}

inline cplx pow(const int& a, const cplx& b)
{
  return pow(double(a),std::complex<double>(b));
}

inline cplx fabs(const cplx& z)
{
  return (RealPart(z)<0.0) ? -z:z;
}

#define surr_TEENY (1.e-24) /* machine zero compared to nominal magnitude of
			       the real part */

inline cplx asin(const cplx& z)
{
  // derivative trouble if ImagPart(z) = +/- 1.0
  return cplx( asin(RealPart(z)),
	       ImagPart(z)/sqrt(1.0-RealPart(z)*RealPart(z)+surr_TEENY));
}

inline cplx acos(const cplx& z)
{
  // derivative trouble if ImagPart(z) = +/- 1.0
  return cplx( acos(RealPart(z)),
	       -ImagPart(z)/sqrt(1.0-RealPart(z)*RealPart(z)+surr_TEENY));
}

#undef surr_TEENY

inline cplx atan(const cplx& z)
{
  return cplx( atan(RealPart(z)),
	       ImagPart(z)/(1.0+RealPart(z)*RealPart(z)));
}

inline cplx atan2(const cplx& z1, const cplx& z2)
{
  return cplx( atan2(RealPart(z1),RealPart(z2)),
	       (RealPart(z2)*ImagPart(z1)-RealPart(z1)*ImagPart(z2))
	       /(RealPart(z1)*RealPart(z1)+RealPart(z2)*RealPart(z2)));
}

inline cplx ceil(const cplx& z)
{
  return cplx( ceil(RealPart(z)), 0.0 );
}

inline cplx floor(const cplx& z)
{
  return cplx( floor(RealPart(z)), 0.0 );
}

inline cplx ldexp(const cplx& z, const int& i)
{
  return cplx( ldexp(RealPart(z),i), ldexp(ImagPart(z),i) );
}

#endif
