#ifndef _KNOTKIT_PERIODICITY_H
#define _KNOTKIT_PERIODICITY_H

#include <knotkit.h>
#include <sstream>
#include <string>
#include <vector>
#include <utility>

extern bool verbose;
extern const char* knot;

extern std::string periodicity_test;

const std::vector<int> primes_list = {5, 7, 11, 13, 17, 19};

const unsigned eval_index = 1;
const unsigned invert_index = 2;

template<class T>
class periodic_congruence_checker {
protected:
  using polynomial = multivariate_laurentpoly<T>;
  using monomial = multivariate_laurent_monomial;
  
  int prime;
  unsigned index;

  polynomial prepare_polynomial(const polynomial& pol) const {
    polynomial inv = invert_variable(pol, index);
    return pol - inv;
  }

public:
  periodic_congruence_checker(int pprime = 5,
			      unsigned ind = invert_index) :
    prime(pprime),
    index(ind)
    {}

  virtual ~periodic_congruence_checker() {};

  const polynomial reduce(const polynomial& pol) const;
  
  bool operator() (const polynomial& pol) const {
    return reduce(prepare_polynomial(pol)) == 0;
    // return reduce(prepare_polynomial(pol)) == 0;
  }
};

template<class T> 
const multivariate_laurentpoly<T>
periodic_congruence_checker<T>::reduce(const multivariate_laurentpoly<T>& pol) const {
  polynomial res;
  for(typename map<monomial, T>::const_iter i = pol.coeffs; i; i++) {
    int c = i.key().m[index] % (2 * prime);
    if(c < 0)
      c += (2 * prime);
    monomial mon = monomial(VARIABLE, index, c);
    res += polynomial(i.val(), mon);
  }
  return res;
}

class Przytycki_periodicity_checker {
  using polynomial = multivariate_laurentpoly<Z>;
  using monomial = multivariate_laurent_monomial;

  polynomial jones_pol;

 public:
  Przytycki_periodicity_checker(polynomial j) : jones_pol(j) {}
    
  ~Przytycki_periodicity_checker() {}

  bool check(int period) const;
  
  std::string operator() (int period) const;
};

class Kh_periodicity_checker {
  using polynomial = multivariate_laurentpoly<Z>;
  using monomial = multivariate_laurent_monomial;

  unsigned ev_index;
  unsigned index;

  polynomial khp, leep, quot;
  polynomial mul;

  std::string knot_name;

  void compute_knot_polynomials(knot_diagram& kd);
  void compute_quot();
  std::pair<polynomial, polynomial> compute_quotient_and_remainder(const polynomial& p,
								   int period) const;
  std::map<polynomial, std::pair<Z,Z>>
  compute_bounds(const polynomial& p, int period) const;
  polynomial get_basis_polynomial(int exp) const {
    return (polynomial(1, VARIABLE, index, exp) * mul).evaluate(-1, ev_index) -
      invert_variable((polynomial(1, VARIABLE, index, exp) * mul).evaluate(-1, ev_index), index);
  }
  polynomial get_basis_polynomial(monomial mon) const;
  std::vector<polynomial> compute_basis_polynomials(int period) const;
  bool check(const polynomial& q, const polynomial& r, int period) const;

 public:
  Kh_periodicity_checker(knot_diagram& kd, std::string knot_n) :
    knot_name(knot_n) {
    ev_index = 1;
    index = 2;
    mul = polynomial(Z(1))
      + polynomial(1, VARIABLE, ev_index) *
      polynomial(1, VARIABLE, index, 2);

    compute_knot_polynomials(kd);
    compute_quot();
  }

  ~Kh_periodicity_checker() {}

  std::string operator () (int period) const;

  polynomial get_KhP() const {
    return khp;
  }

  polynomial get_LeeP() const {
    return leep;
  }
};

#endif // _KNOTKIT_PERIODICITY_H
