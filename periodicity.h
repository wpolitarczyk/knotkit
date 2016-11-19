#ifndef _KNOTKIT_PERIODICITY_H
#define _KNOTKIT_PERIODICITY_H

#include <knotkit.h>
#include <sstream>
#include <string>
#include <vector>
#include <utility>
#include <list>
#include <algorithm>

extern bool verbose;
extern const char* knot;

extern std::string periodicity_test;

const std::vector<int> primes_list = {5, 7, 11, 13, 17, 19};

enum class Test_type {
  Przytycki_criterion,
  BKP_criterion,
  all
};

template<class T>
class periodic_congruence_checker {
public:
  using polynomial = multivariate_laurentpoly<T>;
  using monomial = multivariate_laurent_monomial;
  
  int prime;
  unsigned index;

  polynomial prepare_polynomial(const polynomial& pol) const {
    return (pol - invert_variable(pol, index));
  }

  bool reduce(const polynomial& pol) const;

public:
  periodic_congruence_checker(int pprime = 5,
			      unsigned ind = 2) :
    prime(pprime),
    index(ind)
    {}

  ~periodic_congruence_checker() {};

  bool operator() (const polynomial& pol) const {
    return reduce(prepare_polynomial(pol));
  }
};

template<class T> 
bool periodic_congruence_checker<T>::reduce(const multivariate_laurentpoly<T>& pol) const {
  polynomial res;
  for(typename map<monomial, T>::const_iter i = pol.coeffs; i; i++) {
    int c = i.key().m[index] % (2 * prime);
    if(c < 0)
      c += (2 * prime);
    monomial mon = monomial(VARIABLE, index, c);
    res += polynomial(i.val(), mon);
  }
  if(verbose)
    std::cout << "reduced = " << res << "\n";
  return res == 0;
}

template<unsigned p>
class Przytycki_periodicity_checker {
  using polynomial = multivariate_laurentpoly<Zp<p>>;
  using monomial = multivariate_laurent_monomial;

  periodic_congruence_checker<Zp<p>> cong_checker;

 public:
  Przytycki_periodicity_checker() : cong_checker(p) {}
    
  ~Przytycki_periodicity_checker() {}
  
  std::string operator() (const polynomial& pol) const {
    std::ostringstream out;
    out << knot << ": period = " << p
	<< ": "
	<< (cong_checker(pol) ? "Maybe" : "No");
    return out.str();
  }
};

class Kh_periodicity_checker {
  using polynomial = multivariate_laurentpoly<Z>;
  using monomial = multivariate_laurent_monomial;

  unsigned ev_index;
  unsigned index;

  polynomial khp, leep, quot;
  polynomial mul;

  void compute_knot_polynomials(knot_diagram& kd);
  void compute_quot();
  std::pair<polynomial, polynomial> compute_quotient_and_remainder(const polynomial& p,
								   int period) const;
  std::list<polynomial> generate_candidates(const polynomial& q) const;
  bool check(const polynomial& q, const polynomial& r, int period) const;

 public:
  Kh_periodicity_checker(knot_diagram& kd) {
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
};

#endif // _KNOTKIT_PERIODICITY_H
