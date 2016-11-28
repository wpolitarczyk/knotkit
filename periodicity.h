#ifndef _KNOTKIT_PERIODICITY_H
#define _KNOTKIT_PERIODICITY_H

#include <knotkit.h>
#include <sstream>
#include <string>
#include <vector>

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

  bool reduce(const polynomial& pol) const;

public:
  periodic_congruence_checker(int pprime = 5,
			      unsigned ind = invert_index) :
    prime(pprime),
    index(ind)
    {}

  virtual ~periodic_congruence_checker() {};

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
  // if(verbose)
  //   std::cout << "reduced = " << res << "\n";
  return res == 0;
}

class Przytycki_periodicity_checker {
  using polynomial = multivariate_laurentpoly<Z>;
  using monomial = multivariate_laurent_monomial;

  polynomial jones_pol;

  bool check(int period) const;

 public:
  Przytycki_periodicity_checker(polynomial j) : jones_pol(j) {}
    
  ~Przytycki_periodicity_checker() {}
  
  std::string operator() (int period) const;
};

template<class T>
class polynomial_iterator {
  using polynomial = multivariate_laurentpoly<T>;
  using monomial = multivariate_laurent_monomial;

  std::vector<monomial> monomials;
  std::vector<T> bounds;
  std::vector<T> current_pos;
  unsigned level;

  void check_current_pos();

public:
  enum class start_pos { begin, end };
  
  polynomial_iterator(const polynomial& init,
		      start_pos sp = start_pos::begin);
  polynomial_iterator(const polynomial_iterator& pi) =default;
  polynomial_iterator(polynomial_iterator&& pi) =default;

  ~polynomial_iterator() {}

  polynomial_iterator& operator = (const polynomial_iterator& pi) =default;
  polynomial_iterator& operator = (polynomial_iterator&& pi) =default;

  polynomial_iterator& operator ++();

  bool operator == (const polynomial_iterator& pi) const {
    if(level == monomials.size() || pi.level == pi.monomials.size()) {
      return level == pi.level &&
	monomials == pi.monomials &&
	bounds == pi.bounds;
    }
    else {
      return level == pi.level &&
	bounds == pi.bounds &&
	monomials == pi.monomials &&
	current_pos == pi.current_pos;
    }
  }
  bool operator != (const polynomial_iterator& pi) const {
    return !(*this == pi);
  }
  polynomial operator*() const;

  T get_count() const {
    Z res = 1;
    for(auto& v : bounds)
      res *= (v + 1);
    return res;
  }

  std::string write_self() const;
  friend inline std::ostream& operator << (std::ostream& os, const polynomial_iterator& pi) {
    return os << pi.write_self();
  }
};

template<class T>
std::ostream& operator << (std::ostream& os, const polynomial_iterator<T>& pi) {
  return os << *pi;
}

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
  // std::list<polynomial> generate_candidates(const polynomial& q) const;
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

  polynomial get_KhP() const {
    return khp;
  }

  polynomial get_LeeP() const {
    return leep;
  }
};

#endif // _KNOTKIT_PERIODICITY_H
