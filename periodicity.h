#ifndef _KNOTKIT_PERIODICITY_H
#define _KNOTKIT_PERIODICITY_H

#include <knotkit.h>
#include <sstream>
#include <string>
#include <vector>
#include <utility>
#include <tuple>

extern bool verbose;

constexpr std::array<int, 6> primes_list = {5, 7, 11, 13, 17, 19};

constexpr unsigned eval_index = 1;
constexpr unsigned invert_index = 2;

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
  periodic_congruence_checker(int p = 5,
			      unsigned ind = invert_index) :
    prime(p),
    index(ind)
    {}

  virtual ~periodic_congruence_checker() {};

  const polynomial reduce(const polynomial& pol) const;
  
  bool operator() (const polynomial& pol) const {
    return reduce(prepare_polynomial(pol)) == 0;
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
  std::string knot_name;

 public:
  Przytycki_periodicity_checker(polynomial j, std::string knot_n) :
    jones_pol(j), knot_name(knot_n) {}
    
  ~Przytycki_periodicity_checker() {}

  bool check(int period) const;
  
  std::string operator() (int period) const;
};

class Kh_bounds_iterator {
  using polynomial = multivariate_laurentpoly<Z>;
  using monomial = multivariate_laurent_monomial;
  using polynomial_tuple = std::vector<std::tuple<polynomial, polynomial, polynomial>>;
  using bounds_vector = std::map<multivariate_laurentpoly<Z>, std::pair<Z, Z>>;
  
  bounds_vector bv;
  int period;
  std::map<polynomial, Z> current_state;
  std::map<polynomial, std::pair<Z,Z>>::iterator level;

public:
  Kh_bounds_iterator(bounds_vector v, int p) :
    bv(v), period(p) {
    for(auto& v: bv) {
      current_state[v.first] = v.second.first;
    }
    level = bv.begin();
  }
  ~Kh_bounds_iterator() {}

  bool advance();
  polynomial get_polynomial() const;
};

class Kh_periodicity_checker {
  using polynomial = multivariate_laurentpoly<Z>;
  using monomial = multivariate_laurent_monomial;
  using polynomial_tuple = std::vector<std::tuple<polynomial, polynomial, polynomial>>;
  using bounds_vector = std::map<multivariate_laurentpoly<Z>, std::pair<Z, Z>>;

  unsigned ev_index;
  unsigned index;

  enum class Test_Result { MAYBE, NO, NO_NONTRIVIAL_DECOMP };

  polynomial khp, leep;
  std::vector<polynomial> quot, mul, quotients, remainders;

  std::string knot_name;
  std::string field;

  template<typename R>
  std::vector<polynomial> compute_knot_polynomials(knot_diagram& kd);
  void compute_quot(const std::vector<polynomial>& lee_ss_polynomials);
  polynomial_tuple
  compute_quotient_and_remainder(const std::vector<polynomial>& p, int period) const;
  bounds_vector
  compute_bounds(const polynomial_tuple& p, int period) const;
  Test_Result check(const polynomial_tuple& polynomials, int period) const;

 public:
  Kh_periodicity_checker(knot_diagram& kd, std::string knot_n, std::string f);

  ~Kh_periodicity_checker() {}

  std::string operator () (int period) const;

  polynomial get_KhP() const {
    return khp;
  }

  polynomial get_LeeP() const {
    return leep;
  }
};

void check_periodicity(knot_diagram& kd, const std::string knot_name,
		       int period = 5, const std::string field = "Z2");

#endif // _KNOTKIT_PERIODICITY_H
