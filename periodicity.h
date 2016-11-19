#ifndef _KNOTKIT_PERIODICITY_H
#define _KNOTKIT_PERIODICITY_H

#include <algebra/algebra.h>
#include <vector>

extern bool verbose;
extern const char* knot;

int period = 5;
std::vector<int> primes = {2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31};

template<class T>
class congruence_checker {
  using polynomial = multivariate_laurentpoly<T>;
  using monomial = multivariate_laurent_monomial;

  int prime;
  unsigned index;
  unsigend ev_index;

  polynomial prepare_polynomial(const polynomial& pol) const {
    return (pol - invert_variable(pol, index)).evaluate(-1, ev_index);
  }

  bool reduce (const polynomial& pol) const;

public:
  congruence_checker(int pprime = 5,
		     unsigned ind = 2,
		     unsigned ev_ind = 1) :
    prime(pprime), index(ind), ev_index(ev_ind) {}

  congruence_checker(const congruence_checker& cc) =default;
  congruence_checker(congruence_checker&& cc) =default;
  
  ~congruence_checker() {};

  congruence_checker& operator= (const congruence_checker&& cc) =default;
  congruence_checker& operator= (congruence_checker&& cc) =default;

  bool operator () (const polynomial& pol) const {
    return reduce(prepare_polynomial(pol));
  }
};

enum class Test_type {
  przytycki_criterion,
  BKP_criterion,
  all
};

class periodicity_checker {
  using polynomial = multivariate_laurentpoly;
  using monomial = multivariate_laurent_monomial;

  const unsigned ev_index = 1;
  const unsigned index = 2;

  int prime;

  polynomial<Z> khp;
  polynomial<Z> leep;

  polynomial<Z> quotient;
  polynomial<Z> remainder;
  
  polynomial<Z> mul;

  Test_type type;

  void write_test_result(bool result) const {
    std::cout << knot;
    if(result)
      std::cout << " may be ";
    else
      std::cout << " is not ";
    std::cout << prime << "-periodic.\n";

  }

  polynomial<Z> divide(const polynomial<Z>& kh, const polynomial<Z>& lee) const;

  std::pair<polynomial<Z>, polynomial<Z>>
  divide_p(const polynomial<Z>& q) const;
  
public:
  periodicity_checker(polynomial<Z> kh,
		      polynomial<Z> lee,
		      int p = 5,
		      Test_type t = Test_type::all) :
    khp(kh),
    leep(lee),
    prime(p),
    type(t) {
    // check if prime is in primes

    // compute mul
    mul = polynomial(Z(1)) + polynomial(1, VARIABLE, ev_index) * polynomial(1, VARIABLE, index, 2);

    // compute quotient and remainder
    std::pair<polynomial<Z>, polynomial<Z>>
      p = divide_p(divide(khp, leep));
    quotient = std::get<0>(p);
    remainder = std::get<1>(p);
  }

  periodicity_checker(const periodicity_checker& pc) =default;
  periodicity_checker(periodicity_checker&& pc) =default;

  periodicity_checker& operator= (const periodicity_checker& pc) =default;
  periodicity_checker& operator= (periodicity_checker&& pc) =default;

  void operator () () const {
    if(type == Test_type::przytycki_criterion ||
       type == Test_type::all) {
      std::cout << "Przytycki's test: ";
      write_test_result(przytycki_test());
    }
    if(type == Test_type::BKP_criterion ||
       type == Test_type::all) {
      std::cout << "BKP test: ";
      write_test_result(BKP_test);
    }
  }

  bool przytycki_test() const {
    return true;
  }

  bool BKP_test() const {
    return false;
  }
}
#endif // _KNOTKIT_PERIODICITY_H
