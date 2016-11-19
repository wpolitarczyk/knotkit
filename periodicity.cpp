#include <periodicity.h>

using polynomial = multivariate_laurentpoly;
using monomial = multivariate_laurent_monomial;

template<class T>
bool congruence_checker::reduce(const polynomial<T>& pol) const {
  polynomial res;
  for(typename map<monomial, T>::const_iter i = p.coeffs; i; i++) {
    int c = i.key().m[index] % (2 * prime);
    if(c < 0)
      c += (2 * prime);
    monomial mon = monomial(VARIABLE, index, c);
    res += polynomial(i.val(), mon);
  }
  if(verbose)
    std::cout << "res = " << res << "\n";
  return res == 0;
}

polynomial<Z> periodicity_checker::divide(const polynomial<Z> kh,
					  const polynomial<Z> lee) const {
  polynomial <Z> diff = kh - lee;
  polynomial <Z> quotient;
  while(diff != 0) {
    pair<monomial, Z> m = diff.head();
    if(m.first.m[1] == 1) {
      pair<monomial, Z> m1 = diff.tail();
      while(m1.first.m.card() == 1 && m1.first.m[2]) {
	quotient += polynomial(m1.second, m1.first);
	polynomial p = polynomial(m1.second, m1.first) * mul;
	diff -= p;
	m1 = diff.tail();
      }
    }
    quotient += polynomial(m.second, m.first);
    polynomial p = polynomial(m.second, m.first) * mul;
    diff -= p;
  }
  if(verbose) {
    std::cout << "Decomposition of the Khovanov polynomial = "
	      << lee << " + ("
	      << mul << ") * ("
	      << quotient << ")\n";
  }
  return quotient;
}

std::pair<polynomial<Z>, polynomial<Z>>
periodicity_checker::divide_p(const polynomial& q) const {
  polynomial quotient, remainder = 0;
  for(map<monomial, Z>::iter i = q.coeffs; i; i++) {
    std::tuple<Z,Z> div = i.val().divide_with_remainder(prime - 1);
    quotient += polynomial(std::get<0>(div), i.key());
    remainder += polynomial(std::get<1>(div), i.key());
  }
  return std::make_pair(quotient, remainder);
}

////// Functions that verify periodicity of knots and links

void periodicity_congruence(const multivariate_laurentpoly<Z>& lee,
			    const multivariate_laurentpoly<Z>& remainder,
			    const multivariate_laurentpoly<Z>& quotient,
			    int prime) {
  using polynomial = multivariate_laurentpoly<Z>;
  using monomial = multivariate_laurent_monomial;
  // prepare polynomials
  polynomial m = (polynomial(1) + polynomial(1, VARIABLE, 1) * polynomial(1, VARIABLE, 2, 2));
  polynomial r = lee + m * (remainder - quotient);
  // first check if quotient is zero
  if(quotient == 0) {
    // only one case to check
    if(verbose)
      std::cout << "All coefficients of the quotient "
	           "are smaller than "
		<< (prime -1) << "...\n";
    if(check_congruence((r - invert_variable(r,2)).evaluate(-1,1), prime))
      std::cout << knot << " may be " << prime << "-periodic...\n";
    else
      std::cout << knot << " is not " << prime << "-periodic..." << "\n";
  }
  // quotient not zero
  else {
    if(verbose)
      std::cout << "Decomposition of the quotient:" << "\n"
		<< remainder << " + "
		<< (prime - 1) << " * ("
		<< quotient << ")\n";
    //generate and check all cases
    std::vector<std::pair<monomial,Z>> v;
    Z number_of_cases = 1;
    for(map<monomial,Z>::const_iter i = quotient.coeffs; i; i++) {
      v.push_back(std::make_pair(i.key(), i.val()));
      number_of_cases *= (i.val() + 1);
      // std::cout << i.val() << "\n";
    }
    if(verbose) {
      std::cout << "There are "
		<< number_of_cases
		<< " cases to check...\n";
    }
    Z counter = 0;
    Z candidates = 0;
    // std::cout << "v.size() = " << v.size() << "\n";
    for(Z level = 0; level < v.size(); level += 1) {
      polynomial pol_temp = m * polynomial(prime, std::get<0>(v[level.get_ui()]));
      // std::cout << "level = " << level << " / "
      // 		<< (v.size() - 1) << "\n";
      int i = 0;
      if(level != Z(0))
	i++;
      for( ; Z(i) < std::get<1>(v[level.get_ui()]) + 1; i++) {
	// std::cout << "i = " << i << " / "
	// 	  << std::get<1>(v[level.get_ui()]) << "\n";
	polynomial p = r + polynomial(i) * pol_temp;
	if(level == 0) {
	  if(check_congruence((p - invert_variable(p, 2)).evaluate(-1,1), prime)) {
	    candidates += 1;
	    if(verbose)
	      std::cout << "Found a candidate..." << "\n";
	  }
	  counter += 1;
	  if(verbose) {
	      std::cout << counter
			<< " / "
			<< number_of_cases
			<< " cases checked...\n";
	  }
	} // level = 0
	else {
	  for(int level2 = 0; Z(level2) < level; level2++) {
	    // std::cout << "level2 = " << level2 << " / "
	    // 	      << level << "\n";
	    Z n_temp = std::get<1>(v[level2]);
	    polynomial mon_temp = m * polynomial(prime, std::get<0>(v[level2]));
	    for(int j = 0; Z(j) < n_temp + 1; j++) {
	      p += pol_temp;
	      // std::cout << "j = " << j << " / "
	      // 		<< n_temp << "\n";
	      if(check_congruence((p - invert_variable(p, 2)).evaluate(-1,1), prime)) {
		candidates += 1;
		if(verbose) {
		  std::cout << "Found a candidate..." << "\n";
		}
	      }
	      counter += 1;
	      if(verbose)
		std::cout << counter
			  << " / "
			  << number_of_cases
			  << " cases checked...\n"; 
	    } // loop over j
	  } // loop over level2
	} // level not zero
      } // loop over i
    } // loop over level
    if(candidates == 0) {
    std::cout << knot << " is not "
	      << prime << "-periodic...\n"; 
    }
    else {
      std::cout << knot << " may be "
		<< prime << "-periodic,\n"
		<< "found " << candidates
		<< " decompositions of khp.";
    }
  }
}
