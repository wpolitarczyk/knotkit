#include <periodicity.h>
#include <simplify_chain_complex.h>

void Kh_periodicity_checker::compute_knot_polynomials(knot_diagram& kd) {

  unsigned m = kd.num_components ();
  if (m != 1) {
    std::cerr << "warning: this implementation of the criterion works for knots only...";
    exit (EXIT_FAILURE);
  }
      
  cube<Z2> c (kd, 0);
  ptr<const module<Z2> > C = c.khC;
      
  mod_map<Z2> d = c.compute_d (1, 0, 0, 0, 0);
  for (unsigned i = 1; i <= kd.n_crossings; i ++)
    d = d + c.H_i (i);
  assert (d.compose (d) == 0);

  // computing Khovanov homology
  if(verbose)
    std::cout << "Computing Khovanov homology" << std::endl;
  {
    chain_complex_simplifier<Z2> s (C, d, maybe<int>(1), maybe<int>(0));
    C = s.new_C;
    d = s.new_d;
    khp = C->free_poincare_polynomial();
    if(verbose)
      std::cout << "KhP = " << khp << "\n";
  }
  
  // computing Lee homolgy
  if(verbose)
    std::cout << "Computing Lee homology" << std::endl;
  {
    chain_complex_simplifier<Z2> s(C, d, maybe<int>(1), maybe<int>(2));
    C = s.new_C;
    d = s.new_d;
    leep = C->free_poincare_polynomial();
    if(d != 0) {
      std::cout << "For now, you can only use this criterion on Kh-thin knots." << std::endl;
      exit(EXIT_FAILURE);
    }
    if(verbose) {
      std::cout << "LeeP = " << leep << "\n";
    }
  }
}

void Kh_periodicity_checker::compute_quot() {
  polynomial diff = khp - leep;
  while(diff != 0) {
    pair<monomial, Z> m = diff.head();
    if(m.first.m[1] == 1) {
      pair<monomial, Z> m1 = diff.tail();
      while(m1.first.m.card() == 1 && m1.first.m[2]) {
	quot += polynomial(m1.second, m1.first);
	polynomial p = polynomial(m1.second, m1.first) * mul;
	diff -= p;
	m1 = diff.tail();
      }
    }
    quot += polynomial(m.second, m.first);
    polynomial p = polynomial(m.second, m.first) * mul;
    diff -= p;
  }
}

std::pair<multivariate_laurentpoly<Z>, multivariate_laurentpoly<Z>>
Kh_periodicity_checker::compute_quotient_and_remainder(const polynomial& quot,
						       int period) const {
  polynomial quotient, remainder;
  for(map<monomial, Z>::const_iter i = quot.coeffs; i; i++) {
    std::tuple<Z,Z> div = i.val().divide_with_remainder(period - 1);
    quotient += polynomial(std::get<0>(div), i.key());
    remainder += polynomial(std::get<1>(div), i.key());
  }
  if(verbose) {
    std::cout << "Decomposition of Khp = " << std::endl
	      << leep << " + ("
	      << mul << ") * ("
	      << remainder;
    if(quotient != 0) {
      std::cout << " + " << (period - 1)
		<< " * (" << quotient
		<< ")";
    }
    std::cout << ")" << std::endl;
   
  }
  return std::make_pair(quotient, remainder);
}

std::list<multivariate_laurentpoly<Z>>
Kh_periodicity_checker::generate_candidates(const polynomial &q) const {
  std::list<polynomial> result;
  Z size = 0;
  map<monomial, Z>::const_iter i = q.coeffs;

  for(int j = 0; Z(j) <= i.val(); j++) {
    result.push_back(polynomial(Z(j), i.key()));
    size += 1;
  }

  i++;
  
  for( ; i; i++) {
    for(int j = 1; Z(j) <= i.val(); j++) {
      std::list<polynomial> temp_list;
      for_each(result.begin(), result.end(),
	       [&result, &temp_list, j, i](const polynomial& p){
		 temp_list.push_back(p + polynomial(Z(j), i.key()));
	       });
      result.splice(result.end(), temp_list);
      size += 1;
    }
  }
  return result;
}

bool Kh_periodicity_checker::check(const polynomial& q,
				   const polynomial& r,
				   int period) const {
  periodic_congruence_checker<Z> pcc(period);
  polynomial t = leep + mul * r;
  if(q == 0) {
    return pcc(t.evaluate(-1,1));
  }
  else {
    // generate all polynomials
    if(verbose) {
      std::cout << "Generating all candidates..." << std::endl;
    }
    std::list<polynomial> candidates = generate_candidates(q);
    // Z i = 0;
    // for(auto& l : candidates) {
    //   i += 1;
    //   std::cout << i << ": " << l << std::endl;
    // }
    // and check each one of them
    if(verbose)
      std::cout << "Checking congruences..." << std::endl;
    polynomial m = mul;
    return any_of(candidates.begin(), candidates.end(),
		  [&pcc, &t, &m, &period](const polynomial& p){
		    return pcc((t + Z(period - 1) * m * p).evaluate(-1,1));
      });
  }
}

std::string Kh_periodicity_checker::operator () (int period) const {
  std::ostringstream out;
  std::pair<polynomial, polynomial> q_r = compute_quotient_and_remainder(quot, period);
  bool res = check(std::get<0>(q_r), std::get<1>(q_r), period);
  out << knot << ": period = " << period << ": "
      << (res ? "Maybe" : "No");
  return out.str();
}
