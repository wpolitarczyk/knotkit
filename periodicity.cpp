#include <periodicity.h>
#include <simplify_chain_complex.h>

bool Przytycki_periodicity_checker::check(int period) const {
  switch(period) {
  case 5: {
    periodic_congruence_checker<Zp<5>> pcc(5);
    return pcc(jones_pol);
  }
  case 7: {
    periodic_congruence_checker<Zp<7>> pcc(7);
    return pcc(jones_pol);
  }
  case 11: {
    periodic_congruence_checker<Zp<11>> pcc(11);
    return pcc(jones_pol);
  }
  case 13: {
    periodic_congruence_checker<Zp<13>> pcc(13);
    return pcc(jones_pol);
  }
  case 17: {
    periodic_congruence_checker<Zp<17>> pcc(17);
    return pcc(jones_pol);
  }
  case 19: {
    periodic_congruence_checker<Zp<19>> pcc(19);
    return pcc(jones_pol);
  }
  }
  return false;
}

std::string Przytycki_periodicity_checker::operator () (int period) const {
  std::ostringstream res;
  res << knot << ": period = " << period << ": "
      << (check(period) ? "Maybe" : "No");
  return res.str();
}

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
	if(diff != 0)
	  m1 = diff.tail();
	else break;
      }
      if(diff != 0)
	m = diff.head();
      else
	break;
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

std::map<multivariate_laurentpoly<Z>, std::pair<Z,Z>>
Kh_periodicity_checker::compute_bounds(const polynomial& p, int period) const {
  std::map<polynomial, std::pair<Z, Z>> bounds;
  periodic_congruence_checker<Z> pcc(period);
  for(map<monomial, Z>::const_iter i = p.coeffs; i; ++i) {
    monomial mon;
    int exp = 0;
    if(i.key().m % ev_index) {
      exp = i.key().m[ev_index];
      for(map<unsigned, int>::const_iter j = i.key().m; j; ++j) {
	if(j.key() != ev_index) {
	  int v = j.val() % (2 * period);
	  if(v < 0) v += (2 * period);
	  mon *= monomial(VARIABLE, j.key(), v);
	}
      }
    }
    else {
      for(map<unsigned, int>::const_iter j = i.key().m; j; ++j) {
	int v = j.val() % (2 + period);
	if (v < 0) v += (2 * period);
	mon *= monomial(VARIABLE, j.key(), v);
      }
    }
    // std::cout << polynomial(i.val() * pow(-1, exp), mon) << "\n";
    Z v_temp = i.val() * pow(-1, exp);
    polynomial p_temp = (polynomial(1, mon) * mul).evaluate(-1, ev_index);
    p_temp = pcc.reduce(p_temp - invert_variable(p_temp, index));
    if(v_temp >= 0)
      bounds[p_temp].second += (v_temp * period);
    else
      bounds[p_temp].first += (v_temp * period);
  }
  
  // for(std::map<polynomial, std::pair<Z,Z>>::iterator mi = bounds.begin(); mi != bounds.end(); ++mi) {
  //   std::cout << "Monomial: " << mi->first << "\n";
  //   std::cout << "Max: " << std::get<1>(mi->second)
  // 	      << ", Min: " << std::get<0>(mi->second) << "\n";
  // }
  return bounds;
}

std::vector<multivariate_laurentpoly<Z>>
Kh_periodicity_checker::compute_basis_polynomials(int period) const {
  std::vector<polynomial> res;
  periodic_congruence_checker<Z> pcc(period);
  for(int i = 1; i < period; i += 2) {
    res.push_back(pcc.reduce(get_basis_polynomial(i)));
  }
  return res;
}

multivariate_laurentpoly<Z> Kh_periodicity_checker::get_basis_polynomial(monomial mon) const {
  return (polynomial(Z(1), mon) * mul).evaluate(-1, ev_index) -
    invert_variable((polynomial(Z(1), mon) * mul).evaluate(-1, ev_index), index);
}

bool Kh_periodicity_checker::check(const polynomial& q,
				   const polynomial& r,
				   int period) const {
  periodic_congruence_checker<Z> pcc(period);
  polynomial t = (leep + mul * (r - q)).evaluate(-1, ev_index);
  t = pcc.reduce(t - invert_variable(t, index));
  if(pcc(t)) {
    return true;
  }
  else if(q == 0)
    return false;
  // std::cout << t << std::endl;
  // std::cout << q << "\n";
  std::map<polynomial, std::pair<Z,Z>> bounds = compute_bounds(q, period);
  for(std::map<polynomial, std::pair<Z,Z>>::iterator it = bounds.begin();
      it != bounds.end(); ++it) {
    polynomial mon = it->first;
  }
  std::vector<polynomial> basis_polynomials = compute_basis_polynomials(period);
  polynomial p = pcc.reduce(get_basis_polynomial(2 * period - 1));
  for(Z i = bounds[p].first; i <= bounds[p].second; i += 5) {
    polynomial p_temp = t + polynomial(i, VARIABLE, index, 0) * p;
    // std::cout << "i = " << i << "\n";
    // std::cout << "p_temp = " << p_temp << "\n";
    if(p_temp == 0)
      return true;
    for(std::vector<polynomial>::iterator it = basis_polynomials.begin(); it != basis_polynomials.end(); ++it) {
      pair<monomial, Z> m = p_temp.coeffs.head();
      monomial mon = m.first;
      Z c = m.second;
      polynomial pp = pcc.reduce(get_basis_polynomial(mon));
      if(pp == *it) {
	if(c < bounds[pp].first || c > bounds[pp].second)
	  break;
	else {
	  // std::cout << "pp = " << pp << "\n";
	  p_temp -= polynomial(c, VARIABLE, index, 0) * pp;
	  // std::cout << "p_temp = " << p_temp << "\n";
	  if(p_temp == 0)
	    return true;
	}
      }
    }
  }
  
  return false;
}

std::string Kh_periodicity_checker::operator () (int period) const {
  std::ostringstream out;
  // first check Przytycki's criterion
  Przytycki_periodicity_checker P_pc(evaluate_with_copy<Z>(khp, -1, ev_index));
  if(!P_pc.check(period)) {
    out << knot_name << ": period =  " << period << ": No (Przytycki's criterion).";
  }
  else {
    std::pair<polynomial, polynomial> q_r = compute_quotient_and_remainder(quot, period);
    bool res = check(std::get<0>(q_r), std::get<1>(q_r), period);
    out << knot_name << ": period = " << period << ": "
	<< (res ? "Maybe" : "No");
  }
  return out.str();
}
