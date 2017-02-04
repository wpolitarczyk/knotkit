#include <periodicity.h>
#include <simplify_chain_complex.h>
#include <algorithm>
#include <utility>

using polynomial_tuple = std::vector<std::tuple<multivariate_laurentpoly<Z>, multivariate_laurentpoly<Z>, multivariate_laurentpoly<Z>>>;

using bounds_vector = std::map<multivariate_laurentpoly<Z>, std::pair<Z, Z>>;

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

bool Kh_bounds_iterator::advance() {
  if(level == bv.end())
    return false;
  for(auto bv_it = bv.begin(); bv_it != level; ++bv_it) {
    if(current_state[bv_it->first] < (bv_it->second).second) {
      current_state[bv_it->first] += period;
      for(auto bv_it_2 = bv.begin(); bv_it_2 != bv_it; ++bv_it_2) {
	current_state[bv_it_2->first] = bv_it_2->second.first;
      }
      return true;
    }
  }
  if(current_state[level->first] < bv[level->first].second) {
    current_state[level->first] += period;
    for(auto bv_it = bv.begin(); bv_it != level; ++bv_it) {
      current_state[bv_it->first] = bv_it->second.first;
    }
    return true;
  }
  ++level;
  if(level == bv.end())
    return false;
  current_state[level->first] += period;
  for(auto bv_it = bv.begin(); bv_it != level; ++bv_it) {
    current_state[bv_it->first] = bv_it->second.first;
  }
  
  return true;
}

multivariate_laurentpoly<Z> Kh_bounds_iterator::get_polynomial() const {
  polynomial p;
  for(auto& cs : current_state) {
    p += cs.second * cs.first;
  }
  return p;
}

std::vector<multivariate_laurentpoly<Z>>
Kh_periodicity_checker::compute_knot_polynomials(knot_diagram& kd) {
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
    std::cerr << "Computing Khovanov homology" << std::endl;
  std::vector<polynomial> lee_ss_polynomials;
  int k = 0;
  for(;;) {
    chain_complex_simplifier<Z2> s(C, d, maybe<int>(1), maybe<int>(2*k));
    C = s.new_C;
    d = s.new_d;
    lee_ss_polynomials.push_back(C->free_poincare_polynomial());
    if(k != 0)
      mul.push_back(polynomial(Z(1)) + polynomial(Z(1), VARIABLE, 1, 1) * polynomial(Z(1), VARIABLE, 2, 2 * k));
    if(d == 0)
      break;
    k++;
  }
  
  khp = *lee_ss_polynomials.begin();
  leep = *lee_ss_polynomials.rbegin();

  if(verbose) {
    std::cerr << "KhP = " << khp << "\n";
    std::cerr << "LeeP = " << leep << "\n";
  }
  // for(unsigned i = 0; i < lee_ss_polynomials.size(); ++i) {
  //   std::cerr << "lee_ss_polynomials[" << i << "]= "
  // 	      << lee_ss_polynomials[i] << "\n";
  //   std::cerr << "mul[" << i << "] = " << mul[i] << "\n";
  // }
  return lee_ss_polynomials;
}

void Kh_periodicity_checker::compute_quot(const std::vector<polynomial>& lee_ss_polynomials) {
  // quot.push_back(polynomial(Z(0)));
  for(unsigned i = 1; i < lee_ss_polynomials.size(); ++i) {
    polynomial diff = lee_ss_polynomials[i-1] - lee_ss_polynomials[i];
    polynomial q = 0;
    // std::cerr << "diff = " << diff << "\n";
    // std::cerr << "mul = " << mul[i-1] << "\n";
    while(diff != 0) {
      pair<monomial, Z> m = diff.head();
      if(m.first.m[1] == 1) {
    	pair<monomial, Z> m1 = diff.tail();
    	while(m1.first.m.card() == 1 && m1.first.m[2]) {
    	  q += polynomial(m1.second, m1.first);
    	  polynomial p = polynomial(m1.second, m1.first) * mul[i-1];
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
      q += polynomial(m.second, m.first);
      polynomial p = polynomial(m.second, m.first) * mul[i-1];
      diff -= p;
    }
    quot.push_back(q);
  }
  // for(unsigned i = 0; i < quot.size(); ++i) {
  //   std::cerr << "quot[" << i << "] = " << quot[i] << "\n";
  // }
}

polynomial_tuple
Kh_periodicity_checker::compute_quotient_and_remainder(const std::vector<polynomial>& quot, int period) const {
  polynomial_tuple decomposed_khp;
  for(unsigned i = 0; i < quot.size(); ++i) {
    polynomial quotient, remainder;
    for(map<monomial, Z>::const_iter j = quot[i].coeffs; j; j++) {
      std::tuple<Z,Z> div = j.val().divide_with_remainder(period - 1);
      quotient += polynomial(std::get<0>(div), j.key());
      remainder += polynomial(std::get<1>(div), j.key());
    }
    decomposed_khp.push_back(std::make_tuple(quotient, remainder, std::move(mul[i])));
  }
  if(verbose) {
    std::cerr << "Decomposition of Khp = " << std::endl
  	      << leep;
    for(auto& p: decomposed_khp) {
      polynomial quotient, remainder, mul;
      tie(quotient, remainder, mul) = p;
      std::cerr << " + (" << mul << ") * ("
		<< remainder;
      if(quotient != 0)
	std::cerr << " + " << (period - 1)
		  << "*(" << quotient << ")";
      std::cerr << ")";
    }
    std::cerr << "\n";
  }
  return decomposed_khp;
}

bounds_vector
Kh_periodicity_checker::compute_bounds(const polynomial_tuple& p_tuple, int period) const {
  periodic_congruence_checker<Z> pcc(period);
  bounds_vector bounds_v;
  for(auto& p: p_tuple) {
    polynomial quotient, remainder, mul;
    tie(quotient, remainder, mul) = p;
    for(map<monomial, Z>::const_iter i = quotient.coeffs; i; ++i) {
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
	  int v = j.val() % (2 * period);
	  if(v < 0) v += (2 * period);
	  mon *= monomial(VARIABLE, j.key(), v);
	}
      }
      // std::cerr << polynomial(i.val() * pow(-1,exp), mon) << "\n";
      Z v_temp = i.val() * pow(-1, exp);
      polynomial p_temp = (polynomial(1, mon) * mul).evaluate(-1, ev_index);
      p_temp = pcc.reduce(p_temp - invert_variable(p_temp, index));
      // std::cerr << "p_temp = " << p_temp << "\n";
      // std::cerr << "v_temp = " << v_temp << "\n";
      // std::cerr << "min_exp = " << min_exp << "\n";
      if(bounds_v.count(p_temp)) {
	if(v_temp >= 0)
	  bounds_v[p_temp].second += (v_temp * period);
	else
	  bounds_v[p_temp].first += (v_temp * period);
      }
      else if(bounds_v.count(p_temp)) {
	if(v_temp >= 0)
	  bounds_v[p_temp].first -= (v_temp * period);
	else
	  bounds_v[p_temp].second -= (v_temp * period);
      }
      else {
	bounds_v.emplace(p_temp,
			 std::make_pair<Z,Z>((v_temp < 0 ? (v_temp * period) : Z(0)), (v_temp >= 0 ? (v_temp * period) : Z(0))));
      }
    }
  }
  
  if(verbose) {
    for(auto& t: bounds_v) {
      Z neg, pos;
      tie(neg, pos) = t.second;
      std::cerr << "polynomial = " << t.first << "\n";
      std::cerr << "min = " << neg << ", max = " << pos << "\n";
    }
  }
  return bounds_v;
}

Test_Result Kh_periodicity_checker::check(const polynomial_tuple& polynomials,
				   int period) const {
  periodic_congruence_checker<Z> pcc(period);
  polynomial t = polynomial(COPY, leep);
  for(auto& p : polynomials) {
    polynomial quotient, remainder, mul;
    tie(quotient, remainder, mul) = p;
    t += mul * (remainder - quotient);
    //std::cerr << "t = " << t << "\n";
  }
  polynomial s = t.evaluate(-1, ev_index);
  s = pcc.reduce(s - invert_variable(s, index));
  if(pcc(s)) {
    return Test_Result::MAYBE;
  }
  else if(all_of(polynomials.begin(), polynomials.end(),
		 [](std::tuple<polynomial, polynomial, polynomial> t)
		 { return get<0>(t) == 0; }))
    return Test_Result::NO;
  bounds_vector bounds = compute_bounds(polynomials, period);

  if(verbose)
    std::cerr << "s = " << s << "\n";
  Kh_bounds_iterator Kh_b_it(bounds, period);
  if(verbose)
    std::cerr << "current_state = " << Kh_b_it.get_polynomial() << "\n";
  if(Kh_b_it.get_polynomial() == s)
    return Test_Result::MAYBE;
  while(Kh_b_it.advance()) {
    if(verbose)
      std::cerr << "current_state = " << Kh_b_it.get_polynomial() << "\n";
    if(s == Kh_b_it.get_polynomial())
      return Test_Result::MAYBE;
  }
  
  return Test_Result::NO_NONTRIVIAL_DECOMP;
}

std::string Kh_periodicity_checker::operator () (int period) const {
  std::ostringstream out;
  // first check Przytycki's criterion
  Przytycki_periodicity_checker P_pc(evaluate_with_copy<Z>(khp, -1, ev_index));
  if(!P_pc.check(period)) {
    out << knot_name << ": period =  " << period << ": No (Przytycki's criterion).";
  }
  else {
    auto q_r = compute_quotient_and_remainder(quot, period);
    Test_Result res = check(q_r, period);
    out << knot_name << ": period = " << period << ": "
  	<< (res == Test_Result::MAYBE ? "Maybe" :
	    (res == Test_Result::NO ? "No" : "No (Nontrivial decomposition)."));
  }
  return out.str();
}
