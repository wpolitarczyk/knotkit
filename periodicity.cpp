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
}

std::string Przytycki_periodicity_checker::operator () (int period) const {
  std::ostringstream res;
  res << knot << ": period = " << period << ": "
      << (check(period) ? "Maybe" : "No");
  return res.str();
}

template<class T>
polynomial_iterator<T>::polynomial_iterator(const multivariate_laurentpoly<T>& pol,
					 start_pos sp) {
  for(typename map<monomial, T>::const_iter i = pol.coeffs; i; ++i) {
    monomials.push_back(i.key());
    bounds.push_back(i.val());
    current_pos.push_back(Z(0));
  }

  if(sp == start_pos::begin) {
    level = 0;
  }
  else {
    level = bounds.size();
  }
#ifndef NDEBUG
  check_current_pos();
#endif
}

#ifndef NDEBUG
template<class T>
void polynomial_iterator<T>::check_current_pos() {
  assert(bounds.size() == monomials.size());
  assert(bounds.size() == current_pos.size());
  assert(level <= current_pos.size());
  for(unsigned i = 0; i < current_pos.size(); i++) {
    if(i < level) {
      assert((current_pos[i] <= bounds[i]) &&
	     Z(0) <= (current_pos[i]));
    }
    else if(i == level) {
      if(level > 0)
	assert(current_pos[i] <= bounds[i] && Z(0) < current_pos[i]);
      else
	assert((current_pos[i] <= bounds[i]) && Z(0) <= (current_pos[i]));	
    }
    else
      assert(current_pos[i] == Z(0));
  }
}
#endif // NDEBUG

template<class T>
polynomial_iterator<T>&
polynomial_iterator<T>::operator ++ () {
#ifndef NDEBUG
  check_current_pos();
#endif
  if(level == monomials.size())
    return *this;

  unsigned i = 0;
  while(i <= level) {
#ifndef NDEBUG
    check_current_pos();
#endif
    if(current_pos[i] < bounds[i]) {
      current_pos[i] += 1;
      break;
    }
    else {
      if(i == level) {
	if(level < monomials.size() - 1) {
	  current_pos[i] = 0;
	  current_pos[i+1] += 1;
	  level++;
	  break;
	}
	else {
	  level++;
	  break;
	}
      }
      else {
	current_pos[i] = 0;
      }
    }
    i++;
  }
  return *this;
}

template<class T>
multivariate_laurentpoly<T> polynomial_iterator<T>::operator *() const {
  polynomial res;
  for(unsigned i = 0; i <= level; i++) {
    res += polynomial(current_pos[i], monomials[i]);
  }
  return res;
}

template<class T>
std::string polynomial_iterator<T>::write_self() const {
  std::ostringstream res;
  res << "level = " << level << std::endl
      << "monomials:" << std::endl;
  for(auto& mon : monomials)
    res << mon << std::endl;
  res << "bounds: " << std::endl;
  for(auto& b : bounds)
    res << b << std::endl;
  res << "current_pos: " << std::endl;
  for(auto& pos : current_pos)
    res << pos << std::endl;
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

bool Kh_periodicity_checker::check(const polynomial& q,
				   const polynomial& r,
				   int period) const {
  periodic_congruence_checker<Z> pcc(period);
  polynomial t = leep + mul * r;
  if(q == 0) {
    return pcc(t.evaluate(-1,1));
  }
  if(verbose)
    std::cout << "Checking congruences..."; 
  polynomial_iterator<Z> pi(q);
  polynomial_iterator<Z> pi_end(q, polynomial_iterator<Z>::start_pos::end);
  Z count = pi.get_count();
  if(verbose)
    std::cout << count << " candidates..." << std::endl;
  Z step = 1;
  if(count >= 32)
    step = count / 32;
  Z c = 0;
  while(pi != pi_end) {
    //std::cout << "pi: " << std::endl << pi;
    polynomial temp = t + polynomial(period - 1) * mul * (*pi);
    if(pcc(temp.evaluate(-1,1))) {
      std::cout << "Candidates:" << std::endl
  		<< "EKhP_1 = " << temp << std::endl
  		<< "EKhP_" << (period - 1) << " = "
  		<< (polynomial(period - 1) * mul *(r - *pi))
  		<< std::endl;
      return true;
    }
    ++pi;
    c += 1;
    if(verbose && c % step == 0)
      std::cout << c << "/" << count << "..." << std::endl;
  }
  return false;
}

std::string Kh_periodicity_checker::operator () (int period) const {
  std::ostringstream out;
  std::pair<polynomial, polynomial> q_r = compute_quotient_and_remainder(quot, period);
  bool res = check(std::get<0>(q_r), std::get<1>(q_r), period);
  out << knot << ": period = " << period << ": "
      << (res ? "Maybe" : "No");
  return out.str();
}
