#ifndef KNOTKIT_ALGEBRA_Z2_H
#define KNOTKIT_ALGEBRA_Z2_H
#include <algebra/Z.h>
#include <iostream>
#include <tuple>
#include <cassert>

class Z2
{
 public:
  using linear_combination = ::linear_combination<Z2>;
  using linear_combination_const_iter = ::linear_combination_const_iter<Z2>;
  
 private:
  bool v;
  
 public:
  Z2 () : v(0) { }
  Z2 (int x) : v((bool)(x & 1)) { }
  Z2 (unsigned x) : v((bool)(x & 1)) { }
  Z2 (bool v_) : v(v_) { }
  Z2 (const Z2& x) : v(x.v) { }
  Z2 (Z2&& x) : v(std::move(x.v)) {
    x.v = 0;
  }
  Z2 (reader &r) { v = r.read_bool (); }
  Z2 (const Z& z) : v((z % 2).get_ui()) {}
  ~Z2 () { }
  
  Z2& operator = (const Z2& x) { v = x.v; return *this; }
  Z2& operator = (Z2&& x) { v = x.v; x.v = 0; return *this; }
  Z2& operator = (int x) { v = (bool)(x & 1); return *this; }
  
  bool operator == (const Z2& x) const { return v == x.v; }
  bool operator != (const Z2& x) const { return v != x.v; }
  
  bool operator == (int x) const { return v == (bool)(x & 1); }
  bool operator != (int x) const { return !operator == (x); }
  
  bool is_unit () const { return v; }
  
  Z2 operator + (const Z2& x) const { return Z2 (v ^ x.v); }
  Z2 operator - (const Z2& x) const { return Z2 (v ^ x.v); }
  Z2 operator - () const { return *this; }
  Z2 operator * (const Z2& x) const { return Z2 (v & x.v); }
  Z2 operator / (const Z2& x) const { assert (x.v); return *this; }
  
  Z2 recip () const { assert (v); return *this; }
  
  Z2& operator += (const Z2& x) { v ^= x.v; return *this; }
  Z2& operator -= (const Z2& x) { v ^= x.v; return *this; }
  Z2& operator *= (const Z2& x) { v &= x.v; return *this; }
  Z2& operator /= (const Z2& x) { assert (x.v); return *this; }
  
  // *this += z1*z2
  Z2& muladdeq (const Z2& z1, const Z2& z2)
  {
    return operator += (z1 * z2);
  }
  
  bool divides (Z2 x) const { return v || !x.v; }
  bool operator | (const Z2 x) const { return divides (x); }
  
  Z2 div (Z2 x) const
  {
    assert (x.divides (*this));
    return *this;
  }
  
  Z2 gcd (Z2 x) const
  {
    assert (v && x.v);
    return Z2 (1);
  }
  
  std::tuple<Z2, Z2, Z2> extended_gcd (Z2 x) const
  {
    if (v)
      return std::make_tuple (Z2 (1), Z2 (1), Z2 (0));
    else
      return std::make_tuple (Z2 (1), Z2 (0), Z2 (1));
  }

  friend std::ostream& operator << (std::ostream& os, const Z2& x) {
    return os << (int)x.v;
  }
  static void show_ring () { printf ("Z2"); }
  void write_self (writer &w) const { w.write_bool (v); }
  void show_self () const { std::cout << *this; }
  void display_self () const { std::cout << *this << "\n"; }
};
#endif //KNOTKIT_ALGEBRA_Z2_H
