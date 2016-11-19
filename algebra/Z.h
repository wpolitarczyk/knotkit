#ifndef _KNOTKIT_ALGEBRA_Z_H
#define _KNOTKIT_ALGEBRA_Z_H
#include<gmpxx.h>
#include<tuple>
#include<memory>
#include<iostream>
#include<cassert>

class Z
{
 public:
  using linear_combination = ::linear_combination<Z>;
  using linear_combination_const_iter = ::linear_combination_const_iter<Z>;

 private:  
  std::shared_ptr<mpz_class> impl;
  void write_state() const {
    std::cout << "I store the following value " << *this << "\n";
    std::cout << "Number of objects pointing to the same value " << impl.use_count() << "\n";
    /* std::cout << "I point to " << impl.get() << "\n"; */
    /* std::cout << "My size " << sizeof(*this) << "\n"; */
    /* std::cout << "Size of std::shared_ptr<mpz_class> " << sizeof(impl) << "\n"; */
    /* std::cout << "Size of mpz_class " << sizeof(*impl) << "\n"; */
  }
  
 public:
  Z() : impl(new mpz_class) {
#ifdef DEBUG_Z
    std::cout << "Z()" << "\n";
    write_state();
#endif
  }
  Z(mpz_t z) : impl(new mpz_class(z)) {
#ifdef DEBUG_Z
    std::cout << "Z(mpz_t)" << "\n";
    write_state();
#endif
  }
  Z(mpz_class z) : impl(new mpz_class(z)) {
#ifdef DEBUG_Z
    std::cout << "Z(mpz_class)" << "\n";
    write_state();
#endif
  }
  Z(int init) : impl(new mpz_class(init)) {
#ifdef DEBUG_Z
    std::cout << "Z(int)" << "\n";
    write_state();
#endif
  }
  Z(const Z& z) : impl(z.impl) {
#ifdef DEBUG_Z
    std::cout << "Z(const Z& z)" << "\n";
    write_state();
#endif
  }
  Z(Z&& z) : impl(std::move(z.impl)) {
#ifdef DEBUG_Z
    std::cout << "Z(Z&& z)" << "\n";
    write_state();
#endif
    z.impl = std::make_shared<mpz_class>(mpz_class(0));
  }
  ~Z() {
#ifdef DEBUG_Z
    std::cout << "~Z()" << "\n";
    write_state();
    if(impl.use_count() == 1)
      std::cout << "Destroying..." << "\n";
#endif
  }

  Z& operator = (const Z& z) {
    impl = z.impl;
#ifdef DEBUG_Z
    std::cout << "Z& operator = (const Z&)" << "\n";
    write_state();
#endif
    return *this;
  }
  
  Z& operator = (const int x) {
    impl = std::make_shared<mpz_class>(mpz_class(x));
#ifdef DEBUG_Z
    std::cout << "Z& operator = (int)" << "\n";
    write_state();
#endif
    return *this;
  }
  Z& operator = (Z&& z) {
    impl = std::move(z.impl);
    z.impl = std::make_shared<mpz_class>(mpz_class(0));
#ifdef DEBUG_Z
    std::cout << "Z& operator = (Z&&)" << "\n";
    write_state();
#endif
    return *this;
  }

  bool operator == (const Z& z) const {
    return *impl.get() == *z.impl.get();
  }
  bool operator != (const Z& z) const {
    return !operator == (z);
  }

  bool operator == (const int y) const {
    return *impl.get() == mpz_class(y);
  }
  bool operator != (const int y) const {
    return ! operator == (y);
  }

  bool operator < (const Z& z) const {
    return *impl.get() < *z.impl.get();
  }

  bool is_unit () const
  {
    return *this == 1 || *this == -1;
  }

  Z operator + (const Z& z) const {
    return Z(*impl.get() + *z.impl.get());
  }

  Z operator - () const {
    return Z(-*impl.get());
  }

  Z operator - (const Z& z) const {
    return Z(*impl.get() - *z.impl.get());
  }

  Z operator * (const Z& z) const {
    return Z(*impl.get() * *z.impl.get());
  }

  Z operator / (const Z& z) const {
    if(z == 0)
      return *this;
    else {
      assert(z != 0);
      return Z(*impl.get() / *z.impl.get());
    }
  }

  Z recip () const {
    assert(is_unit());
    return *this;
  }

  Z& muladdeq(const Z& z1, const Z& z2) {
    return *this += z1 * z2;
  }

  Z& operator += (const Z& z) {
    *this = *this + z;
    return *this;
  }

  Z& operator -= (const Z& z) {
    *this = *this - z;
    return *this;
  }

  Z& operator *= (const Z& z) {
    *this = *this * z;
    return *this;
  }

  Z& operator /= (const Z& z) {
    *this = *this / z;
    return *this;
  }

  bool divides(const Z& num) const {
    return mpz_divisible_p(num.impl.get()->get_mpz_t(),impl.get()->get_mpz_t());
  }

  bool operator | (const Z& num) const {
    return divides(num);
  }

  Z divide_exact(const Z& denom) const {
    mpz_t q;
    mpz_init(q);
    mpz_divexact(q,impl.get()->get_mpz_t(), denom.impl.get()->get_mpz_t());
    return Z(q);
  }

  std::tuple<Z,Z> divide_with_remainder(const Z& denom) const {
    mpz_t q,r;
    mpz_init(q);
    mpz_init(r);
    mpz_tdiv_qr(q, r, impl.get()->get_mpz_t(), denom.impl.get()->get_mpz_t());
    return std::make_tuple(Z(q), Z(r));
  }

  Z gcd (const Z& z) const {
    mpz_t d;
    mpz_init(d);
    mpz_gcd(d, impl.get()->get_mpz_t(), z.impl.get()->get_mpz_t());
    return Z(d);
  }

  Z lcm (const Z& z) const {
    mpz_t m;
    mpz_init(m);
    mpz_lcm(m, impl.get()->get_mpz_t(), z.impl.get()->get_mpz_t());
    return Z(m);
  }

  std::tuple<Z, Z, Z> extended_gcd(const Z& z) const {
    mpz_t d,s,t;
    mpz_init(d);
    mpz_init(s);
    mpz_init(t);
    mpz_gcdext(d, s, t, impl.get()->get_mpz_t(), z.impl.get()->get_mpz_t());
    return std::make_tuple<Z, Z, Z>(Z(d), Z(s), Z(t));
  }
  
  static void show_ring () { printf ("Z"); }
  void show_self () const {
    std::cout << *this;
  }
  void display_self () const {
    std::cout << *this << "\n";
  }

  friend std::ostream& operator << (std::ostream& os, const Z& z) {
    return os << *z.impl;
  }
  int get_count() const {
    return impl.use_count();
  }
};

#endif // _KNOTKIT_ALGEBRA_Z_H
