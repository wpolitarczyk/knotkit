#ifndef _KNOTKIT_ALGEBRA_Q_H
#define _KNOTKIT_ALGEBRA_Q_H
#include <gmpxx.h>
#include <iostream>
#include <memory>
#include <tuple>
#include <cassert>

class Q
{
 public:
  using linear_combination = ::linear_combination<Q>;
  using linear_combination_const_iter = ::linear_combination_const_iter<Q>;

 private:
  std::shared_ptr<mpq_class> impl;
  void write_state() const {
    std::cout << "I store the following value " << *this << "\n";
    std::cout << "Number of objects pointing to the same value " << impl.use_count() << "\n";
    /* std::cout << "I point to " << impl.get() << "\n"; */
    /* std::cout << "My size " << sizeof(*this) << "\n"; */
    /* std::cout << "Size of std::shared_ptr<mpz_class> " << sizeof(impl) << "\n"; */
    /* std::cout << "Size of mpz_class " << sizeof(*impl) << "\n"; */
  }
  
 public:
  Q() : impl(new mpq_class()) {
#ifdef DEBUG_Q
    std::cout << "Q()" << "\n";
    write_state();
#endif
  }
  Q(mpq_t q) : impl(new mpq_class(q)) {
#ifdef DEBUG_Q
    std::cout << "Q(mpq_t)" << "\n";
    write_state();
#endif
  }
  Q(mpq_class q) : impl(new mpq_class(q)) {
#ifdef DEBUG_Q
    std::cout << "Q(mpq_class)" << "\n";
    write_state();
#endif
  }
  Q(int init) : impl(new mpq_class (init)) {
#ifdef DEBUG_Q
    std::cout << "Q(int)" << "\n";
    write_state();
#endif
  }
  Q(const Q& q) : impl(q.impl) {
#ifdef DEBUG_Q
    std::cout << "Q(const Q&)" << "\n";
    write_state();
#endif
  }
  Q(Q&& q) : impl(std::move (q.impl)) {
    q.impl = std::make_shared<mpq_class>(0);
#ifdef DEBUG_Q
    std::cout << "Q(Q&& q)" << "\n";
    write_state();
#endif
  }
  Q(copy, const Q &q) : impl(new mpq_class (*q.impl.get())) {
#ifdef DEBUG_Q
    std::cout << "Q(copy, const Q&)" << "\n";
    write_state();
#endif
  }
  //Q(reader &r) : impl(new Q_impl (r)) { }
  ~Q() {
#ifdef DEBUG_Q
    std::cout << "~Q()" << "\n";
    write_state();
    if(get_count() == 1)
      std::cout << "Destroying..." << "\n";
#endif
  }
  
  Q &operator = (const Q &q) {
    impl = q.impl;
    return *this;
  }
  Q& operator = (Q&& q) {
    impl = std::move(q.impl);
    q.impl = std::make_shared<mpq_class>(0);
    return *this;
  }
  Q &operator = (int x) {
    impl = std::make_shared<mpq_class>(0);
    return *this;
  }
  
  bool operator == (const Q &q) const {
    return *impl.get() == *q.impl.get();
  }
  bool operator == (const int r) const {
    return *impl.get() == mpq_class(r);
  }
  bool operator != (const Q& q) const {
    return ! operator == (q);
  }
  bool operator != (const int r) const {
    return !operator == (r);
  }
  
  bool operator < (const Q &q) const {
    return *impl.get() < *q.impl.get();
  }
  // bool operator > (const Q& q) const {
  //   return *impl.get() > *q.impl.get();
  // }
  
  bool is_unit () const
  {
    return *this != 0;
  }
  
  Q operator - () const
  {
    return Q(-*impl.get());
  }

  Q recip () const
  {
    mpq_t q;
    mpq_init(q);
    mpq_inv(q, impl.get()->get_mpq_t());
    return Q(q);
  }
  
  Q operator + (const Q& q) const
  {
    return Q(*impl.get() + *q.impl.get());
  }
  
  Q operator - (const Q& q) const
  {
    return Q(*impl.get() - *q.impl.get());
  }
  
  Q operator * (const Q& q) const
  {
    return Q(*impl.get() * (*q.impl.get()));
  }
  
  Q operator / (const Q& q) const
  {
    return Q(*impl.get() / *q.impl.get());
  }
  
  Q &muladdeq (const Q& q1, const Q& q2)
  {
    return operator += (q1 * q2);
  }
  
  Q &operator += (const Q& q)
  {
    return *this = *this + q;
  }
  
  Q &operator -= (const Q& q)
  {
    return *this = *this - q;
  }
  
  Q &operator *= (const Q &q)
  {
    return *this = *this * q;
  }
  
  Q &operator /= (const Q &q)
  {
    assert (q != 0);
    return *this = *this / q;
  }
  
  bool divides (const Q &num) const
  {
    return *this != 0 || num == 0;
  }
  
  bool operator | (const Q &num) const { return divides (num); }
  
  Q div (const Q &d) const { return operator / (d); }
  
  std::tuple<Q, Q, Q> extended_gcd (const Q &q) const
  {
    if (*this != 0)
      return std::tuple<Q, Q, Q> (*this, 1, 0);
    else
      return std::tuple<Q, Q, Q> (q, 0, 1);
  }
  
  Q gcd (const Q &q) const
  {
    assert (*this != 0 || q != 0);
    return 1;
  }
  friend std::ostream& operator << (std::ostream& os, const Q& q) {
    return os << *q.impl.get();
  }
  static void show_ring () { printf ("Q"); }
  void show_self () const { std::cout << *this; }
  void display_self () const { std::cout << *this << "\n"; }
  // void write_self (writer &w) const { write (w, *impl); }
  int get_count() const {
    return impl.use_count();
  }
  Z get_num() const {
    return Z(impl.get()->get_num());
  }
  Z get_den() const {
    return Z(impl.get()->get_den());
  }
};
#endif // _KNOTKIT_ALGEBRA_Q_H
