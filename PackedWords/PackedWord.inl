//
//  PackedWord.inl
//  sp4m
//
//  Created by Matthew Lambert on 10/13/22.
//

#ifndef packed_word_inl
#define packed_word_inl

/**
 * Construct PackedWord from word T.
 *
 * @param v Value to copy
 */
template <std::unsigned_integral T>
PackedWord<T>::PackedWord(typename PackedWord<T>::Scalar v)
    : w_{v}
{}

/**
 * Combine a pair of PackedWords into a single PackedWord. The values used are the bits
 * of the first item of p starting at bit a and the first a bits of the second item of p.
 *  [..|aaaaaa][bb|......] --> [aaaaaabb]
 * @param p Pair of PackedWords
 * @param a Index of first bit (inclusive) of p.first to use, index of last bit (exclusive) of p.second to use
 */
template <std::unsigned_integral T>
PackedWord<T>::PackedWord(Pair p, Index a)
    : w_{(p.first.word() >> a) | (p.second.word() << (CHAR_BIT*sizeof(T) - a))}
{}

/**
 * Construct PackedWord with bits [a,b) set to 1. All other bits set to 0.
 *
 * @param a First index, inclusive, of mask
 * @param b Last index, exclusive, of mask
 */
template <std::unsigned_integral T>
PackedWord<T>::PackedWord(Index a, Index b)
    :w_{((static_cast<T>(1) << (b-a)) - 1) << a}
{}

/**
 * Construct PackedWord from vector of indices where bit is 1.
 *
 * @param v Vector of indices of set bits
 */
template <std::unsigned_integral T>
PackedWord<T>::PackedWord(typename PackedWord<T>::Index_Vec v)
    : w_{static_cast<T>(0)}
{
    for (auto x : v) set_bit(x);
}

/**
 * Construct PackedWord from array of values.
 *
 * @param v Array of values (0, 1), one for each bit of T.
 */
template <std::unsigned_integral T>
PackedWord<T>::PackedWord(typename PackedWord<T>::Scalar_Vec v)
    : w_{static_cast<T>(0)}
{
    std::size_t idx = 0;
    for (auto x : v) set_bit(x, ++idx);
}


/**
 * Assign w to storage.
 *
 * @param w Value to copy into storage.
 */
template <std::unsigned_integral T>
void PackedWord<T>::set(typename PackedWord<T>::Word w)
{
    w_ = w;
}

/**
 * Standard accessor.
 *
 * @return contents of storage
 */
template <std::unsigned_integral T>
typename PackedWord<T>::Word PackedWord<T>::word() const
{
    return w_;
}

/**
 * Set the bit at a given index to the given (0,1) value.
 * The value is explicitly masked to be a single bit.
 *
 * @param v value to set (masked to a single bit)
 * @param b index to set
 */
template <std::unsigned_integral T>
void PackedWord<T>::set_bit(typename PackedWord<T>::Scalar v,
                            typename PackedWord<T>::Index b)
{
    clear_bit(b);
    w_ |= (v & static_cast<T>(1)) << b;
}

/**
 * Set the bit at the given index to 1.
 *
 * @param b index to set
 */
template <std::unsigned_integral T>
void PackedWord<T>::set_bit(typename PackedWord<T>::Index b)
{
    w_ |= (static_cast<T>(1) << b);
}

/**
 * Set the bit at the given index to 0.
 *
 * @param b index to set
 */
template <std::unsigned_integral T>
void PackedWord<T>::clear_bit(typename PackedWord<T>::Index b)
{
    w_ &= ~(static_cast<T>(1) << b);
}

/**
 * Return bit at given index.
 *
 * @param b index to return
 */
template <std::unsigned_integral T>
typename PackedWord<T>::Scalar PackedWord<T>::get_bit(typename PackedWord<T>::Index b) const
{
    return (w_ >> b) & static_cast<T>(1);
}

/**
 * Swap the bits at the two given indices.
 *
 * @param a first index
 * @param b second index
 */
template <std::unsigned_integral T>
void PackedWord<T>::swap_bits(typename PackedWord<T>::Index a, typename PackedWord<T>::Index b)
{
    auto t = get_bit(a);
    set_bit(get_bit(b),a);
    set_bit(t,b);
}

/**
 * Perform a left rotate on storage (w_ << a) | (w_ >> (bitsize(T) - a))
 *
 * Bit 0 becomes bit a.
 * @param a index of rotation
 */
template <std::unsigned_integral T>
void PackedWord<T>::rotl(typename PackedWord<T>::Index a)
{
    w_ = std::rotl(w_, a);
}

/**
 * Perform a right rotate on storage (w_ >> a) | (w_ << (bitsize(T) - a))
 *
 * Bit a becomes bit 0.
 * @param a index of rotation
 */
template <std::unsigned_integral T>
void PackedWord<T>::rotr(typename PackedWord<T>::Index a)
{
    w_ = std::rotr(w_, a);
}

/**
 * Count number of set bits
 *
 * @return number of set bits
 */
template <std::unsigned_integral T>
std::size_t PackedWord<T>::popcount() const
{
    return std::popcount(w_);
}

/**
 * Find first set bit.
 *
 * @return index of first set bit or CHAR_BIT*sizeof(T) if word is 0.
 */
template <std::unsigned_integral T>
std::size_t PackedWord<T>::ffs() const
{
    return std::countr_zero(w_);
}

/**
 * In-place bitwise negation.
 *
 */
template <std::unsigned_integral T>
void PackedWord<T>::negate_in()
{
    w_ = ~w_;
}

/**
 * Extract indices of all set bits.
 *
 * @return vector of indices containing locations of set bits
 */
template <std::unsigned_integral T>
typename PackedWord<T>::Index_Vec PackedWord<T>::scatter() const
{
    typename PackedWord<T>::Index_Vec v(popcount());
    for (auto i = 0; i < PackedWord<T>::BITSIZE; ++i) {
        if (get_bit(i)) v.emplace_back(i);
    }
    return v;
}

/**
 * Convert bitarray of storage to array of Ts.
 *
 * @return array of Ts containing values of each bit
 */
template <std::unsigned_integral T>
typename PackedWord<T>::Scalar_Vec PackedWord<T>::gather() const
{
    typename PackedWord<T>::Index_Vec v;
    for (auto i = 0; i < PackedWord<T>::BITSIZE; ++i) {
        v[i] = get_bit(i);
    }
    return v;
}

/**
 * Reverse contents of storage.
 *
 */
template <std::unsigned_integral T>
void PackedWord<T>::reverse()
{
    T w = static_cast<T>(0);
    for (auto i = 0; i < PackedWord<T>::BITSIZE; ++i) {
        w |= get_bit(i) << (Packed_word<T>::BITSIZE - i - 1);
    }
    w_ = w;
}

/**
 * In-place vector addition
 *
 * @param rhs PackedWord to add
 * @return this
 */
template <std::unsigned_integral T>
PackedWord<T>& PackedWord<T>::operator += (typename PackedWord<T>::Vector rhs)
{
    w_ += rhs.word();
    return *this;
}

/**
 * In-place vector subtraction
 *
 * @param rhs PackedWord to subtract
 * @return this
 */
template <std::unsigned_integral T>
PackedWord<T>& PackedWord<T>::operator -= (typename PackedWord<T>::Vector rhs)
{
    w_ -= rhs.word();
    return *this;
}

/**
 * In-place bitwise-or
 *
 * @param rhs PackedWord to or
 * @return this
 */
template <std::unsigned_integral T>
PackedWord<T>& PackedWord<T>::operator |= (typename PackedWord<T>::Vector rhs)
{
    w_ |= rhs.word();
    return *this;
}

/**
 * In-place bitwise-and
 *
 * @param rhs PackedWord to and
 * @return this
 */
template <std::unsigned_integral T>
PackedWord<T>& PackedWord<T>::operator &= (typename PackedWord<T>::Vector rhs)
{
    w_ &= rhs.word();
    return *this;
}

/**
 * In-place xor
 *
 * @param rhs PackedWord to xor
 * @return this
 */
template <std::unsigned_integral T>
PackedWord<T>& PackedWord<T>::operator ^= (typename PackedWord<T>::Vector rhs)
{
    w_ ^= rhs.word();
    return *this;
}

/**
 * In-place left bitwise shift
 *
 * @param rhs number of indices to shift left
 * @return this
 */
template <std::unsigned_integral T>
PackedWord<T>& PackedWord<T>::operator <<= (typename PackedWord<T>::Index rhs)
{
    w_ <<= rhs;
    return *this;
}

/**
 * In-place right bitwise shift
 *
 * @param rhs number of indices to shift right
 * @return this
 */
template <std::unsigned_integral T>
PackedWord<T>& PackedWord<T>::operator >>= (typename PackedWord<T>::Index rhs)
{
    w_ >>= rhs;
    return *this;
}

/**
 * In-place vector modulo
 *
 * @param rhs modulus
 * @return this
 */
template <std::unsigned_integral T>
PackedWord<T>& PackedWord<T>::operator %= (typename PackedWord<T>::Scalar rhs)
{
    w_ %= rhs;
    return *this;
}

template <std::unsigned_integral T>
PackedWord<T> PackedWord<T>::operator ~ () const
{
    return PackedWord<T>{~w_};
}

/**
 * Friend operator of +.
 *
 * @param lhs PackedWord
 * @param rhs PackedWord
 * @return new PackedWord: lhs + rhs
 */
template <std::unsigned_integral T>
PackedWord<T> operator + (PackedWord<T> lhs, PackedWord<T> rhs)
{
    PackedWord<T> x(lhs);
    x += rhs;
    return x;
}

/**
 * Friend operator of -.
 *
 * @param lhs PackedWord
 * @param rhs PackedWord
 * @return new PackedWord: lhs - rhs
 */
template <std::unsigned_integral T>
PackedWord<T> operator - (PackedWord<T> lhs, PackedWord<T> rhs)
{
    PackedWord<T> x(lhs);
    x -= rhs;
    return x;
}

/**
 * Friend operator of |.
 *
 * @param lhs PackedWord
 * @param rhs PackedWord
 * @return new PackedWord: lhs | rhs
 */
template <std::unsigned_integral T>
PackedWord<T> operator | (PackedWord<T> lhs, PackedWord<T> rhs)
{
    PackedWord<T> x(lhs);
    x |= rhs;
    return x;
}

/**
 * Friend operator of &.
 *
 * @param lhs PackedWord
 * @param rhs PackedWord
 * @return new PackedWord: lhs & rhs
 */
template <std::unsigned_integral T>
PackedWord<T> operator & (PackedWord<T> lhs, PackedWord<T> rhs)
{
    PackedWord<T> x(lhs);
    x &= rhs;
    return x;
}

/**
 * Friend operator of ^.
 *
 * @param lhs PackedWord
 * @param rhs PackedWord
 * @return new PackedWord: lhs ^ rhs
 */
template <std::unsigned_integral T>
PackedWord<T> operator ^ (PackedWord<T> lhs, PackedWord<T> rhs)
{
    PackedWord<T> x(lhs);
    x ^= rhs;
    return x;
}

/**
 * Friend operator of <<.
 *
 * @param lhs PackedWord
 * @param rhs index
 * @return new PackedWord: lhs << rhs
 */
template <std::unsigned_integral T>
PackedWord<T> operator << (PackedWord<T> lhs, typename PackedWord<T>::Index rhs)
{
    PackedWord<T> x(lhs);
    x <<= rhs;
    return x;
}

/**
 * Friend operator of >>.
 *
 * @param lhs PackedWord
 * @param rhs index
 * @return new PackedWord: lhs >> rhs
 */
template <std::unsigned_integral T>
PackedWord<T> operator >> (PackedWord<T> lhs, typename PackedWord<T>::Index rhs)
{
    PackedWord<T> x(lhs);
    x += rhs;
    return x;
}

/**
 * Friend operator of ==. Strict equality check of underlying storage of both operands.
 *
 * @param lhs PackedWord
 * @param rhs PackedWord
 * @return bool
 */
template <std::unsigned_integral T>
bool operator == (PackedWord<T> lhs, PackedWord<T> rhs)
{
    return lhs.word() == rhs.word();
}

/**
 * Fills storage with as many copies of v as possible, where we use the first s bits of v.
 * Fill occurs from least to most significant bits (starting at bit 0)
 * Unused bits (i.e., if s doesn't divide BITSIZE) are cleared.
 *
 * @param v scalar value to copy
 * @param s number of bits of v to use
 */
template <std::unsigned_integral T>
void PackedWord<T>::fill(typename PackedWord<T>::Scalar v, typename PackedWord<T>::Index s)
{
    w_ = static_cast<T>(0);
    T mask = (static_cast<T>(1) << s) - 1;
    T nv = static_cast<T>(v) & mask;
    
    for (auto i = 0; i <= PackedWord<T>::BITSIZE - s; i += s) {
        w_ |= (nv << i);
    }
}

/**
 * Inverse of fill operation. Fill a vector s-bit entries of our storage.
 * Only return full words.
 * In essence, ret[0] = bits 0 .. s-1, ret[1] = bits s .. 2s-1, ... ret[n] = bits n..n+s-1
 *
 * @param s number of bits per entry to use
 * @return std::vector of Ts of whole s-bit entries from the PackedWord
 */
template <std::unsigned_integral T>
typename PackedWord<T>::Scalar_Vec PackedWord<T>::deposit(typename PackedWord<T>::Index s)
{
    auto count = PackedWord<T>::BITSIZE / s;
    T mask = (static_cast<T>(1) << s) - 1;
    typename PackedWord<T>::Scalar_Vec x(count);
    for (auto i = 0; i < count; ++i) {
        x.emplace_back((w_ >> (s*i)) & mask);
    }
    return x;
}

/**
 * Inverse of constructor from Pair and Index.
 * Takes the BITSIZE bits of our storage and splits these into a pair of PackedWords,
 * such that bit 0 of our storage is copied to bit a of pair.first
 * if w_ is[12345678]  and a is 3, we construct
 * [---12345][678-----]
 *
 * @param a index of offset
 * @return Pair of PackedWords
 */
template <std::unsigned_integral T>
typename PackedWord<T>::Pair PackedWord<T>::split (typename PackedWord<T>::Index a) const
{
    return std::make_pair(PackedWord<T>(w_ << a),
                          PackedWord<T>(w_ >> (PackedWord<T>::BITSIZE - a)));
}

/**
 * Prints contents to given ostream, 1 bit at a time, from least to most significant.
 *
 * @param os std::ostream reference
 */
template <std::unsigned_integral T>
void PackedWord<T>::print (std::ostream &os) const
{
    for (auto i = 0; i < PackedWord<T>::BITSIZE; ++i) {
        os << get_bit(i);
    }
}

/**
 * Prints contents to given ostream, 1 bit at a time, from least to most significant.
 *
 * @param os std::ostream reference
 * @param rhs PackedWord
 * @return reference to os parameter
 */
template <std::unsigned_integral T>
std::ostream& operator << (std::ostream &os, PackedWord<T> rhs)
{
    rhs.print(os);
    return os;
}

#endif /* packed_word_inl */
