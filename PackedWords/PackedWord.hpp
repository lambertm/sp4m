//
//  PackedWord.hpp
//  sp4m
//
//  Created by Matthew Lambert on 10/13/22.
//

#ifndef packed_word_h
#define packed_word_h

#include <array>
#include <cstdlib>
#include <utility>
#include <vector>
#include <bit>
#include <climits>
#include <cstdint>
#include <concepts>
#include <iostream>

// to do:
// iterators + references
// PackedVector is wrapper over this with additional template param size

// forward declarations of friends
template <std::unsigned_integral T> class PackedWord;
template <std::unsigned_integral T> PackedWord<T> operator + (PackedWord<T> lhs, PackedWord<T> rhs);
template <std::unsigned_integral T> PackedWord<T> operator - (PackedWord<T> lhs, PackedWord<T> rhs);
template <std::unsigned_integral T> PackedWord<T> operator | (PackedWord<T> lhs, PackedWord<T> rhs);
template <std::unsigned_integral T> PackedWord<T> operator & (PackedWord<T> lhs, PackedWord<T> rhs);
template <std::unsigned_integral T> PackedWord<T> operator ^ (PackedWord<T> lhs, PackedWord<T> rhs);
template <std::unsigned_integral T> PackedWord<T> operator << (PackedWord<T> lhs, typename PackedWord<T>::Index rhs);
template <std::unsigned_integral T> PackedWord<T> operator >> (PackedWord<T> lhs, typename PackedWord<T>::Index rhs);
template <std::unsigned_integral T> bool operator == (PackedWord<T> lhs, PackedWord<T> rhs);
template <std::unsigned_integral T> std::ostream& operator << (std::ostream &os, PackedWord<T> rhs);

template <std::unsigned_integral T>
class PackedWord {
public: //
    static const std::size_t BITSIZE = CHAR_BIT * sizeof(T);
    
public: // typedefs
    using Scalar = T; // expresses single value
    using Word = T; // expresses entire storage
    using Vector = PackedWord<Scalar>;
    using Index = std::size_t;
    using Pair = std::pair<Vector,Vector>;
    using Scalar_Vec = std::array<Scalar, CHAR_BIT * sizeof(T)>;
    using Index_Vec = std::vector<std::size_t>;
    
private: // data members
    T w_;

public: // defaults
    PackedWord() = default;
    PackedWord(const PackedWord &other) = default;
    PackedWord& operator = (const PackedWord &other) = default;
    PackedWord& operator = (PackedWord &&other) = default;
    ~PackedWord() = default;
    
public: // helpful constructors
    PackedWord(Scalar v);
    PackedWord(Pair p, Index a);
    PackedWord(Index a, Index b);
    PackedWord(Index_Vec v);
    PackedWord(Scalar_Vec v);
    
public: // standard ops
    void set(Word w);
    Word word() const;
    void set_bit(Scalar v, Index b);
    void set_bit(Index b);
    void clear_bit(Index b);
    Scalar get_bit(Index b) const;
    void swap_bits(Index a, Index b);
    void rotl(Index a);
    void rotr(Index a);
    std::size_t popcount() const;
    Index ffs() const;
    void negate_in();
    
public: // other
    Index_Vec scatter() const;
    Scalar_Vec gather() const;
    void reverse();
    
public: // operators
    PackedWord& operator += (Vector rhs);
    PackedWord& operator -= (Vector rhs);
    PackedWord& operator |= (Vector rhs);
    PackedWord& operator &= (Vector rhs);
    PackedWord& operator ^= (Vector rhs);
    
    PackedWord& operator <<= (Index rhs);
    PackedWord& operator >>= (Index rhs);
    PackedWord& operator %= (Scalar rhs);
    
    PackedWord operator ~() const;
    
    friend PackedWord operator + <>(PackedWord lhs, PackedWord rhs);
    friend PackedWord operator - <>(PackedWord lhs, PackedWord rhs);
    friend PackedWord operator | <>(PackedWord lhs, PackedWord rhs);
    friend PackedWord operator & <>(PackedWord lhs, PackedWord rhs);
    friend PackedWord operator ^ <>(PackedWord lhs, PackedWord rhs);
    
    friend PackedWord operator << <>(PackedWord lhs, Index rhs);
    friend PackedWord operator >> <>(PackedWord lhs, Index rhs);
        
    friend bool operator == <>(PackedWord lhs, PackedWord rhs);
    
    
    

public: // builders
    void fill(Scalar v, Index s);
    Scalar_Vec deposit(Index s);
    Pair split (Index a) const;
    
public: // io
    void print (std::ostream &os) const;
    friend std::ostream& operator << <>(std::ostream &os, PackedWord rhs);
};

#include "PackedWord.inl"

#endif /* packed_word_h */
