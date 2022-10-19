//
//  PackedWordRef.hpp
//  sp4m
//
//  Created by Matthew Lambert on 10/14/22.
//

#ifndef PackedWordRef_h
#define PackedWordRef_h

#include <array>
#include <cstdlib>
#include <utility>
#include <vector>
#include <bit>
#include <climits>
#include <cstdint>
#include <concepts>
#include <iostream>

template <std::unsigned_inegral T>
class PackedWordRef
{
public:
    using Word = PackedWord<T>;
    using Index = std::size_t;
    using Ref = PackedWordRef;
    using Scalar = typename Word::Scalar;
    
private:
    Word &v;
    Index idx;
    
public:
    PackedWordRef() = delete;
    PackedWordRef (Word &w);
    PackedWordRef& operator = (Scalar x);
    PackedWordRef (PackedWordRef &&other);
    PackedWordRef (const PackedWordRef &other);
    ~PackedWordRef() = default;
    
public:

};


#endif /* PackedWordRef_h */
