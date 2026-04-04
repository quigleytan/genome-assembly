/*
 * KmerTypes.h
 * Summary:
 * - Shared type definitions for kmer and node encoding.
 * - NodeId is the encoded representation of a k-1 mer used as a graph node key.
 * - Centralizes the encoding type so DeBruijnGraph and EulerianPath stay in sync.
 * Important notes:
 * - __uint128_t is supported on GCC and Clang (x86-64) but not MSVC.
 * - Max k value is 64 with 2-bit encoding on __uint128_t.
 */

#ifndef KMER_TYPES_H
#define KMER_TYPES_H

using NodeId = __uint128_t;

#endif //KMER_TYPES_H