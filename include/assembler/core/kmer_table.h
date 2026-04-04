/*
 * kmer_table.h
 * Summary:
 * - Fixed-size open addressing hash table used to count unique k-mers in
 *   DNA sequences.
 * - Uses linear probing for collision resolution and flexible size initialization
 *   to minimize collisions and ensure reasonable allocation.
 * Important notes:
 * - Uses 2-bit encoding for kmers, functions as the key and the kmer identity.
 * - Stores number of instances as size_t.
 * - Based functions off of OpenAddressing from CS 2240 Project 5.
 */

#ifndef KMER_TABLE_H
#define KMER_TABLE_H

#include <cstdint>
#include <vector>
#include <algorithm>
#include <stdexcept>

#include "kmer_types.h"
#include "encoding/open_addressing_table.h"

class KmerTable : public HashTable<NodeId, size_t> {

public:

    /**
     * @brief Constructor for QuadraticProbingTable.
     *
     * Uses the length of the input sequence and k to decide between the physical and
     * combinational limits of unique kmers using: min(4^k, dnaLength - k + 1).
     *
     * Ensures the max load is < 0.5 by initializing table size to the
     * next prime > 2 * maxUnique.
     *
     * @param dnaLength Length of the DNA sequence (bases).
     * @param k         Desired k-mer size (bases).
     */
    KmerTable(size_t dnaLength, size_t k) {
        if (k > dnaLength) {
            items.resize(nextPrime(1));
            numItems = 0;
            return;
        }
        size_t maxPossible = (k >= 32) ? SIZE_MAX : (size_t(1) << (2 * k));
        size_t maxUnique = std::min(maxPossible, dnaLength - k + 1);
        size_t tableSize = (maxUnique > SIZE_MAX / 2) ? nextPrime(maxUnique) : nextPrime(2 * maxUnique);
        items.resize(tableSize);
        numItems = 0;
        this->k = k;
    }

    /**
     * @breif Getter method for k.
     * @return k The kmer size for this table.
     */
    [[nodiscard]] size_t getK() const {
        return k;
    }

protected:

    size_t k;

    /**
     * @brief Rehash overridden method.
     *
     * Overrides the rehash function to disable rehashing as the table is initialized
     * to handle the expected number of unique kmers.
     * Throws a runtime_error if the table was sized incorrectly and rehash is invoked.
     */
    void rehash() override {
        throw std::runtime_error(
            "kmer_table exceeded capacity — increase totalBases estimate at construction"
        );
    }

    /**
     * @brief onDuplicate overridden method.
     *
     * Overrides the onDuplicate function to increment the k-mer count when a duplicate key
     * is inserted.
     *
     * @param value Reference to the value associated with a duplicate key.
     */
    void onDuplicate(size_t &value) override {
        value++; // Increment count for duplicate k-mer
    }

    /**
     * @brief onInitial overridden method.
     *
     * Overrides the onInitial function to start the k-mer count at 1 when a new key
     * is inserted.
     *
     * @param value Reference to the value associated with a duplicate key.
     */
    void onInitial(size_t &value) override {
        value = 1;
    }

};
#endif //KMER_TABLE_H