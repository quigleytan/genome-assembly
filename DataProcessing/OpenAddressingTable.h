/*
 * OpenAddressingTable.h
 * Summary:
 * - General linear probing hash table.
 * Important notes:
 * - Find and insert return a pointer to the item.
 */

#ifndef OPEN_ADDRESSING_TABLE_H
#define OPEN_ADDRESSING_TABLE_H

#include <string>
#include <vector>
#include <cstdint>

template <typename Key, typename Value>

class OpenAddressingTable {
protected:
    // Enumerated types declaration
    enum state {EMPTY, TAKEN};

    // Struct representing items in the table
    struct Item {
        Key key{};
        Value value{};
        state status = EMPTY;
    };

    // Variables
    size_t numItems; // Number of items in the table
    std::vector<Item> items;

    /**
     * @brief Computes a hash value for the given key.
     *
     * For __uint128_t inputs, the key is split into high and low 64-bit halves, each
     * mixed with bitwise operations and multipliers to reduce collisions.
     * For other key types, `std::hash` is used. The result is modulo the table size.
     *
     * @tparam Key Type of the key.
     * @param key  Key to hash.
     * @return Hash value in the range [0, table size).
     */
    virtual size_t hashKey(const Key& key) const {
        if constexpr (std::is_same_v<Key, __uint128_t>) {
            uint64_t lo = static_cast<uint64_t>(key);
            uint64_t hi = static_cast<uint64_t>(key >> 64);
            // Mix each half independently first
            hi ^= hi >> 33;
            hi *= 0xff51afd7ed558ccdULL;
            hi ^= hi >> 33;
            lo ^= lo >> 33;
            lo *= 0xc4ceb9fe1a85ec53ULL;
            lo ^= lo >> 33;
            // Combine
            return (lo ^ (hi * 0x9e3779b97f4a7c15ULL)) % items.size();
        } else {
            return std::hash<Key>{}(key) % items.size();
        }
    }

    /**
     * @brief Raw insert method for rehashing.
     *
     * Bypasses onInitial and onDuplicate hooks entirely.
     * Assumes the key does not already exist in the table.
     */
    void rawInsert(const Key& key, const Value& value) {
        size_t index = hashKey(key);
        size_t i = 0;
        while (items[index].status == TAKEN) {
            i++;
            index = (index + i) % items.size();
        }
        items[index].key    = key;
        items[index].value  = value;
        items[index].status = TAKEN;
        ++numItems;
    }

    /**
     * @brief Rebuilds and resizes the hash table to accommodate increased requirements.
     *
     * This method is called when the load factor exceeds 0.5, resizing the table
     * and reinserting all prior items in the table.
     */
    virtual void rehash() {
        std::vector<Item> preRehash = items;
        items.clear();
        numItems = 0;
        items.resize(nextPrime(preRehash.size() * 2));

        for (const auto& item : preRehash) {
            if (item.status == TAKEN) {
                rawInsert(item.key, item.value); // ← no hooks, no corruption
            }
        }
    }

    // HOOK FUNCTIONS FOR KmerTable

    /**
     * @brief Exists to be overridden for kmer increment.
     * @param value
     */
    virtual void onDuplicate(Value& value) {
        // Default behavior: do nothing (keep existing value)
    }

    /**
     * @brief Exists to be overridden for initial kmer insertion.
     * @param value
     */
    virtual void onInitial(Value& value) {
        // Default behavior: do nothing (keep existing value)
    }

    /**
     * @brief Finds the next prime number greater than or equal to n
     * @param n Starting integer to find the next prime
     * @return The next prime number.
     */
    static size_t nextPrime(size_t n) {
        auto isPrime = [](int x) {
            if (x < 2) return false;
            for (int i = 2; i*i <= x; ++i)
                if (x % i == 0) return false;
            return true;
        };
        while (!isPrime(n)) ++n;
        return n;
    }

public:

    class iterator {

        public:
            /**
             * @brief Constructor for an iterator that allows for simple and complete scanning of the table.
             * @param inputTable Input table to be iterated through.
             * @param inputIndex The index to start iteration at.
             */
            iterator(const OpenAddressingTable* inputTable, size_t inputIndex)
                : table(inputTable), index(inputIndex) {
                advanceToNextValid();
            }

            /**
             * @brief Advances the iterator to the next occupied item.
             * @return Reference to the updated iterator.
             */
            iterator& operator++() {
                ++index;
                advanceToNextValid();
                return *this;
            }

            /**
             * @brief Compares two iterators for inequality.
             * @param other Iterator to compare against.
             * @return True if the iterators refer to different positions.
             */
            bool operator!=(const iterator& other) const {
                return index != other.index;
            }

            /**
             * @brief Returns the item at the current iterator position.
             * @return Reference to the current item.
             */
            const Item& operator*() const {
                return table->items[index];
            }

        private:
            const OpenAddressingTable* table;
            size_t index;

            /**
             * @brief Advances through the vector until the next item is found.
             */
            void advanceToNextValid() {
                while (index < table->items.size() &&
                       table->items[index].status != TAKEN) {
                    ++index;
                }
            }

        };

        /**
         * @brief Returns an iterator to the first occupied table entry.
         * @return Iterator to the beginning of the table.
         */
        iterator begin() const {
            return iterator(this, 0);
        }

        /**
         * @brief Returns an iterator one past the end of the table.
         * @return Iterator representing the end position.
         */
        iterator end() const {
            return iterator(this, items.size());
        }

    // OPEN ADDRESSING TABLE METHODS
    virtual ~OpenAddressingTable() = default;

    /**
     * @brief Main constructor for OpenAddressingTable.
     *
     * Initializes the hash table with a specified initial size, which is adjusted to the next prime number to
     * help reduce collisions.
     * The number of items is initialized to zero.
     *
     * @param initialSize initial size for construction of the hash table.
     */
    OpenAddressingTable(size_t initialSize = 101) {
        items.resize(nextPrime(initialSize));
        numItems = 0;
    }

    // Getter method

    /**
     * @brief Returns the number of items in the table.
     * @return numItems The number of items in the table.
     */
    [[nodiscard]] size_t getNumItems() const {
        return numItems;
    }

    /**
     * @brief Inserts an item into the table with the specified key.
     *
     * Inserts an item into the table using a linear probing strategy to resolve collisions.
     * If the key already exists, the onDuplicate method is called to handle the duplicate key scenario.
     *
     * @param insertKey
     * @return Value& Reference to the value associated with the inserted key.
     */
    virtual std::pair<Value&, bool> insert(const Key& insertKey) {

        // Makes sure references are correct after rehashing
        if (numItems + 1 > items.size()/2) {
            rehash();
        }

        size_t index = hashKey(insertKey);
        size_t i = 0;
        const size_t maxProbes = items.size();

        while (items[index].status == TAKEN && items[index].key != insertKey) {
            i++;
            index = (index + i) % items.size(); // Linear probing strategy.

            if (i >= maxProbes) // Announces error in table filling.
                throw std::runtime_error("Hash table is full — rehash failed or table undersized");
        }

        if (items[index].status == TAKEN) {
            // Return reference to an existing value
            onDuplicate(items[index].value);
            return {items[index].value, false};
        }
        items[index].key = insertKey;
        items[index].value = Value{};
        items[index].status = TAKEN;
        ++numItems;

        onInitial(items[index].value);

        return {items[index].value, true};
    }

    /**
     * @brief Searches the hash table for an entry with the specified key.
     *
     * CONST METHOD: Does not allow mutation of pointer information.
     * Performs a lookup using the same probing strategy as insertion (linear probing)
     * starting from the hashed index of the key.
     *
     * @param searchKey Key to search for.
     * @return Pointer to the value associated with key if found; otherwise nullptr.
     */
    virtual const Value* find(const Key& searchKey) const {
        size_t index = hashKey(searchKey);
        size_t i = 0;
        while (items[index].status != EMPTY) {
            if (items[index].status == TAKEN && items[index].key == searchKey)
                return &items[index].value;
            i++;
            index = (index + i) % items.size();
        }
        return nullptr;
    }

    /**
     * @brief Searches the hash table for an entry with the specified key.
     *
     * NON-CONST OVERRIDE METHOD: Context dependent, allows for value mutation.
     * Performs a lookup using the same probing strategy as insertion (linear probing)
     * starting from the hashed index of the key.
     *
     * @param searchKey Key to search for.
     * @return Pointer to the value associated with key if found; otherwise nullptr.
     */
    virtual Value* find(const Key& searchKey) {
        size_t index = hashKey(searchKey);
        size_t i = 0;
        while (items[index].status != EMPTY) {
            if (items[index].status == TAKEN && items[index].key == searchKey)
                return &items[index].value;
            i++;
            index = (index + i) % items.size();
        }
        return nullptr;
    }

};
#endif //OPEN_ADDRESSING_TABLE_H