#ifndef PRINT_UTIL_HPP
#define PRINT_UTIL_HPP

#include <array>
#include <cstddef>
#include <deque>
#include <iostream>
#include <iterator>
#include <map>
#include <queue>
#include <set>
#include <stack>
#include <unordered_map>
#include <unordered_set>
#include <utility>
#include <vector>

// Vector
template <typename T>
std::ostream& operator<<(std::ostream& os, const std::vector<T>& v) {
    os << "[";
    for (std::size_t i = 0; i < v.size(); ++i) {
        os << v[i];
        if (i < v.size() - 1) os << ", ";
    }
    os << "]";
    return os;
}

// Set
template <typename T>
std::ostream& operator<<(std::ostream& os, const std::set<T>& s) {
    os << "{";
    for (auto it = s.begin(); it != s.end(); ++it) {
        os << *it;
        if (std::next(it) != s.end()) os << ", ";
    }
    os << "}";
    return os;
}

// Unordered Set
template <typename T>
std::ostream& operator<<(std::ostream& os, const std::unordered_set<T>& s) {
    os << "{";
    for (auto it = s.begin(); it != s.end(); ++it) {
        os << *it;
        if (std::next(it) != s.end()) os << ", ";
    }
    os << "}";
    return os;
}

// Pair
template <typename T1, typename T2>
std::ostream& operator<<(std::ostream& os, const std::pair<T1, T2>& p) {
    os << "(" << p.first << ", " << p.second << ")";
    return os;
}

// Map
template <typename K, typename V>
std::ostream& operator<<(std::ostream& os, const std::map<K, V>& m) {
    os << "{";
    for (auto it = m.begin(); it != m.end(); ++it) {
        os << it->first << ": " << it->second;
        if (std::next(it) != m.end()) os << ", ";
    }
    os << "}";
    return os;
}

// Unordered Map
template <typename K, typename V>
std::ostream& operator<<(std::ostream& os, const std::unordered_map<K, V>& m) {
    os << "{";
    for (auto it = m.begin(); it != m.end(); ++it) {
        os << it->first << ": " << it->second;
        if (std::next(it) != m.end()) os << ", ";
    }
    os << "}";
    return os;
}

// Array
template <typename T, std::size_t N>
std::ostream& operator<<(std::ostream& os, const std::array<T, N>& arr) {
    os << "[";
    for (std::size_t i = 0; i < N; ++i) {
        os << arr[i];
        if (i < N - 1) os << ", ";
    }
    os << "]";
    return os;
}

// Deque
template <typename T>
std::ostream& operator<<(std::ostream& os, const std::deque<T>& d) {
    os << "[";
    for (std::size_t i = 0; i < d.size(); ++i) {
        os << d[i];
        if (i < d.size() - 1) os << ", ";
    }
    os << "]";
    return os;
}

#endif
