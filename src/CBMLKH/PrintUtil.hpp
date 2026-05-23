#ifndef PRINT_UTIL_HPP
#define PRINT_UTIL_HPP

#include <array>
#include <deque>
#include <iostream>
#include <map>
#include <queue>
#include <set>
#include <stack>
#include <unordered_map>
#include <unordered_set>
#include <vector>

using namespace std;

// Vector
template <typename T>
ostream& operator<<(ostream& os, const vector<T>& v) {
    os << "[";
    for (size_t i = 0; i < v.size(); ++i) {
        os << v[i];
        if (i < v.size() - 1) os << ", ";
    }
    os << "]";
    return os;
}

// Set
template <typename T>
ostream& operator<<(ostream& os, const set<T>& s) {
    os << "{";
    for (auto it = s.begin(); it != s.end(); ++it) {
        os << *it;
        if (next(it) != s.end()) os << ", ";
    }
    os << "}";
    return os;
}

// Unordered Set
template <typename T>
ostream& operator<<(ostream& os, const unordered_set<T>& s) {
    os << "{";
    for (auto it = s.begin(); it != s.end(); ++it) {
        os << *it;
        if (next(it) != s.end()) os << ", ";
    }
    os << "}";
    return os;
}

// Pair
template <typename T1, typename T2>
ostream& operator<<(ostream& os, const pair<T1, T2>& p) {
    os << "(" << p.first << ", " << p.second << ")";
    return os;
}

// Map
template <typename K, typename V>
ostream& operator<<(ostream& os, const map<K, V>& m) {
    os << "{";
    for (auto it = m.begin(); it != m.end(); ++it) {
        os << it->first << ": " << it->second;
        if (next(it) != m.end()) os << ", ";
    }
    os << "}";
    return os;
}

// Unordered Map
template <typename K, typename V>
ostream& operator<<(ostream& os, const unordered_map<K, V>& m) {
    os << "{";
    for (auto it = m.begin(); it != m.end(); ++it) {
        os << it->first << ": " << it->second;
        if (next(it) != m.end()) os << ", ";
    }
    os << "}";
    return os;
}

// Array
template <typename T, size_t N>
ostream& operator<<(ostream& os, const array<T, N>& arr) {
    os << "[";
    for (size_t i = 0; i < N; ++i) {
        os << arr[i];
        if (i < N - 1) os << ", ";
    }
    os << "]";
    return os;
}

// Deque
template <typename T>
ostream& operator<<(ostream& os, const deque<T>& d) {
    os << "[";
    for (size_t i = 0; i < d.size(); ++i) {
        os << d[i];
        if (i < d.size() - 1) os << ", ";
    }
    os << "]";
    return os;
}

#endif
