/* 
 * File:   sorted_vector.h
 * Author: eannila
 *
 * Created on 22. syyskuuta 2016, 9:41
 */

#ifndef SORTED_VECTOR_H
#define	SORTED_VECTOR_H

#include <vector>
#include <functional>

template<class T>
class SortedVector {
private:
    std::vector<T> vec;

public:
    std::pair<typename std::vector<T>::iterator, bool> insert(T const& c);
    typename std::vector<T>::iterator begin();
    typename std::vector<T>::iterator end();
    typename std::vector<T>::iterator find(T const& elem);
    bool empty();
    int size();
};

template<class T>
std::pair<typename std::vector<T>::iterator, bool> SortedVector<T>::insert(T const& elem) {

    // append copy of passed element
    if (vec.empty()) {
        vec.push_back(elem);
        return std::make_pair(vec.begin(), true);
    } else {
        auto it = std::lower_bound(vec.begin(), vec.end(), elem);
        if (*it != elem) {
            it = vec.insert(it, elem);
            return std::make_pair(it, true);
        }else{
            return std::make_pair(it, false);
        }
    }
}

template<class T>
typename std::vector<T>::iterator SortedVector<T>::begin() {
    return vec.begin();
}

template<class T>
typename std::vector<T>::iterator SortedVector<T>::end() {
    return vec.end();
}

template<class T>
typename std::vector<T>::iterator SortedVector<T>::find(T const& elem) {
    auto it = std::lower_bound(vec.begin(), vec.end(), elem);
    if (*it == elem) {
        return it;
    }
    return vec.end();
}

template<class T>
bool SortedVector<T>::empty() {
    return vec.empty();
}

template<class T>
int SortedVector<T>::size() {
    return vec.size();
}


#endif	/* SORTED_VECTOR_H */

