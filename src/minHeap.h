/*
 * minHeap.h
 *
 *  Created on: Aug 24, 2016
 *      Author: elias
 */

#ifndef SRC_MINHEAP_H_
#define SRC_MINHEAP_H_
#include <vector>
#include <functional>
#include <ctime>

template<class T, typename T2>
class MinHeap {
private:
    std::vector<T> heap;
    T2 comparator;

public:
    MinHeap(T2 comp);
    void push(T const& c);
    T top();
    void pop();
    void update();
    bool empty();
};

template<class T, typename T2>
MinHeap<T, T2>::MinHeap(T2 comp) {
    this->comparator = comp;
}

template<class T, typename T2>
void MinHeap<T, T2>::push(T const& elem) {
    // append copy of passed element
    heap.push_back(elem);
    std::push_heap(heap.begin(), heap.end(), comparator);
}

template<class T, typename T2>
T MinHeap<T, T2>::top() {
    return heap.front();
}

template<class T, typename T2>
void MinHeap<T, T2>::pop() {
    std::pop_heap(heap.begin(), heap.end(), comparator);
    heap.pop_back();

}

template<class T, typename T2>
void MinHeap<T, T2>::update() {
    std::make_heap(heap.begin(), heap.end(), comparator);
}

template<class T, typename T2>
bool MinHeap<T, T2>::empty() {
    return heap.empty();

}

#endif /* SRC_MINHEAP_H_ */
