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
    double heaptime = 0;
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
    std::clock_t begin = std::clock();
    heap.push_back(elem);
    std::push_heap(heap.begin(), heap.end(), comparator);
    std::clock_t end = std::clock();
    this->heaptime += double(end - begin) / CLOCKS_PER_SEC;
}

template<class T, typename T2>
T MinHeap<T, T2>::top() {
    std::clock_t begin = std::clock();
    return heap.front();
    std::clock_t end = std::clock();
    this->heaptime += double(end - begin) / CLOCKS_PER_SEC;
}

template<class T, typename T2>
void MinHeap<T, T2>::pop() {
    std::clock_t begin = std::clock();

    std::pop_heap(heap.begin(), heap.end(), comparator);
    heap.pop_back();
    std::clock_t end = std::clock();
    this->heaptime += double(end - begin) / CLOCKS_PER_SEC;
}

template<class T, typename T2>
void MinHeap<T, T2>::update() {
    std::clock_t begin = std::clock();
    std::make_heap(heap.begin(), heap.end(), comparator);
    std::clock_t end = std::clock();
    this->heaptime += double(end - begin) / CLOCKS_PER_SEC;
}

template<class T, typename T2>
bool MinHeap<T, T2>::empty() {
    std::clock_t begin = std::clock();
    return heap.empty();
    std::clock_t end = std::clock();
    this->heaptime += double(end - begin) / CLOCKS_PER_SEC;
}

#endif /* SRC_MINHEAP_H_ */
