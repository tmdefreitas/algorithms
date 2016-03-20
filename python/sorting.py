from __future__ import division

from random import randint
import unittest

#############################################################
## 
## Heapsort
##
## Sort an array using a heap O(n log(n))
## see Intro to Algorithms, Ch. 6.1 - 6.4, pp 151
#############################################################
def heap_sort(A):
    """
    Sorts input list A in place using heap operations
    """
    _max_heap(A)
    heap_size = len(A)
    for idx in range(len(A)-1, 0, -1):
        A[0], A[idx] = A[idx], A[0]
        heap_size -= 1
        _heapify(A, 0, heap_size)
    return A

def _max_heap(A):
    heap_size = len(A)
    for i in range(len(A), -1, -1):
        _heapify(A, i, heap_size)
    return A

def _heapify(A, i, heap_size):
    """
    Helper for heap sort which floats down the value at i if 
    it is smaller than the heaps at 2i and 2i + 1. Assumes the 
    heaps rooted at 2i and 2i + 1 are valid max heaps
    """
    size=heap_size
    left = 2*(i+1) - 1
    right = 2*(i+1)

    ## Pick which of the three is largest, which must be the new root
    if left < size and A[left] > A[i]:
        largest = left
    else:
        largest = i
    if right < size and A[right] > A[largest]:
        largest = right

    ##If i is not largest, replace, and call max-heapify on changed child
    if largest != i:
        A[i], A[largest] = A[largest], A[i] ##python is cool!
        _heapify(A, largest, heap_size)


class TestHeapsort(unittest.TestCase):
    def test_heapsort(self):
        self.assertEqual(heap_sort([5,2, 12, 80, 19]), [2, 5, 12, 19, 80])
        self.assertEqual(heap_sort([]), [])
        self.assertEqual(heap_sort([10,9,8,7,6,5,4,3,2,1]), [1,2,3,4,5,6,7,8,9,10])

#############################################################
## 
## Quicksort, randomized quicksort
##
## Sort an array using a heap O(n log(n))
## see Intro to Algorithms, Ch. 6.1 - 6.4, pp 151
#############################################################

#Wraps the real quicksort function
def quicksort(A):
    return _quicksort(A, 0, len(A)-1)

## Randomized quicksort, pick a random pivot element
def r_quicksort(A):
    return _r_quicksort(A, 0, len(A)-1)

def _quicksort(A, p, r):
    if p < r:
        q = _partition(A, p, r)
        _quicksort(A, p, q-1)
        _quicksort(A, q+1, r)
    return A

def _partition(A, p, r):
    pivot = A[r]
    i = p - 1
    for j in range(p, r):
        if A[j] <= pivot:
            i += 1
            A[i], A[j] = A[j], A[i]
    A[i+1], A[r] = A[r], A[i+1]
    return i + 1

def _r_quicksort(A, p, r):
    if p < r:
        q = _r_partition(A, p, r)
        _r_quicksort(A, p, q-1)
        _r_quicksort(A, q + 1, r)
    return A

def _r_partition(A, p, r):
    i = randint(p,r)
    A[r], A[i] = A[i], A[r]
    return _partition(A, p, r)

class TestQuicksort(unittest.TestCase):
    def test_quicksort(self):
        self.assertEqual(quicksort([5,2, 12, 80, 19]), [2, 5, 12, 19, 80])
        self.assertEqual(quicksort([]), [])
        self.assertEqual(quicksort([10,9,8,7,6,5,4,3,2,1]), 
                                   [1,2,3,4,5,6,7,8,9,10])

    def test_rquicksort(self):
        self.assertEqual(r_quicksort([5,2, 12, 80, 19]), [2, 5, 12, 19, 80])
        self.assertEqual(r_quicksort([]), [])
        self.assertEqual(r_quicksort([10,9,8,7,6,5,4,3,2,1]), 
                                     [1,2,3,4,5,6,7,8,9,10])


if __name__ == '__main__':
    unittest.main()
