from __future__ import division

import unittest


## Sort an array using a heap O(n log(n))
## see Intro to Algorithms, Ch. 6.1 - 6.4, pp 151
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



if __name__ == '__main__':
	unittest.main()
