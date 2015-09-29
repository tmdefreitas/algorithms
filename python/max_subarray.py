import unittest


## Maximum common subarray problem O(n log(n))
## see Intro to Algorithms, Ch. 4.1, pp 68
class MaxSubArray(object):

	def __init__(self, A):
		self.A = A
		self.low, self.high, self.max = self.max_subbarray(0, len(self.A)-1)


	def max_crossing_subarray(self, low, mid, high):
		left_sum = None
		total = 0
		max_left = None
		for i in range(mid,low-1,-1):
			
			total += self.A[i]
			if total > left_sum:
				left_sum = total
				max_left = i

		right_sum = None
		total = 0
		max_right = None
		for j in range(mid+1,high+1):
			
			total += self.A[j]
			if total > right_sum:
				right_sum = total
				max_right = j

		##Could be that there is no max_crossing subarray		
		result = (None, None, None)
		if left_sum is not None and right_sum is not None:
			result = (max_left, max_right, left_sum+right_sum)

		return result

	def max_subbarray(self,low,high):
		if low == high:
			return (low, high, self.A[low]) ##base case
		else:
			mid = (low + high) / 2 ##Midpoint in array (integer division)
			left_low, left_high, left_sum = self.max_subbarray(low,mid) ##left subproblem
			right_low, right_high, right_sum = self.max_subbarray(mid+1, high) ##right subproblem
			cross_low, cross_high, cross_sum = self.max_crossing_subarray(low,mid,high) ##crossing

			
			result = None
			if left_sum >= right_sum and left_sum >= cross_sum:
				result = left_low, left_high, left_sum
			elif right_sum >= left_sum and right_sum >= cross_sum:
				result = right_low, right_high, right_sum
			else: 
				result =  cross_low, cross_high, cross_sum
			return result


class TestMaxSubArray(unittest.TestCase):

	def test_one_element(self):
		self.assertEqual(MaxSubArray([3]).max, 3)

	def test_one_negative(self):
		self.assertEqual(MaxSubArray([-2]).max, -2)

	def test_discontinuous(self):
		self.assertEqual(MaxSubArray([-2, 3, 7, -2, 8, -3]).max, 16)

	def test_all_negative(self):
		self.assertEqual(MaxSubArray([-2, -8, -19, -12, -5, -90]).max, -2)

	def test_all_positive(self):
		self.assertEqual(MaxSubArray([4, 2, 10, 19, 6]).max, 41)




if __name__ == '__main__':
	unittest.main()
