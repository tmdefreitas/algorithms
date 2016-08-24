#! /usr/bin/python3
import argparse
import inspect

"""
Solving the project Euler problems for fun in python 3
"""


def problem1():
    """ Find the sum of all the multiples of 3 or 5 below 1000. """
    return sum(i for i in range(1000) if i % 3 == 0 or i % 5 == 0)

def problem2():
    """ Find the sum of the even fibonacci numbers that do not exceed
    4e6 """

    def _fib_iter(n=4000000):
        """  Generator for fibonacci numbers less than n  """
        fib1 = 1
        fib2 = 2
        # Yield the first two fibonacci numbers
        yield fib1
        yield fib2
        fib_next = fib1 + fib2
        while fib_next < n:
            # iteratively gen
            yield fib_next
            fib1 = fib2
            fib2 = fib_next
            fib_next = fib1 + fib2

    return sum(i for i in _fib_iter() if i % 2 == 0)




if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('problem', type=int)
    args = parser.parse_args()
    p = "problem" + str(args.problem)
    if p in globals():
        # Print the result of that problem
        print(globals()[p]())
    else:
        raise NotImplementedError(p + " has not been implemented!")


