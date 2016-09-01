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

def problem3():
    """ What is the largest prime factor of the number 600851475143?"""
    def _prime_factorization(n):
        """Returns the list of prime factors of a number n"""
        factors = []
        f = 2
        # Use trial division to add factors
        while f**2 <= n:
            while (n % f) == 0:
                factors.append(f)
                n //= f
            f += 1

        if n > 1:
            factors.append(n)

        return factors

    return max(_prime_factorization(600851475143))

def problem5():
    """What is the smallest positive number that is evenly divisible
    by all the numbers from 1 to 20?"""
    # To solve this, we enumerate the factors of the number in queston.
    # First, each of the primes less than 20 must be a factor
    factors = [2, 3, 5, 7, 11, 13, 17, 19]

    #Next we add factors for each of the numers that can't be factored 
    # using just unique primes
    # 4 = 2*2
    factors.append(2)
    # 8 = 2*2*2
    factors.append(2)
    # 9 = 3*3
    factors.append(3)
    #16 = 2*2*2*2
    factors.append(2)

    return reduce(lambda x,y: x*y, factors)

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


