#! /usr/bin/python3
import argparse
import inspect

"""
Solving the project Euler problems for fun in python 3
"""


def problem1():
    """ Find the sum of all the multiples of 3 or 5 below 1000. """
    return sum(i for i in range(1000) if i % 3 == 0 or i % 5 == 0)




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


