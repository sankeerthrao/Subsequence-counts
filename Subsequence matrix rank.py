#This is the code for the rank deficiency of the subsequence incidence matrices.
import math
import numpy as np

#We use this to compute the rank.

from numpy.linalg import matrix_rank as mr

#We use this function for subsequence counts - it only works when the lengths are off by 1 as we need.  

def subseq_count(x,y):
    count = 0
    for i in range(1, len(y)+1):
        if x == y[:i-1]+y[i:]:
            count = count+1
    return count

#We construct two arrays A,B.

#This constructs an array with all strings of length k and weight m
def constweight(prefix,A,k,m):
    if m == 0:
        A.append(prefix + "0"*k)
        return
    if k == m:
        A.append(prefix + "1"*k)
        return
    constweight(prefix + "0",A,k-1, m)
    constweight(prefix + "1",A,k-1, m-1)

#We construct the incidence matrix here.

def incidence_matrix(k,m):
    A = []
    B = []
    constweight("",A,k-1,m)
    constweight("",A,k-1,m-1)
    constweight("",B,k,m)
    M = np.zeros((binom(k,m),binom(k,m)))
    for i in range(0,len(A)):
        for j in range(0,len(B)):
            M[i][j] = subseq_count(A[i],B[j])
    return M

def binom(n,k):
    return int(math.factorial(n)/(math.factorial(k)*math.factorial(n-k)))

#We compute the rank deficiency here - inputs are the length of the strings and the hamming weight. 
def rank_deficiency(k,m):
    return binom(k,m)-mr(incidence_matrix(k,m))

for k in range (3,10):
    for m in range (1,k):
        print("At length ", k ," and weight", m ," the rank deficiency is ", rank_deficiency(k,m))
