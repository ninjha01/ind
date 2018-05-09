import mmh3
import statistics
import random
from progress.bar import Bar
import edit_distance as ed

'''
Input of A: a set A.
Input of B: a set B.
Public Input: a (sufficiently long) common random string

k and l are accuracy bounds

Output: |Diff(A, B)| is calculated as the k-median of l-mean of random sampling of (d_A - d_B)^2

1. For j in (1..k)
   a. For i in (1..l)
      i. A and B use the (same) random string to randomly pick a function h from the family of hash functions.
      ii. A computes d_A = Sum_{s \in A} h(s) and B computes d_b = Sum_{s \in B} h(s), both independently.
      iii. A and B securely compute D_i = (d_A - d_b)^2
   b. A and B securely compute D^_j = Sum(i -> l) D_i / l
2. A and B securely compute the median Z of D^_1 ... D^_k
3. Output Z
'''
def main():
    print("Reading ref genome from: ./genome_ref.txt")
    ref_gen = open("./genome_ref.txt", "r").read()
    print("Reading genome A from: ./genome_A.txt")
    A_gen = open("./genome_A.txt", "r").read()
    print("Reading genome B from: ./genome_B.txt")
    B_gen = open("./genome_B.txt", "r").read()
    A_min_edits = ed.compute_min_edits(A_gen, ref_gen)
    B_min_edits = ed.compute_min_edits(B_gen, ref_gen)
    
    pub_rand_str = ["19528620395885535", "ajajcwdeiqoksmzlhtes", "jha", "mahmoody", "cryptography"]
    epsilon = float(input("Enter privacy budget epsilon: ") or "0.1")
    l = int(1/(epsilon**2))
    k = 2000
    D_hat = []
    bar = Bar('Computing...     ', max=k) 
    for j in range(1,k+1):
        D = []
        for i in range(1,l+1):
            seed = get_seed(pub_rand_str[random.randint(0,4)])
            d_A = compute_d_X(A_min_edits, seed)
            d_B = compute_d_X(B_min_edits, seed)
            D_i = sec_compute_D_i(d_A, d_B)
            D.append(D_i)
        D_hat_j = sec_compute_D_hat_j(D, l)
        D_hat.append(D_hat_j)
        bar.next()
    median = statistics.median(D_hat)
    bar.finish()
    print("|A - B|: " + str(median))

def get_seed(str):
    return mmh3.hash(str, random.randint(1, 1000))

def compute_d_X(setX, seed):
    sum = 0
    for s in setX:
        hashed = mmh3.hash(s, seed)
        if(hashed % 2 == 0):
            sum += 1
        else:
            sum += -1
    return sum

##FUTURE WORK: implement secure compute
def sec_compute_D_i(d_A, d_B):
    val = (d_A - d_B)**2
    return val

##FUTURE WORK: implement secure compute
def sec_compute_D_hat_j(D, l):
    sum = 0
    for i in range(0,l):
        val = D[i] / l
        sum += val
    return sum

#don't worry about forward declaring
if __name__=="__main__":
   main()
