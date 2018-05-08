import mmh3
import statistics
import random
from progress.bar import Bar
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
    A_min_edits = compute_min_edits(A_gen, ref_gen)
    B_min_edits = compute_min_edits(B_gen, ref_gen)
    
    pub_rand_str = ["19528620395885535", "ajajcwdeiqoksmzlhtes", "jha", "mahmoody", "cryptography"]
    epsilon = float(input("Enter privacy budget ε: ") or "0.1")
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

##TODO: implement secure compute
def sec_compute_D_i(d_A, d_B):
    val = (d_A - d_B)**2
    return val

##TODO: implement secure compute
def sec_compute_D_hat_j(D, l):
    sum = 0
    for i in range(0,l):
        val = D[i] / l
        sum += val
    return sum
##########################################################################################################
##########################################################################################################
def compute_min_edits(ref_str, trans_str):
    leven = levenshtein(ref_str, trans_str)
    edits = leven[1]
    parsed = parse_edits(edits, ref_str, trans_str)
    return parsed
    
def parse_edits(edits, ref_str, trans_str):
    parsed = []
    for line in edits:
        op = line['type']
        loc = int(line['i']) + 1
        trans_loc = int(line['j']) - 1
        change = trans_str[trans_loc]
        if(op == "deletion"):
            parsed.append(str(loc) + op[0:3] + '')
        elif(op == "insertion" or  op == "substitution"):
            parsed.append(str(loc) + op[0:3] + change)
    return parsed

def compute_edit_set_diff(A_edits, B_edits):
    set_diff = []
    for x in A_edits:
        if x not in B_edits:
            set_diff.append(x)
    for x in B_edits:
        if x not in A_edits:
            set_diff.append(x)

    return set_diff
### Code taken from (Curzona's Levenshtein Gist)[https://gist.github.com/curzona/9435822] and slightly modified for readability ########

# Calculates the levenshtein distance and the edits between two strings
def levenshtein(ref_str, op_str, key=hash):
  rows = costmatrix(ref_str, op_str, key)
  edits = backtrace(ref_str, op_str, rows, key)
 
  return rows[-1][-1], edits
 
# Generate the cost matrix for the two strings
# Based on http://en.wikibooks.org/wiki/Algorithm_Implementation/Strings/Levenshtein_distance#Python
def costmatrix(ref_str, op_str, key=hash):
  rows = []
 
  previous_row = range(len(op_str) + 1)
  rows.append(list(previous_row))
 
  for i, c1 in enumerate(ref_str):
    current_row = [i + 1]
    for j, c2 in enumerate(op_str):
      insertions = previous_row[j + 1] + 1
      deletions = current_row[j] + 1
      substitutions = previous_row[j] + (key(c1) != key(c2))
      current_row.append(min(insertions, deletions, substitutions))
    previous_row = current_row
 
    rows.append(previous_row)
 
  return rows
 
# Trace back through the cost matrix to generate the list of edits
def backtrace(ref_str, op_str, rows, key=hash):
  i, j = len(ref_str), len(op_str)
 
  edits = []
 
  while(not (i == 0  and j == 0)):
    prev_cost = rows[i][j]
 
    neighbors = []
 
    if(i!=0 and j!=0):
      neighbors.append(rows[i-1][j-1])
    if(i!=0):
      neighbors.append(rows[i-1][j])
    if(j!=0):
      neighbors.append(rows[i][j-1])
 
    min_cost = min(neighbors)
 
    if(min_cost == prev_cost):
      i, j = i-1, j-1
      edits.append({'type':'match', 'i':i, 'j':j})
    elif(i!=0 and j!=0 and min_cost == rows[i-1][j-1]):
      i, j = i-1, j-1
      edits.append({'type':'substitution', 'i':i, 'j':j})
    elif(i!=0 and min_cost == rows[i-1][j]):
      i, j = i-1, j
      edits.append({'type':'deletion', 'i':i, 'j':j})
    elif(j!=0 and min_cost == rows[i][j-1]):
      i, j = i, j-1
      edits.append({'type':'insertion', 'i':i, 'j':j})
 
  edits.reverse()
 
  return edits

##########################################################################################################

#don't worry about forward declaring
if __name__=="__main__":
   main()
