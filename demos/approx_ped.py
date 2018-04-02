def main():
    ref_gen = input("Enter Reference Genome or enter nothing for demo: \n") or "ATTGCCCGA"

    A_gen = input("Enter Genome A or enter nothing for demo: \n") or "GTTGGATAA"
    A_min_edits = compute_min_edits(A_gen, ref_gen)

    B_gen = input("Enter Genome B or enter nothing for demo: \n") or "GTTCGATGA"
    B_min_edits = compute_min_edits(B_gen, ref_gen)
    
    print("Reference Genome: " + ref_gen)
    print("Genome A: " + A_gen)
    print("Genome B: " + B_gen)
    print("Min edits from Ref -> A: ")
    print(A_min_edits)
    print("Min edits from Ref -> B:")
    print(B_min_edits)

    A_B_diff = compute_edit_set_diff(A_min_edits, B_min_edits)
    print("Set Difference: " + str(A_B_diff))
    print("With Cardinality: " + str(len(A_B_diff)))

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
            parsed.append((loc, op[0:3], ''))
        elif(op == "insertion" or  op == "substitution"):
            parsed.append((loc, op[0:3], change))
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

