def main():
    print("Reading ref genome from: ./genome_ref.txt")
    ref_gen = open("./genome_ref.txt", "r").read()
    print("Reading genome A from: ./genome_A.txt")
    A_gen = open("./genome_A.txt", "r").read()
    print("Reading genome B from: ./genome_B.txt")
    B_gen = open("./genome_B.txt", "r").read()
    B_min_edits = compute_min_edits(B_gen, ref_gen)
    A_min_edits = compute_min_edits(A_gen, ref_gen)
    A_B_diff = compute_set_diff(A_min_edits, B_min_edits)
    print("|A - B|: " + str(len(A_B_diff)))

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

def compute_set_diff(set_A, set_B):
    set_diff = []
    for a in set_A:
        if a not in set_B:
            set_diff.append(a)
    for b in set_B:
        if b not in set_A:
            set_diff.append(b)

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

