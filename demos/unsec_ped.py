import edit_distance as ed

def main():
    print("Reading ref genome from: ./genome_ref.txt")
    ref_gen = open("./genome_ref.txt", "r").read()
    print("Reading genome A from: ./genome_A.txt")
    A_gen = open("./genome_A.txt", "r").read()
    print("Reading genome B from: ./genome_B.txt")
    B_gen = open("./genome_B.txt", "r").read()
    B_min_edits = ed.compute_min_edits(B_gen, ref_gen)
    A_min_edits = ed.compute_min_edits(A_gen, ref_gen)
    A_B_diff = ed.compute_edit_set_diff(A_min_edits, B_min_edits)
    print("|A - B|: " + str(len(A_B_diff)))
    
#don't worry about forward declaring
if __name__=="__main__":
   main()
