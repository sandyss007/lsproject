import re
from collections import defaultdict, deque

#  Hybridization and SOPT extraction ---
def get_bh_data(lines):  #This function extracts bond hybridization data from log file
    results = []
    i = 0
    while i < len(lines):
        current_line = lines[i].strip()
        if "(Occupancy)   Bond orbital / Coefficients / Hybrids" in current_line:# this line searches for keywords if tht particular keywords are found then it starts extraction
            i += 1
            while i < len(lines):
                current_line = lines[i].strip()
                if "Bond orbital" in current_line and "(Occupancy)" not in current_line:
                    break
                if "CR" in current_line or "RY" in current_line:  #avoids "CR" and "RY" orbitals
                    i += 1 
                    continue
                if re.match(r"\d+\.\s+\(\d+\.\d+\)", current_line):#pattern recognition by using regex
                    results.append("\n" + current_line)
                elif re.search(r"[spd]\(\s*\d+\.\d+%\)", current_line): # pattern recognition for spd orbital and its occupancy
                    hybrid_line = current_line
                    while i + 1 < len(lines) and "f" in lines[i + 1]:  # usually f orbital is in the next line so this line recognizes it
                        i += 1 
                        hybrid_line += " " + lines[i].strip()
                    results.append(hybrid_line)#appends both spd orbital occupancy and if found f orbital occupancy too
                i += 1
        else:
            i += 1
    return results

def get_sopt_data(lines):# This function extracts second order pertubation table from the log file
    results = []
    found_section = False # assigned a False flag at the begining when the code reaches sopt data flag turns to true and starts extracting
    for line in lines:
        if "SECOND ORDER PERTURBATION THEORY" in line: # searches for this exact keyword if it is found then starts extracting
            found_section = True
            continue
        if found_section:
            if "NATURAL BOND ORBITALS" in line: # searches for this keyword when it is found stops extracting the data 
                break
            if "CR" not in line and "RY" not in line: # avoids "CR" and "RY" orbital
                results.append(line.strip())
    return results

def save_extracted_data(input_file, output_file):  # saves all the extracted data in the output file
    with open(input_file, "r", encoding="utf-8") as f:
        lines = f.readlines()

    bh_results = get_bh_data(lines)
    sopt_results = get_sopt_data(lines)

    with open(output_file, "w") as out:
        out.write("### Bond hybridization Data ###\n")
        out.write("\n".join(bh_results) + "\n\n")
        out.write("### SOPT Energy Data ###\n")
        out.write("\n".join(sopt_results))
    print("Data saved to", output_file)

# --- Integrated long-range NCI analysis ---
def analyze_long_range_nci(pdb_file, interactions_file, output_file):# takes pdb file as input and converts the connectivity data into connectivity map
    def extract_connectivity(pdb_file):
        connectivity = defaultdict(list)
        with open(pdb_file) as f:
            for line in f:
                if line.startswith("CONECT"): # identify the keyword "CONECT" to recognize connectivity data between atoms
                    parts = list(map(int, line.strip().split()[1:]))
                    atom = parts[0]
                    for connected in parts[1:]:
                        if connected not in connectivity[atom]:
                            connectivity[atom].append(connected)
                        if atom not in connectivity[connected]:
                            connectivity[connected].append(atom)
        return connectivity

    def bfs_path(connectivity, start, end): # uses bfs algorithm to identify the shortest path between two atoms 
        visited = set() # Creates a set to store atoms we've already checked and also avoids looping over the same atom again and again
        queue = deque([[start]]) # creates a queue using deque.This queue will hold possible paths (lists of atom numbers)  starting from the start atom
        while queue:
            path = queue.popleft() #Removes the first path from the queue
            node = path[-1] #Gets the last atom in the current path
            if node == end: #Checks if we’ve reached the destination atom.
                return path
            if node not in visited:
                visited.add(node) #Mark this atom as visited, so we don’t check it again later.
                for neighbor in connectivity[node]:
                    if neighbor not in path:
                        queue.append(path + [neighbor])
        return None

    def parse_interactions(file): #searches for lp-bd* and bd-bd* interactions by using regex
        interactions = []
        lp_bd_pattern = re.compile(
            r'LP\s*\(\s*\d+\).*?(\d+)\s+.*?BD\*\s*\(\s*\d+\).*?(\d+)-\s*\w+\s+\d+\s+([-+]?\d*\.\d+)'
        )
        bd_bd_pattern = re.compile(
            r'BD\s*\(\s*\d+\)\s+\w+\s+(\d+)-\s+\w+\s+(\d+)\s+.*?BD\*\s*\(\s*\d+\)\s+\w+\s+(\d+)-\s+\w+\s+\d+\s+([-+]?[0-9]*\.?[0-9]+)'
        )

        with open(file) as f: # helps to identify donor,acceptor atoms and also the E2 energy in the extracted data
            for line in f:
                line = line.strip()
                match_lp = lp_bd_pattern.search(line)
                match_bd = bd_bd_pattern.search(line)
                if match_lp:
                    donor, acceptor, energy = int(match_lp.group(1)), int(match_lp.group(2)), float(match_lp.group(3))
                    interactions.append((donor, acceptor, energy, line))
                elif match_bd:
                    donor = int(match_bd.group(2))
                    acceptor = int(match_bd.group(3))
                    energy = float(match_bd.group(4))
                    interactions.append((donor, acceptor, energy, line))
        return interactions

    connectivity = extract_connectivity(pdb_file)
    interactions = parse_interactions(interactions_file)
    seen = set()# this helps to 

    with open(output_file, "w") as out:
        for donor, acceptor, energy, original_line in interactions:  #Loops through each interaction.Skips if we’ve already seen this donor-acceptor pair.
            pair = (donor, acceptor)
            if pair in seen:
                continue
            seen.add(pair)
            if donor in connectivity and acceptor in connectivity: # identifies all the long range nci 
                path = bfs_path(connectivity, donor, acceptor)
                if path and len(path) >= 4:
                    path_str = "->".join(str(n) for n in path)
                    out.write(f"{original_line}  Path: {path_str}\n")

# --- New function to categorize long-range NCI interactions ---
def categorize_long_range_nci(file): # categories and classifies all the long range interactions
    print("\nClassifying interactions in long-range NCI output...")
    lp_pattern = re.compile(r'n-σ\*|n-π\*')
    bd_pattern = re.compile(r'σ-σ\*|σ-π\*|π-π\*|π-σ\*')

    categorized = {
        "n-σ*": [],
        "n-π*": [],
        "σ-σ*": [],
        "σ-π*": [],
        "π-π*": [],
        "π-σ*": [],
    }

    with open(file) as f:
        for line in f:
            for label in categorized:
                if label in line:
                    categorized[label].append(line.strip())
                    break

    with open(file, 'a') as out:
        out.write("\n\n### Categorized Interactions in Long-Range NCI ###\n")
        for label, lines in categorized.items():
            if lines:
                out.write(f"\n{label} Interactions\n" + "-"*50 + "\n")
                for l in lines:
                    out.write(f"{l}\n")
                out.write("-"*50 + "\n")

    print("Categorized interactions appended to", file)

# --- Interaction categorization ---
def find_lp_interactions(input_file, output_file):# this function helps to find lp interactions
    pattern = re.compile(r'LP\s*\(\s*\d+\).*BD\*\s*\(\s*(\d+)\).*([-+]?[0-9]*\.?[0-9]+)') # by using this pattern it recognizes lp interactions in the sopt 
    with open(input_file) as infile, open(output_file, 'w') as out:
        out.write("Lone Pair Interactions\n" + "-"*50 + "\n")
        for line in infile:
            match = pattern.search(line)
            if match:
                bond_type, energy = match.groups()
                label = "n-σ*" if bond_type == '1' else "n-π*"
                out.write(f"{label}: {line.strip()}\n")
        out.write("-"*50 + "\n")

def find_bd_interactions(input_file, output_file): # this function helps to find bp interactions
    pattern = re.compile(r'BD\s*\(\s*(\d+)\).*BD\*\s*\(\s*(\d+)\).*[-+]?[0-9]*\.?[0-9]+\s+([-+]?[0-9]*\.?[0-9]+)') # by using this pattern it reognizes bp interactions in the sopt
    with open(input_file) as infile, open(output_file, 'a') as out:
        out.write("\nBond Interactions\n" + "-"*50 + "\n")
        for line in infile:
            match = pattern.search(line)
            if match:
                donor, acceptor , _ = match.groups()
                if donor == '1' and acceptor == '1':  # classify the types of bp interations
                    label = "σ-σ*"
                elif donor == '1' and acceptor == '2':
                    label = "σ-π*"
                elif donor == '2' and acceptor == '2':
                    label = "π-π*"
                else:
                    label = "π-σ*"
                out.write(f"{label}: {line.strip()}\n")
        out.write("-"*50 + "\n")

# --- User Interface ---
def run_user_interface(): # asks user to  choose the type of interactions to be extracted ,and also set the threshold energy for the extraction
    print("\nChoose what to analyze:")
    print("1) lp-bp* interactions")
    print("2) bp-bp* interactions")
    choice = input("Enter your choice (1, 2 ): ").strip()

    if choice == '3':
        categorize_long_range_nci("lncP2.txt")
        return

    if choice not in ['1', '2']:
        print("Invalid choice")
        return

    if choice == '1':
        print("\nChoose interaction type:")
        print("a) n-σ* interactions")
        print("b) n-π* interactions")
        sub_choice = input("Enter a or b: ").strip().lower()
        if sub_choice == 'a':
            search_for = "n-σ*"
        elif sub_choice == 'b':
            search_for = "n-π*"
        else:
            print("Invalid choice")
            return
    else:
        print("\nChoose interaction type:")
        print("c) σ-σ* interactions")
        print("d) σ-π* interactions")
        print("e) π-π* interactions")
        print("f) π-σ* interactions")
        sub_choice = input("Enter c, d, e or f: ").strip().lower()
        if sub_choice == 'c':
            search_for = "σ-σ*"
        elif sub_choice == 'd':
            search_for = "σ-π*"
        elif sub_choice == 'e':
            search_for = "π-π*"
        elif sub_choice == 'f':
            search_for = "π-σ*"
        else:
            print("Invalid choice")
            return

    try:
        minimum_energy = float(input("\nEnter minimum energy threshold: "))
    except:
        print("Please enter a valid number")
        return

    results = []
    with open("interactions2.txt") as out:
        for line in out:
            if search_for in line:
                energy_match = re.search(r'([-+]?[0-9]*\.?[0-9]+)\s+[-+]?[0-9]*\.?[0-9]+\s+[-+]?[0-9]*\.?[0-9]+$', line) # regex pattern to identify the energy which is used in setting the threshold energy
                if energy_match and float(energy_match.group(1)) >= minimum_energy:  # asks user to set threshold energy to filter out based on energy
                    results.append(line)

    with open("final_output2.txt", "w") as out:
        out.writelines(results)
    print("\nResults saved to final_output2.txt")

# --- Main Program ---
input_file = input("Enter your NBO log file name: ").strip()
save_extracted_data(input_file, "extracted_data2.txt")
find_lp_interactions("extracted_data2.txt", "interactions2.txt")
find_bd_interactions("extracted_data2.txt", "interactions2.txt")

analyze_long_range_nci("N2.pdb", "interactions2.txt", "lncP2.txt")
run_user_interface()

