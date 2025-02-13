import math

def calculate_distance(a1, a2): return math.sqrt(sum((b - a) ** 2 for a, b in zip(a1, a2)))

def calculate_angle(a1, a2, a3):
    v1, v2 = [x1 - x2 for x1, x2 in zip(a1, a2)], [x3 - x2 for x3, x2 in zip(a3, a2)]
    return math.degrees(math.acos(max(-1.0, min(1.0, sum(x * y for x, y in zip(v1, v2)) / (math.sqrt(sum(x ** 2 for x in v1)) * math.sqrt(sum(x ** 2 for x in v2)))))))

def calculate_dihedral(a1, a2, a3, a4):
    def vector_subtract(v1, v2):
        return [x1 - x2 for x1, x2 in zip(v1, v2)]
    def cross_product(v1, v2):
        return [
            v1[1] * v2[2] - v1[2] * v2[1],
            v1[2] * v2[0] - v1[0] * v2[2],
            v1[0] * v2[1] - v1[1] * v2[0]
        ]
    def dot_product(v1, v2):
        return sum(x * y for x, y in zip(v1, v2))
    def magnitude(v):
        return math.sqrt(sum(x ** 2 for x in v))
    # Compute vectors between points
    v1 = vector_subtract(a2, a1)
    v2 = vector_subtract(a3, a2)
    v3 = vector_subtract(a4, a3)
    # Normal vectors to planes
    n1 = cross_product(v1, v2)
    n2 = cross_product(v2, v3)
    # Calculate dihedral angle
    cos_theta = dot_product(n1, n2) / (magnitude(n1) * magnitude(n2))
    
    return math.degrees(math.acos(cos_theta))
           
def parse_pdb_file(pdb_file):
    with open(pdb_file) as f:
        return {f"{line[12:16].strip()}{line[6:11].strip()}": tuple(map(float, [line[30:38], line[38:46], line[46:54]])) for line in f if line.startswith(("ATOM", "HETATM"))}

def main():
    pdb_file, calc_file = "/home/ibab/project/DATA_FROM_JN/real_files/N2.pdb", "/home/ibab/project/DATA_FROM_JN/real_files/1.txt"
    atoms = parse_pdb_file(pdb_file)
    for calc in (line.split() for line in open(calc_file) if line.strip()):
        try:
            if calc[0] == "bond":
                print(f"Bond {calc[1]}-{calc[2]}: {calculate_distance(atoms[calc[1]], atoms[calc[2]]):.4f} Å")
            elif calc[0] == "angle":
                print(f"Angle {calc[1]}-{calc[2]}-{calc[3]}: {calculate_angle(atoms[calc[1]], atoms[calc[2]], atoms[calc[3]]):.2f}°")
            elif calc[0] == "dihedral":
                print(f"Dihedral {calc[1]}-{calc[2]}-{calc[3]}-{calc[4]}: {calculate_dihedral(atoms[calc[1]], atoms[calc[2]], atoms[calc[3]], atoms[calc[4]]):.2f}°")
        except KeyError as e:
            print(f"Error: Atom {e.args[0]} not found.")

if __name__ == "__main__":
    main()





