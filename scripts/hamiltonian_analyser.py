
file_name = "output.txt"

def main():
    with open(file_name, 'r') as f:
        lines = f.readlines()
    
    # hamiltonian = list()
    hamiltonian_0 = 0
    for i, line in enumerate(lines):
        elements = line.split()
        if i == 0:
            continue
        # print(elements)
        one_hamiltonian = float(elements[3])
        if i == 1:
            hamiltonian_0 = one_hamiltonian
        numerator = abs(one_hamiltonian - hamiltonian_0)
        relative_drift = numerator / abs(hamiltonian_0)
        print(f"Hamiltonian: {one_hamiltonian:.6f}, Relative Drift: {relative_drift:.6e}")

        # hamiltonian.append(float(elements[3]))
    # print(hamiltonian)
    


main()
        