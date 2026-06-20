import matplotlib.pyplot as plt


file_name = "output.txt"

def main():
    with open(file_name, 'r') as f:
        lines = f.readlines()
    
    hamiltonian = []
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
        print(f"step {i}: Hamiltonian: {one_hamiltonian:.6f}, Relative Drift: {relative_drift:.6e}")

        hamiltonian.append(one_hamiltonian)   
    plt.plot(hamiltonian, label='Hamiltonian')
    plt.xlabel('Time Step')
    plt.ylabel('Hamiltonian')
    plt.title('Hamiltonian vs Time Step')
    plt.legend()
    plt.grid()
    plt.savefig('my_plot.png') 

    # print(hamiltonian)
    


main()
        