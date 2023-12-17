import sys
import matplotlib.pyplot as plt
import numpy as np


# Function to read data from a file
def read_data(filename):
    with open(filename, 'r') as file:
        lines = file.readlines()
        n = len(lines[0].split())
        data = [float(value) for value in lines[0].split()]
    return n, np.array(data)



# Check if a command-line argument (step number) is provided
if len(sys.argv) != 2:
    print("Usage: python script.py <size> ")
    sys.exit(1)

try:
    num_ranks = int(sys.argv[1])

except ValueError:
    print("Invalid step number. Please provide a valid integer.")
    sys.exit(1)

# Read 'n' and data from the first file and update the arrays
x_filename = f"mpi_x_data_size{num_ranks}_rank0.txt"
n, x_data_0 = read_data(x_filename)

# Initialize empty arrays to store data
x_data = np.zeros((num_ranks, n))
q_data = np.zeros((num_ranks, n))

# Loop through ranks to read and store data for the specified step
for rank in range(num_ranks):
    x_filename = f"mpi_x_data_size{num_ranks}_rank{rank}.txt"
    q_filename = f"mpi_q_data_size{num_ranks}_rank{rank}.txt"

    dummy,x_data[rank, :] = read_data(x_filename)
    dummy,q_data[rank, :] = read_data(q_filename)

# Plotting
plt.figure(figsize=(10, 6))
for rank in sorted(range(num_ranks), reverse=True):
    plt.plot(x_data[rank, :], q_data[rank, :], label=f'Core {rank}')

plt.title(f'Q Data vs X Data')
plt.xlabel('X Value')
plt.ylabel('Q Value')
plt.legend()
plt.show()