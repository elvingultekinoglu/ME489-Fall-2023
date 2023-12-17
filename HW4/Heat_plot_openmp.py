import matplotlib.pyplot as plt

# Read data from x_data.txt
with open('openmp_x_data.txt', 'r') as file:
    x_data_lines = [line.strip().split() for line in file]

# Flatten the list and convert values to float
x_data = [float(value) for line in x_data_lines for value in line]

# Read data from q_data.txt
with open('openmp_q_data.txt', 'r') as file:
    q_data_lines = [line.strip().split() for line in file]

# Flatten the list and convert values to float for the y-axis data
q_data = [float(value) for line in q_data_lines for value in line]

# Plotting the data
plt.plot(x_data, q_data)

# Adding labels and title
plt.xlabel('X-axis Label')
plt.ylabel('Y-axis Label')
plt.title('Data Plot from x_data.txt and q_data.txt')

# Show the plot
plt.show()
