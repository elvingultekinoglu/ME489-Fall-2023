import sys

def update_dat_file(nx_value, ny_value):
    file_path = "input.dat"

    # Read the content of the file
    with open(file_path, 'r') as file:
        lines = file.readlines()

    # Find the indices of [NX] and [NY] lines
    nx_index = lines.index('[NX]\n')
    ny_index = lines.index('[NY]\n')

    # Update the values based on command-line arguments
    lines[nx_index + 1] = str(nx_value) + '\n'
    lines[ny_index + 1] = str(ny_value) + '\n'

    # Write the modified content back to the file
    with open(file_path, 'w') as file:
        file.writelines(lines)

if __name__ == "__main__":
    # Check if command-line arguments are provided
    if len(sys.argv) != 2:
        print("Usage: python3 updated_nx.py <new_NX_value>")
        sys.exit(1)

    try:
        # Get the new value for [NX] from the command-line argument
        nx_value = int(sys.argv[1])

        # Update the .dat file
        update_dat_file(nx_value, nx_value)

        #print("File updated successfully.")
    except ValueError:
        print("Invalid input. Please provide a valid integer for [NX].")
        sys.exit(1)
