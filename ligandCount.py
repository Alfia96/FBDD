def count_ligands(filename):
  """
  This function counts the number of ligands in a file where each ligand name starts with "zinc" and ends with "end".

  Args:
      filename (str): The name of the file containing ligand names.

  Returns:
      int: The total number of ligands in the file.
  """

  # Initialize counter
  ligand_count = 0

  # Open the file in read mode
  with open(filename, 'r') as file:
    # Read the file line by line
    for line in file:
      #print(line, end="")
      # Check if the line starts with "zinc" and ends with "end" (ignoring case)
      if "zinc".lower() in line.lower():
        ligand_count += 1
        #print(f"Found ligand: {line}")

  # Return the total ligand count
  return ligand_count

# Specify the filename containing the merged ligand names
filename = "final.sdf"  # Replace with your actual filename

# Call the function to count ligands
total_ligands = count_ligands(filename)

# Print the results
print(f"Total number of ligands found: {total_ligands}")

