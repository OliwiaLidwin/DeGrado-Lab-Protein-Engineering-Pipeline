import os
import re
import matplotlib.pyplot as plt

# 
# Renames the file using an extracted value.
# Creates a histogram of directory data.
#
# Change the definition of extract_number_from_filename to your desired naming.
# Change the range of data you want your histogram to cover in a different color. (Ex. 0.15 to 0.74 for resname 38E)
# Change line 67 accoringly for the directory.
# 

def extract_number_from_filename(filename):
    """Extracts the decimal number from the filename."""
    match = re.search(r'_no_CG_(\d+\.\d+)', filename)
    return float(match.group(1)) if match else None

def truncate_filename(filename):
    """Truncates the filename to exclude everything after '_no_CG_'."""
    return filename.split('_no_CG_')[0]

def plot_bar_chart(directory):
    """Reads files from the directory, creates a bar chart, and highlights specific ranges."""
    filenames = os.listdir(directory)
    numbers = []

    for filename in filenames:
        number = extract_number_from_filename(filename)
        if number is not None:
            truncated_name = truncate_filename(filename)
            numbers.append((truncated_name, number))

    # Sort the numbers by the extracted distance
    numbers.sort(key=lambda x: x[1])

    # Create a list of colors based on the range
    colors = ['red' if 0.15 <= num <= 0.74 else 'blue' for _, num in numbers]
    truncated_names = [name for name, _ in numbers]
    distances = [num for _, num in numbers]

    # Plot bar chart
    plt.figure(figsize=(12, 8))
    
    bars = plt.bar(truncated_names, distances, color=colors, alpha=0.7)
    
    # Annotate the bars with the distances
    for bar, distance in zip(bars, distances):
        plt.text(bar.get_x() + bar.get_width() / 2, bar.get_height(),
                 f'{distance:.4f}', ha='center', va='bottom', fontsize=8, rotation=90)
    
    # Labeling the axes
    plt.xlabel('Filename')
    plt.ylabel('Distance to Plane (Angstrom)')
    plt.title('Serine Side Chain Hydroxyl Distance to Plane')
    plt.xticks(rotation=90)  # Rotate x-axis labels for better readability

    # Create custom legend
    from matplotlib.patches import Patch
    legend_elements = [Patch(facecolor='blue', edgecolor='blue', label='Outside Range (0.15 - 0.74)'),
                       Patch(facecolor='red', edgecolor='red', label='Within Range (0.15 - 0.74)')]
    plt.legend(handles=legend_elements)

    # Save and show the plot
    plt.tight_layout()
    plt.savefig('bar_chart_COMBS1.pdf')
    plt.show()

# Define the directory containing the files
directory = 'C:/Users/oliwi/Downloads/All_Labeled_w_Distance'

# Plot the bar chart
plot_bar_chart(directory)
