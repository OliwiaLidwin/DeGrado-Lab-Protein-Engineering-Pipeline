from pymol import cmd
import numpy as np

# Use the pymol GUI to visualize a plane made from three atoms.
# Edit line 66 to change to the elements of interest.

def plane_from_points(a1, a2, a3):
    """
    Calculate and display a plane from three points.
    """
    coords = []
    for atom in [a1, a2, a3]:
        model = cmd.get_model(atom)
        if len(model.atom) == 0:
            print(f"Error: Atom selection '{atom}' returned no atoms.")
            return
        coords.append(model.atom[0].coord)  # Extract coordinates of the first atom in the model

    if len(coords) != 3:
        print("Error: Selection should return exactly three atoms.")
        return

    coords = np.array(coords)
    p1, p2, p3 = coords

    # Calculate the plane normal vector
    v1 = p2 - p1
    v2 = p3 - p1
    normal = np.cross(v1, v2)

    # Calculate the plane equation: ax + by + cz = d
    a, b, c = normal
    d = np.dot(normal, p1)

    # Define grid size and resolution
    size = 20.0  # Extended size
    resolution = 1.0  # Adjust resolution as needed

    # Generate grid points
    x = np.arange(-size, size, resolution)
    y = np.arange(-size, size, resolution)
    x, y = np.meshgrid(x, y)

    # Calculate corresponding z values
    z = (d - a*x - b*y) / c

    # Flatten arrays for PyMOL
    x = x.flatten()
    y = y.flatten()
    z = z.flatten()

    # Create a temporary object to hold the plane grid
    cmd.delete('plane_grid')
    cmd.create('plane_grid', 'all')  # Create a placeholder object

    # Add points to the object
    for xi, yi, zi in zip(x, y, z):
        cmd.pseudoatom(object='plane_grid', pos=(xi, yi, zi))

    # Display the grid as spheres
    cmd.show_as('spheres', 'plane_grid')
    cmd.color('blue', 'plane_grid')

# Example usage
plane_from_points('resname 38E and name N1', 'resname 38E and name C7', 'resname 38E and name C6')