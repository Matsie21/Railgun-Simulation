import matplotlib.pyplot as plt
import numpy as np
from deserialize import ExportReader

# Convert our reader into something matplotlib understands
def get_np_arrays(reader, x, y):
    xlist = []
    ylist = []

    point = reader.read_point()
    while point is not None:
        xlist.append(x(point)) # Store the requested value on the x-axis
        ylist.append(y(point)) # Store the requested value on the y-axis
        point = reader.read_point()

    return np.array(xlist), np.array(ylist)

# Helper function to generate a plot based on two functions: one for the x-axis and one for the y-axis
# The functions get a datapoint as input and must output a number
# This allows us to easily plot any pair of two values on a simple graph
def show_plot(reader, x, y):
    xpoints, ypoints = get_np_arrays(reader, x, y) # Convert our points into something matplotlib understands
    plt.plot(xpoints, ypoints) # Generate the plot with matplotlib
    plt.show() # Show the plot

def main():
    # Increase the resolution of all graphs generated my matplotlib as the default resolution is too low
    plt.rcParams['figure.dpi'] = 300

    # Create our exportreader and tell it that the export file is located at specified path
    reader = ExportReader(r"..")
    # Generate a plot with data from reader
    show_plot(
        reader,
        lambda point: point.t, # x-axis
        lambda point: point.I # y-axis
    )

if __name__ == "__main__":
    main()
