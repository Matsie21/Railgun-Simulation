import matplotlib.pyplot as plt
import numpy as np
from deserialize import ExportReader

def get_np_arrays(reader, x, y):
    xlist = []
    ylist = []

    point = reader.read_point()
    while point is not None:
        xlist.append(x(point))
        ylist.append(y(point))
        point = reader.read_point()

    return np.array(xlist), np.array(ylist)

def show_plot(reader, x, y):
    xpoints, ypoints = get_np_arrays(reader, x, y)
    plt.plot(xpoints, ypoints)
    plt.show()

def main():
    reader = ExportReader("../cpp_simulation/buildDir/out.bin")
    show_plot(
        reader,
        lambda point: point.t,
        lambda point: point.loc_x
    )

if __name__ == "__main__":
    main()