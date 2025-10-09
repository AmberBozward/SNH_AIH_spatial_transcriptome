
import numpy as np
import pandas as pd
import matplotlib
import time

from matplotlib.path import Path
from matplotlib.widgets import PolygonSelector


class SelectFromCollection:
    """
    Select indices from a matplotlib collection using `PolygonSelector`.

    Selected indices are saved in the `ind` attribute. This tool fades out the
    points that are not part of the selection (i.e., reduces their alpha
    values). If your collection has alpha < 1, this tool will permanently
    alter the alpha values.

    Note that this tool selects collection objects based on their *origins*
    (i.e., `offsets`).

    Parameters
    ----------
    ax : `~matplotlib.axes.Axes`
        Axes to interact with.
    collection : `matplotlib.collections.Collection` subclass
        Collection you want to select from.
    alpha_other : 0 <= float <= 1
        To highlight a selection, this tool sets all selected points to an
        alpha value of 1 and non-selected points to *alpha_other*.
    """

    def __init__(self, ax, collection, alpha_other=0.3):
        self.canvas = ax.figure.canvas
        self.collection = collection
        self.alpha_other = alpha_other

        self.xys = collection.get_offsets()
        self.Npts = len(self.xys)

        # Ensure that we have separate colors for each object
        self.fc = collection.get_facecolors()
        if len(self.fc) == 0:
            raise ValueError('Collection must have a facecolor')
        elif len(self.fc) == 1:
            self.fc = np.tile(self.fc, (self.Npts, 1))

        self.poly = PolygonSelector(ax, self.onselect, draw_bounding_box=True)
        self.ind = []

    def onselect(self, verts):
        path = Path(verts)
        self.ind = np.nonzero(path.contains_points(self.xys))[0]
        self.fc[:, -1] = self.alpha_other
        self.fc[self.ind, -1] = 1
        self.collection.set_facecolors(self.fc)
        self.canvas.draw_idle()

    def disconnect(self):
        self.poly.disconnect_events()
        self.fc[:, -1] = 1
        self.collection.set_facecolors(self.fc)
        self.canvas.draw_idle()


if __name__ == '__main__':
    import matplotlib.pyplot as plt

    # read file
    df = pd.read_csv('/Users/bozwarag/Desktop/centroidsss.csv', sep='\t')

    # paramters
    sample = "Donor_6446"
    fov_min = 0
    fov_max = 1000
    time_stamp = time.time()

    # subset on sample
    df = df[df["SAMPLE"] == sample]

    # subset on FOV
    df = df[df["FOV"] >= fov_min]
    df = df[df["FOV"] <= fov_max]

    # filter columns
    df = df.drop('SAMPLE', axis=1)

    # get x, y and cell
    grid_x = np.array(df["X"])
    grid_y = np.array(df["Y"])
    cell = np.array(df["CELL"])

    # colour map
    cmap = matplotlib.colors.ListedColormap(["#88CCEE", "#CC6677", "#DDCC77", "#117733", "#332288", "#AA4499", "#44AA99", "#999933", "#882255", "#661100", "#6699CC", "#888888"], name='from_list', N=None)

    # plot
    fig, ax = plt.subplots()
    pts = ax.scatter(grid_x, grid_y, c = cell, s=2, cmap=cmap, )
    ax.set_facecolor('grey')
    selector = SelectFromCollection(ax, pts)

    # instructions
    print("Select points in the figure by enclosing them within a polygon.")
    print("Press the 'esc' key to start a new polygon.")
    print("Try holding the 'shift' key to move all of the vertices.")
    print("Try holding the 'ctrl' key to move a single vertex.")

    # show interactive selector
    plt.show()

    # to disconnect
    selector.disconnect()

    # After figure is closed print the coordinates of the selected points
    print('\nSelected points:')
    print(selector.xys[selector.ind])

    # get the points
    points = selector.xys[selector.ind]

    # put the cells into a dictionary by coordiantes
    cells_dict = {}

    for i in range(1,len(df)):

        sub_df = df.iloc[[i]]
        X = sub_df.iloc[0]['X']
        Y = sub_df.iloc[0]['Y']
        cell = sub_df.index[0]

        cells_dict[str(X) + "_" + str(Y)] = cell


    # interrogate dictionary and save cells
    out_file = open("C:/Users/JohnJ/Desktop/analysis/Ye_Oo/Amber_B_CosMx_2024/analysis/combined_all/polygons/" + sample + "_" + str(fov_min) + "_" + str(fov_max) + "_" + str(time_stamp) + ".txt", "w")
    for point in points:
        try:
            cell = cells_dict[str(point[0]) + "_" + str(point[1])]
            out_file.write(cell + "\n")
        except:
            continue

    out_file.flush()


