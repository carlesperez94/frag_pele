#!/usr/bin/env python
#title           :selectOnPlot.py
#description     :Generates a scatterplot where you can draw and select specific dots.
#author          :Carles Perez Lopez
#date            :20190219
#python_version  :3.6.5
#==============================================================================

from __future__ import absolute_import, division, print_function, unicode_literals
import os
import argparse
import glob
import multiprocessing as mp
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import shutil
from matplotlib.widgets import LassoSelector
from matplotlib.path import Path
from AdaptivePELE.atomset import RMSDCalculator, atomset
import AdaptivePELE.utilities.utilities as adapt_tools


def parseArguments():
    """
        Parse command line arguments
        :returns: str, str str, str, str, bool, int, str, str, str,    -- path to adaptive's results,
            column to plot in the X axis, column to plot in the Y axis, column to plot in the Z axis,
            path to the output's folder, whether to use a summary.csv already created, number of
            processors, report prefix, trajectory prefix, separator of csvs

    """
    desc = "Generates a scatterplot of Adaptive's results given two or three columns (X, Y, and Z if set).\n" \
           "This plot allows the selection of desired points by drawing. Structures will be selected and \n" \
           "stored into an output folder. Additionally, a report file of this selected structures will be created. \n" \
           "To be run for example like: \n" \
           "\">python selectOnPlot.py /home/usr/adaptiveresults -xcol 'Binding Energy' -ycol epoch\""
    parser = argparse.ArgumentParser(description=desc)
    required_named = parser.add_argument_group('required arguments')
    required_named.add_argument("res_path", type=str,
                                help="Path to Adaptive results.")
    parser.add_argument("-xcol", type=str, default="epoch",
                        help="Column name of the report file that will be used in the X axis.")
    parser.add_argument("-ycol", type=str, default="Binding Energy",
                        help="Column name of the report file that will be used in the Y axis.")
    parser.add_argument("-zcol", type=str, default=None,
                        help="If set, column name of the report file that will be used in the Z axis (colorbar).")
    parser.add_argument("-outfol", type=str, default=None,
                        help="If set, path to the output's folder. By default it will be created in the Adaptive's \n"
                             "results path. WARNING: Take into account that if the folder already exists it will be \n"
                             "overwritten!!!")
    parser.add_argument("-done", action="store_true",
                        help="If this is not the first time that you run this script, it is strongly recommended to \n"
                             "set this parameter on. If it is set, instead of looking all the reports and create a new\n"
                             "one, the script will use the summary csv of previous usages, saving computational time.")
    parser.add_argument("-cpus", type=int, default=4,
                        help="Number of processors that you want to use in order to save time.")
    parser.add_argument("-report", type=str, default="report_",
                        help="PELE's report prefix.")
    parser.add_argument("-traj", type=str, default="trajectory_",
                        help="Adaptive's trajectory prefix.")
    parser.add_argument("-sep", type=str, default=";",
                        help="Separator string that will be used in the CSV files.")

    args = parser.parse_args()

    return args.res_path, args.xcol, args.ycol, args.zcol, args.outfol, args.done, args.cpus, args.report, args.traj, \
           args.sep


class SelectFromCollection(object):
    """Select indices from a matplotlib collection using `LassoSelector`.

    Selected indices are saved in the `ind` attribute. This tool fades out the
    points that are not part of the selection (i.e., reduces their alpha
    values). If your collection has alpha < 1, this tool will permanently
    alter the alpha values.

    Note that this tool selects collection objects based on their *origins*
    (i.e., `offsets`).

    Parameters
    ----------
    ax : :class:`~matplotlib.axes.Axes`
        Axes to interact with.

    collection : :class:`matplotlib.collections.Collection` subclass
        Collection you want to select from.

    alpha_other : 0 <= float <= 1
        To highlight a selection, this tool sets all selected points to an
        alpha value of 1 and non-selected points to `alpha_other`.
    """

    def __init__(self, ax, collection, alpha_other=0.1):
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

        self.lasso = LassoSelector(ax, onselect=self.onselect)
        self.ind = []

    def onselect(self, verts):
        path = Path(verts)
        self.ind = np.nonzero(path.contains_points(self.xys))[0]
        self.fc[:, -1] = self.alpha_other
        self.fc[self.ind, -1] = 1
        self.collection.set_facecolors(self.fc)
        self.canvas.draw_idle()

    def disconnect(self):
        self.lasso.disconnect_events()
        self.fc[:, -1] = 1
        self.collection.set_facecolors(self.fc)
        self.canvas.draw_idle()


def concat_reports_in_csv(adaptive_results_path, output_file_path, report_prefix="report_",
                          trajectory_prefix="trajectory_", separator_out=";"):
    """
    It search report files in Adaptive's result folder and creates a csv file with everything concatenated, adding the
    epoch and trajectory information.
    :param adaptive_results_path: Path to the results folder of Adaptive.
    :type adaptive_results_path: str
    :param output_file_path: Path of the output file.
    :type output_file_path: str
    :param report_prefix: Prefix of PELE's reports.
    :type report_prefix: str
    :param trajectory_prefix: Prefix of PELE's trajectories.
    :type trajectory_prefix: str
    :param separator_out: Separator string used in the csv file.
    :type separator_out: str
    :return: Creates a csv file.
    """
    dataframe_lists = []
    for adaptive_epoch in range(0, 2000):
        folder = os.path.join(adaptive_results_path, str(adaptive_epoch))
        if os.path.exists(folder):
            report_list = glob.glob("{}/*{}*".format(folder, report_prefix))
            report_list = sorted(report_list, key=lambda x: int(x.split("_")[-1]))
            for n, report in enumerate(report_list):
                pandas_df = pd.read_csv(report, sep="    ", engine="python", index_col=False, header=0)
                pandas_df["epoch"] = adaptive_epoch
                pandas_df["trajectory"] = glob.glob("{}/{}/*{}{}.*".format(adaptive_results_path, adaptive_epoch,
                                                                          trajectory_prefix, n + 1))[0]
                dataframe_lists.append(pandas_df)
        else:
            break
    dataframe = pd.concat(dataframe_lists, ignore_index=True)
    dataframe.to_csv(output_file_path, sep=separator_out, index=False)


def trajectory_and_snapshot_to_pdb(trajectory_path, snapshot, output_path):
    """
    Given an absolute path to a trajectory of Adaptive and a snapshot (MODEL) in xtc format, the function transform it
    into a PDB format.
    :param trajectory_path: Absolute path to a trajectory from Adaptive, in xtc format.
    :type trajectory_path:str
    :param snapshot: model of a trajectory that you want to transform.
    :type snapshot: int
    :param output_path: output path of the new pdb file.
    :type output_path: str
    :return: Creates a PDB file.
    """
    topology_path_splited = trajectory_path.split("/")[0:-2]
    topology_path = os.path.join("/".join(topology_path_splited), "topology.pdb")
    topology_contents = adapt_tools.getTopologyFile(topology_path)
    trajectory = adapt_tools.getSnapshots(trajectory_path, topology=topology_path)
    try:
        single_model = trajectory[snapshot]
        PDB = atomset.PDB()
        PDB.initialise(single_model, topology=topology_contents)
    except IndexError:
        exit("You are selecting the model {} for a trajectory that has {} models, please, reselect the model index "
             "(starting from 0).".format(snapshot, len(trajectory)))
    with open(output_path, "w") as fw:
        fw.write("MODEL     %4d\n" % (snapshot + 1))
        fw.write(PDB.pdb)
        fw.write("ENDMDL\n")
        fw.write("END\n")


def get_pdb_from_xtc(row, pdbs_output_path, column_file="trajectory"):
    """
    Given a row of a dataframe (expected to come from a csv report) and a column name (that must contain the path to
    its correspondent trajectory), this function extract the file in PDB format in an output file.
    :param row: row of a dataframe (Pandas object).
    :type row: pandas.DataFrame
    :param pdbs_output_path: output path for the PDB file.
    :type pdbs_output_path: str
    :param column_file: Column name of the dataframe that contains the path to the trajectory file.
    :type column_file: str
    :return:
    """
    foldername = row[column_file]
    filepath = glob.glob(foldername)[0]
    epoch = filepath.split("/")[-2]
    snapshot = row["numberOfAcceptedPeleSteps"]
    new_file_name = os.path.basename(foldername.split("/")[-1])
    new_file_name = new_file_name.split(".")[0]
    trajectory_and_snapshot_to_pdb(filepath, snapshot, os.path.join(pdbs_output_path, "{}_epoch_{}_snap_{}.pdb".format(
        new_file_name, epoch, snapshot)
                                                                    ))
    print(os.path.join(pdbs_output_path, "{}_epoch_{}_snap_{}.pdb".format(new_file_name, epoch, snapshot)))


def get_pdbs_from_df_in_xtc(df, pdbs_output_path, processors=4, column_file="trajectory"):
    """
    It uses the function "get_pdb_from_xtc" for a whole dataframe using multiprocessing.
    :param df: Dataframe object (Pandas)
    :type df: pandas.DataFrame
    :param pdbs_output_path: Output path for PDB files.
    :type pdbs_output_path: str
    :param processors: Number of processes to do with multiprocessing.
    :type processors: int
    :param column_file: Column name of the dataframe that contains the path to the trajectory file.
    :type column_file: str
    :return:
    """
    pool = mp.Pool(processes=processors)
    multiprocessing_list = []
    for index, row in df.iterrows():
        multiprocessing_list.append(pool.apply_async(get_pdb_from_xtc,
                                                     (row, pdbs_output_path, column_file)))
    for process in multiprocessing_list:
        process.get()


def main(adaptive_results_folder,  column_to_x="epoch", column_to_y="Binding Energy", column_to_z=None,
         output_selection_folder=None, summary_done=False, processors=4, report_pref="report_",
         trajectory_pref="trajectory_", separator=";", column_file="trajectory"):
    """
    Generates a scatterplot of Adaptive's results given two or three columns (X, Y, and Z if set).
    This plot allows the selection of desired points by drawing. Structures will be selected and
    stored into an output folder. Additionally, a report file of this selected structures
    will be created.
    :param adaptive_results_folder: Path to Adaptive results.
    :type adaptive_results_folder: str
    :param column_to_x: Column name of the report file that will be used in the X axis.
    :type column_to_x: str
    :param column_to_y: Column name of the report file that will be used in the Y axis.
    :type column_to_y: str
    :param column_to_z: If set, column name of the report file that will be used in the Z axis (colorbar).
    :type column_to_z: str
    :param output_selection_folder: If set, path to the output's folder. By default it will be created in the
    Adaptive's results path. WARNING: Take into account that if the folder already exists it will be overwritten!!!
    :type output_selection_folder: str
    :param summary_done: If it is set, instead of looking all the reports and create a new one, the script will use
    the summary csv of previous usages, saving computational time."
    :type summary_done: bool
    :param processors: Number of processors that you want to use in order to save time.
    :type processors: int
    :param report_pref: PELE's report prefix.
    :type report_pref: str
    :param trajectory_pref: Adaptive's trajectory prefix.
    :type trajectory_pref: str
    :param separator: Separator string that will be used in the CSV files.
    :type separator: str
    :param column_file: Column name of the dataframe that contains the path to the trajectory file.
    :type column_file: str
    :return:
    """
    summary_csv_filename = os.path.join(adaptive_results_folder, "summary.csv")
    if not summary_done:
        concat_reports_in_csv(adaptive_results_path=adaptive_results_folder, output_file_path=summary_csv_filename,
                              report_prefix=report_pref, trajectory_prefix=trajectory_pref, separator_out=separator)
    dataframe = pd.read_csv(summary_csv_filename, sep=separator, engine='python', header=0)
    fig, ax = plt.subplots()
    if column_to_z:
        pts = ax.scatter(dataframe[column_to_x], dataframe[column_to_y], c=dataframe[column_to_z], s=20)
        plt.colorbar(pts)
    else:
        pts = ax.scatter(dataframe[column_to_x], dataframe[column_to_y], s=20)
    selector = SelectFromCollection(ax, pts)

    def accept(event, output_selection_folder=output_selection_folder):
        if event.key == "enter":
            print("Selected points:")
            df_select = dataframe.loc[selector.ind]
            print(df_select)
            counter = 0
            if not output_selection_folder:
                output_selection_folder = os.path.join(adaptive_results_folder, "selected_from_plot")
            while True:
                try:
                    os.mkdir(output_selection_folder+"_"+str(counter))
                    break
                except FileExistsError:
                    counter += 1
            output_selection_folder = output_selection_folder+"_"+str(counter)
            df_select.to_csv(os.path.join(output_selection_folder, "selection_report.csv"), sep=separator, index=False)
            get_pdbs_from_df_in_xtc(df_select, output_selection_folder, processors=processors, column_file=column_file)
            selector.disconnect()
            ax.set_title("")
            fig.canvas.draw()

    fig.canvas.mpl_connect("key_press_event", accept)
    ax.set_title("Press enter to accept selected points.")
    ax.set_xlabel(column_to_x)
    ax.set_ylabel(column_to_y)
    plt.show()


if __name__ == '__main__':
    res_path, xcol, ycol, zcol, outfol, done, cpus, report, traj, sep = parseArguments()
    main(adaptive_results_folder=res_path, column_to_x=xcol, column_to_y=ycol, column_to_z=zcol,
         output_selection_folder=outfol, summary_done=done, processors=cpus, report_pref=report,
         trajectory_pref=traj, separator=sep)
