#/usr/bin/python2.7

import matplotlib
matplotlib.use('TkAgg')
import matplotlib.pyplot as plt
from matplotlib.widgets  import RectangleSelector
import numpy as np
import os
import errno
import argparse
import pandas as pd
import glob
import re
import sys
import warnings
import AdaptivePELE.analysis.splitTrajectory as st

"""

   Description: Parse all the reports found under 'path' and sort them all
   by the chosen criteria (Binding Energy as default) having into account the
   frequency our pele control file writes a structure through the -ofreq param
   (1 by default). To sort from higher to lower value use -f "max" otherwise
   will rank the structures from lower to higher criteria's values. The number
   of structures will be ranked is controlled by -i 'nstruct' (default 10).

   For any problem do not hesitate to contact us through the email address written below.

"""

__author__ = "Daniel Soler Viladrich"
__email__ = "daniel.soler@nostrumbiodiscovery.com"

# DEFAULT VALUES
ORDER = "min"
CRITERIA = ["Binding", "Energy"]
OUTPUT = "Structure_{}.pdb"
FREQ = 1
REPORT = "report"
TRAJ = "trajectory"
ACCEPTED_STEPS = 'numberOfAcceptedPeleSteps'
OUTPUT_FOLDER = 'BestStructs'
DIR = os.path.abspath(os.getcwd())
STEPS=3


def parse_args():

    parser = argparse.ArgumentParser()
    parser.add_argument("crit1", type=int, help="First Criteria we want to rank and output the strutures for. Must be a column of the report. i.e: Binding Energy")
    parser.add_argument("crit2", type=int, help="Second Criteria we want to rank and output the strutures for. Must be a column of the report. i.e: Binding Energy")
    parser.add_argument("ad_steps", type=int, help="Adaptive Steps i.e: self.ad_steps")
    parser.add_argument("--path", type=str, help="Path to Pele's results root folder i.e: path=/Pele/results/", default=DIR)
    parser.add_argument("--ofreq", "-f", type=int, help="Every how many steps the trajectory were outputted on PELE i.e: self.ad_steps", default=FREQ)
    parser.add_argument("--out", "-o", type=str, help="Output Path. i.e: BindingEnergies_apo", default=OUTPUT_FOLDER)
    parser.add_argument("--numfolders", "-nm", action="store_true", help="Not to parse non numerical folders")
    parser.add_argument("--top", "-t", type=str, help="Topology file for xtc", default=None)
    args = parser.parse_args()

    return args.crit1, args.crit2, args.ad_steps, os.path.abspath(args.path), args.ofreq, args.out, args.numfolders, args.top

def is_adaptive():
  folders = glob.glob("{}/*/".format(DIR))
  folders_numerical = [os.path.basename(os.path.normpath(folder)) for folder in folders]
  if len(folders_numerical) > 1:
    return True
  else:
    return False


class DataHandler(object):

  def __init__(self, metrics, crit1, crit2, index1, index2, steps, adaptive, ad_steps, axis):
    self.metrics = metrics
    self.crit1 = crit1
    self.crit2 = crit2
    self.index1 = index1
    self.index2 = index2
    self.steps = steps
    self.adaptive = adaptive
    self.ad_steps = ad_steps
    self.axis = axis
    self.descompose_values()


  def descompose_values(self):
    self.paths = self.metrics[DIR].tolist()
    self.epochs = [os.path.basename(os.path.normpath(os.path.dirname(Path))) for Path in self.paths]
    self.values1, self.values2 = self.retrieve_values()
    self.file_ids = self.metrics.report.tolist()
    self.step_indexes = self.metrics[self.steps].tolist()

  def on_press(self, event):
    self.limits_start = [self.axis.get_xlim(), self.axis.get_ylim()]
    if not event.inaxes: return
    self.xo, self.yo = event.xdata, event.ydata

   
  def on_release(self, event):
    self.limits_end = [self.axis.get_xlim(), self.axis.get_ylim()]
    if not event.inaxes: return
    self.xf, self.yf = event.xdata, event.ydata
    if self.limits_start == self.limits_end:
      self.compute()

  def compute(self):
    self.retrieve_data()
    if topology:
        self.extract_snapshots_from_xtc(self.data_to_extract, self.steps)
    else:
        self.extract_snapshots_from_pdb(self.data_to_extract, self.steps)

  def retrieve_data(self):
    if (self.xf > self.xo) and (self.yf < self.yo):
      self.data_to_extract = (self.metrics[(self.metrics[self.crit2] > self.yf) & (self.metrics[self.crit2] < self.yo) &
      (self.metrics[self.crit1] > self.xo) & (self.metrics[self.crit1] < self.xf) ])

    elif (self.xf > self.xo) and (self.yf > self.yo):
      self.data_to_extract = (self.metrics[(self.metrics[self.crit2] < self.yf) & (self.metrics[self.crit2] > self.yo) &
      (self.metrics[self.crit1] > self.xo) & (self.metrics[self.crit1] < self.xf) ])

    elif (self.xf < self.xo) and (self.yf > self.yo):
      self.data_to_extract = (self.metrics[(self.metrics[self.crit2] < self.yf) & (self.metrics[self.crit2] > self.yo) &
      (self.metrics[self.crit1] < self.xo) & (self.metrics[self.crit1] > self.xf) ])

    elif (self.xf < self.xo) and (self.yf < self.yo):
      self.data_to_extract = (self.metrics[(self.metrics[self.crit2] > self.yf) & (self.metrics[self.crit2] < self.yo) &
      (self.metrics[self.crit1] < self.xo) & (self.metrics[self.crit1] > self.xf) ])
    else:
      self.data_to_extract = pd.DataFrame(columns=list(self.metrics))


  def extract_snapshots_from_pdb(self, min_values, steps):
      paths = min_values[DIR].tolist()
      epochs = [os.path.basename(os.path.normpath(os.path.dirname(Path))) for Path in paths]
      values1 = min_values[self.crit1].tolist()
      values2 = min_values[self.crit2].tolist()
      file_ids = min_values.report.tolist()
      step_indexes = min_values[steps].tolist()
      files_out = ["epoch{}_trajectory_{}.{}_{}{:.2f}_{}{:.3f}.pdb".format(epoch, report, int(step), self.crit1.replace(" ",""),
        value1, self.crit2.replace(" ",""), value2) \
        for epoch, step, report, value1, value2 in zip(epochs, step_indexes, file_ids, values1, values2)]
      for f_id, f_out, step, path in zip(file_ids, files_out, step_indexes, paths):

          # Read Trajetory from PELE's output
          f_in = glob.glob(os.path.join(os.path.dirname(path), "*trajectory*_{}.pdb".format(f_id)))
          if len(f_in) == 0:
              sys.exit("Trajectory {} not found. Be aware that PELE trajectories must contain the label \'trajectory\' in their file name to be detected".format("*trajectory*_{}".format(f_id)))
          f_in = f_in[0]

          with open(f_in, 'r') as input_file:
              file_content = input_file.read()

          if self.adaptive and (self.steps in [self.crit1, self.crit2]):
            model = (step % self.ad_steps)/out_freq+1
          else:
            model = (step)/out_freq+1

          trajectory_selected = re.search('MODEL\s+%d(.*?)ENDMDL' %int(model), file_content, re.DOTALL)

          # Output Trajectory
          try:
              mkdir_p(output)
          except OSError:
              pass

          traj = []
          with open(os.path.join(output,f_out),'w') as f:
              traj.append("MODEL     %d" %int(model))
              try:
                  traj.append(trajectory_selected.group(1))
              except AttributeError:
                  raise AttributeError("Model not found. Check the -f option.")
              traj.append("ENDMDL\n")
              f.write("\n".join(traj))
          print("MODEL {} has been selected".format(f_out))

  def extract_snapshots_from_xtc(self, min_values, steps):
    paths = min_values[DIR].tolist()
    epochs = [os.path.basename(os.path.normpath(os.path.dirname(Path))) for Path in paths]
    values1 = min_values[self.crit1].tolist()
    values2 = min_values[self.crit2].tolist()
    file_ids = min_values.report.tolist()
    step_indexes = min_values[steps].tolist()
    files_out = ["epoch{}_trajectory_{}.{}_{}{:.2f}_{}{:.3f}.pdb".format(epoch, report, int(step), self.crit1.replace(" ",""),
        value1, self.crit2.replace(" ",""), value2) \
        for epoch, step, report, value1, value2 in zip(epochs, step_indexes, file_ids, values1, values2)]
    for f_id, f_out, step, path in zip(file_ids, files_out, step_indexes, paths):
        f_in = glob.glob(os.path.join(os.path.dirname(path), "*trajectory*_{}.xtc".format(f_id)))
        found = st.main(output, f_in, topology, [step % self.ad_steps/out_freq+1], template=f_out)
        if found:
            print("MODEL {} has been selected".format(f_out))
        else:
            print("MODEL {} not found. Check -f option".format(f_out))


  def retrieve_values(self):

    if (self.steps == self.crit1) & self.adaptive:
      values1_raw = self.metrics[self.crit1]
      for i, (value, epoch) in enumerate(zip(values1_raw, self.epochs)):
        self.metrics.iloc[i, self.index1-1] = (int(epoch) * self.ad_steps + value)

    if (self.steps == self.crit2) & self.adaptive:
      values1_raw = self.metrics[self.crit2]
      for i, (value, epoch) in enumerate(zip(values1_raw, self.epochs)):
        self.metrics.iloc[i, self.index2-1] = (int(epoch) * self.ad_steps + value)

    values1 = self.metrics[self.crit1].tolist()
    values2 = self.metrics[self.crit2].tolist()

    return values1, values2


def line_select_callback(eclick, erelease):
    x1, y1 = eclick.xdata, eclick.ydata
    x2, y2 = erelease.xdata, erelease.ydata



def toggle_selector(event):
    if event.key in ['Q', 'q'] and toggle_selector.RS.active:
        toggle_selector.RS.set_active(False)
    if event.key in ['A', 'a'] and not toggle_selector.RS.active:
        toggle_selector.RS.set_active(False)


def main(criteria1, criteria2, ad_steps, path=DIR, out_freq=FREQ, output=OUTPUT_FOLDER, numfolders=False, topology=None):
    """

      Description: Rank the traj found in the report files under path
      by the chosen criteria. Finally, output the best n_structs.

      Input:

         Path: Path to look for *report* files in all its subfolders.

         Criteria: Criteria to sort the structures.
         Needs to be the name of one of the Pele's report file column.
         (Default= "Binding Energy")

         n_structs: Numbers of structures to create.

         sort_order: "min" if sorting from lower to higher "max" from high to low.

         out_freq: "Output frequency of our Pele control file"

     Output:

        f_out: Name of the n outpu
    """
    adaptive = is_adaptive()

    reports = find_reports(path, numfolders)


    # Retrieve Column Names
    steps, crit1_name, crit2_name = get_column_names(reports, STEPS, criteria1, criteria2)

    # Data Mining
    min_values = parse_values(reports, criteria1, criteria2, steps, crit1_name, crit2_name)



    # Figure
    fig, current_ax = plt.subplots()

    # Plot data
    data = DataHandler(min_values, crit1_name, crit2_name, criteria1, criteria2, steps, adaptive, ad_steps, current_ax)

    # Plot axis  
    plt.plot(data.values1, data.values2, 'ro')
    plt.title('{} vs {}'.format(crit1_name, crit2_name))
    plt.xlabel(crit1_name)
    plt.ylabel(crit2_name)


    # Plot Callbacks
    cidpress = fig.canvas.mpl_connect('button_press_event', data.on_press)
    cidrealese= fig.canvas.mpl_connect('button_release_event', data.on_release)
    toggle_selector.RS = RectangleSelector(current_ax, line_select_callback,
                         drawtype='box', useblit=True,
                         button=[1, 3],  # don't use middle button
                         minspanx=5, minspany=5,
                         spancoords='pixels',
                         interactive=False)

    # Show Plot on screen
    plt.show()

def find_reports(path, numfolders):
    reports = glob.glob(os.path.join(path, "*/*report*"))
    reports = glob.glob(os.path.join(path, "*report*")) if not reports else reports
    reports = filter_non_numerical_folders(reports, numfolders)
    try:
        reports[0]
    except IndexError:
        raise IndexError("Not report file found. Check you are in adaptive's or Pele root folder")
    return reports

def parse_values(reports, criteria1, criteria2,  steps, crit1_name, crit2_name):
    """

       Description: Parse the 'reports' and create a sorted array
       of size n_structs following the criteria chosen by the user.

    """
    if steps == crit1_name:
      INITIAL_DATA = [(DIR, []),
                      (REPORT, []),
                      (steps, []),
                      (crit2_name, [])
                      ]
    elif steps == crit2_name:
      INITIAL_DATA = [(DIR, []),
                    (REPORT, []),
                    (steps, []),
                    (crit1_name, []),
                    ]
    else:
      INITIAL_DATA = [(DIR, []),
                    (REPORT, []),
                    (steps, []),
                    (crit1_name, []),
                    (crit2_name, [])
                    ]
    min_values = pd.DataFrame.from_items(INITIAL_DATA)
    for file in reports:
        report_number = os.path.basename(file).split("_")[-1]
        try:
            data = pd.read_csv(file, sep='    ', engine='python')
        except pd.errors.EmptyDataError:
            warnings.warn("Report {} corrupted".format(file), UserWarning)
            continue
        if steps == crit1_name:
          selected_data = data.iloc[:, [2, criteria2-1]]
        elif steps == crit1_name:
          selected_data = data.iloc[:, [2, criteria1-1]]
        else:
          selected_data = data.iloc[:, [2, criteria1-1, criteria2-1]]
        selected_data.insert(0, DIR, [file]*selected_data[steps].size)
        selected_data.insert(1, REPORT, [report_number]*selected_data[steps].size)
        min_values = pd.concat([min_values, selected_data])
    return min_values 


def filter_non_numerical_folders(reports, numfolders):
    """
    Filter non numerical folders among
    the folders to parse
    """
    if(numfolders):
        new_reports = [report for report in reports if(os.path.basename(os.path.dirname(report)).isdigit())]
        return new_reports
    else:
        return reports

def get_column_names(reports, steps, criteria1, criteria2):
    data = pd.read_csv(reports[0], sep='    ', engine='python')
    data = list(data)

    return data[int(steps)-1], data[criteria1-1], data[criteria2-1]

def mkdir_p(path):
    try:
        os.makedirs(path)
    except OSError as exc:  # Python >2.5
        if exc.errno == errno.EEXIST and os.path.isdir(path):
            pass
        else:
            raise


if __name__ == "__main__":
    criteria1, criteria2, ad_steps, path, out_freq, output, numfolders, topology = parse_args()
    main(criteria1, criteria2, ad_steps, path, out_freq, output, numfolders, topology)
