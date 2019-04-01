from fpdf import FPDF
import re
import sys
import os
import subprocess
from AdaptivePELE.analysis import plotAdaptive
import argparse

CONTACTS = "contacts"
BE = "bindingEnergy" 
SASA = "sasa"
kind_Print, colModifier, traj_range = "PRINT_BE_RMSD", 4, None


def arg_parse():
    parser = argparse.ArgumentParser(description='Build a summary pdf file of the simulation')
    parser.add_argument('control_file', type=str, help='adaptive control file')
    parser.add_argument('--traj', type=str, help='Trajectory file name i.e. run_traj_', default='trajectory_')
    parser.add_argument('--report', type=str, help='Report file name i.e run_report_', default='report_')
    args = parser.parse_args()

    return args

def retrieve_metrics(control_file):
    metrics_names = ["rmsd", "com_distance", "distanceToPoint", BE, SASA]
    discard = ["random", "constrainAtomToPosition", "sphericalBox", "toThisOtherAtom", "constrainThisAtom", "springConstant", "equilibriumDistance", ]
    try:
        with open(control_file, "r") as f:
            metrics = [line for line in f if "type" in (line.strip()) and not any(x in line for x in discard)]
    except IOError:
        raise IOError("Pele control file {} not found. Check the path inside Adaptive control file".format(control_file))
    pattern = r'"([A-Za-z0-9_\./\\-]*)"'
    metrics = re.findall(pattern, "".join(metrics))
    metrics = [metric for metric in metrics if metric in metrics_names]
    return metrics


def write_report(metrics, resname, initial_column=4, traj="trajectory_", report="report_"):

    OUTPUT = 'adaptive_report.pdf'

    pdf = FPDF()
    pdf.add_page() 

    # Title
    pdf.set_font('Arial', 'B', 15)
    pdf.cell(80)
    pdf.cell(30, 10, 'Simulation Report', align='C')

    # Plot Binding SASA
    plot(1+metrics.index(BE)+initial_column, 1+metrics.index(SASA)+initial_column, ".", "BE.png", "BindingE", "sasa", report=report)
    pdf.cell(-100)
    pdf.set_font('Arial', 'B', 11)
    pdf.cell(10, 49, 'Interaction Energy vs Sasa' +34*"\t" + "Total Energy vs Interaction Energy")
    pdf.image("BE.png", 10, 40, 83)

    # Plot Total Biding
    plot(initial_column,  1+metrics.index(BE)+initial_column, ".", "total.png", "totalE", "BindingE", report=report)
    pdf.cell(0)
    pdf.set_font('Arial', 'B', 11)
    pdf.image("total.png", 100, 40, 83)

    # Contacts
    create_contact_plot(".")
    pdf.cell(0)
    pdf.set_font('Arial', 'B', 11)
    pdf.image("{}_threshold.png".format(CONTACTS), 10, 110, 83)
    pdf.image("{}_hist.png".format(CONTACTS), 100, 110, 83)
    

    pdf.add_page()

    #Plot other metrics agains binding
    images_per_page = 1
    images_per_row = 0
    #metrics = ["rmsd", "com_distance", "distanceToPoint", "clusters"]
    metrics_names = ["rmsd", "com_distance", "distanceToPoint"]
    for user_metric in metrics_names:
        if images_per_page > 4:
            pdf.add_page()
            images_per_page = 1
        if user_metric in metrics or user_metric == "clusters":
            pdf, images_per_page, images_per_row = write_metric(
	        pdf, user_metric, BE, metrics, images_per_page, images_per_row, resname, traj=traj, report=report)

    #Output report    
    pdf.output(OUTPUT, 'F')


def plot_clusters(metric1, metric2, metrics, resname, traj, report):
    command = "python -m AdaptivePELE.analysis.clusterAdaptiveRun 200 {} {} {} --report {} --traj {} --png".format(
       get_column(metric1, metrics), get_column(metric2, metrics), resname, report, traj)
    os.system(command)

def get_column(metric, metrics, initial_pos=4):
    return 1+metrics.index(metric)+initial_pos

def write_metric(pdf, metric1, metric2, metrics, images_per_page, images_per_row, resname=None, initial_pos=4, traj="trajectory_", report="report_"):

    #Create Image
    if metric1 == "clusters":
        image_name = "Cluster_analisis/ClusterMap.png".format(metric1)
        plot_clusters("bindingEnergy", "sasa", metrics, resname, report)
    else:
        image_name = "{}.png".format(metric1)
        plot(1+metrics.index(metric1)+initial_pos,  1+metrics.index(metric2)+initial_pos, ".",
          image_name, metric1, metric2, zcol=initial_pos + 1+ metrics.index("sasa"), report=report)

    #Move over pdf
    pdf.image(image_name, pdf.get_x(), pdf.get_y(), 83)
    images_per_row += 1
    if images_per_row == 1:
        pdf.cell(80, 0, "", 0, 0)
    else:
        pdf.set_xy(10, 100)
        pdf.cell(80, 0, "", 0, 2);
        images_per_row = 0
    images_per_page += 1
    
    #Update cursor
    return pdf, images_per_page, images_per_row
    
def plot(Xcol, Ycol, path, name, xcol_name, ycol_name, zcol=5, report="report_"):
    try:
        plot_line = plotAdaptive.generatePrintString(8, Xcol, Ycol, report, "PRINT_BE_RMSD", zcol, None).strip("\n")
    except TypeError:
        raise TypeError("Report not found use the flag --report. i.e --report run_report_ ") 
    command = '''gnuplot -e "set terminal png; set output '{}'; set xlabel '{}'; set ylabel '{}'; {}"'''.format(
        name, xcol_name, ycol_name, plot_line)
    os.system(command)
        
def create_contact_plot(path, filenames=CONTACTS):
    command = "python -m AdaptivePELE.analysis.numberOfClusters -f contacts"
    os.system(command)
   
def retrieve_fields(control_file):
    with open(control_file, "r") as f:
        content = f.readlines()
        control_file_line = [line for line in content if line.strip().startswith('"controlFile"')][0]
        ligand_res_line = [line for line in content if line.strip().startswith('"ligandResname"')][0]
    pele, resname = control_file_line.split(":")[1].strip().strip(",").strip('"'), ligand_res_line.split(":")[1].strip().strip(",").strip('"')
    if not os.path.isfile(pele):
        path = os.path.dirname(os.path.abspath(control_file))
        pele = os.path.join(path, pele)
    return pele, resname



def main(control_file, traj, report):
    print("Search pele control file")
    pele_conf, resname = retrieve_fields(control_file)
    print("Retrieve metrics")
    metrics = retrieve_metrics(pele_conf)
    print("Build report")
    write_report(metrics, resname, 4, traj, report)
    print("Analysis finished succesfully")
 

if __name__ == "__main__":
    args = arg_parse()
    main(args.control_file, args.traj, args.report)
