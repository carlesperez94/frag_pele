import argparse
import AdaptivePELE.analysis.bestStructs as bs
import AdaptivePELE.analysis.backtrackAdaptiveTrajectory as bt

HELP = " \
Command:\n \
----------\n \
python -m AdaptivePELE.analysis.unbinding column_name_of_report --epoch epoch_folder --top topology_file.pdb --n number_of_structures_to_output \n \
i.e: python -m AdaptivePELE.analysis.unbinding sasaLig --epoch 11 --top topology.pdb \
"


def parse_args():
    parser = argparse.ArgumentParser(description=HELP, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('crit', type=str, help=' Column of Pele report to build the exit path with (Normally sasaLig)')
    parser.add_argument('--epoch', type=str, help='Epoch where to finish the exit path', default=".")
    parser.add_argument('--n', type=int, help='Topology file. It must be specified for .xtc files not for pdb', default="1")
    parser.add_argument('--top', type=str, help='Number of exit path to retrieve', default=None)
    parser.add_argument('--outputPath', type=str, help='Output path. i.e: exit_path', default=".")
    parser.add_argument('--outputFile', type=str, help='Output file name. i.e: path.pdb', default=None)
    parser.add_argument('--sort', type=str, help='The end of the exit part will be the max or min value of the metric. default="max". i.e --sort min', default="max")
    args = parser.parse_args()

    return args.crit, args.epoch, args.n, args.top, args.outputPath, args.outputFile, args.sort


def analyze_unbinding(crit1, epoch=".", top=None, n=1, outputPath=".", out_filename=None, sort="max"):
    files_out, epochs, file_ids, step_indexes = bs.main(crit1, path=epoch, topology=top, n_structs=n, output=outputPath, sort_order=sort)
    initial_filename = out_filename
    for epoch, traj, snap in zip(epochs, file_ids, step_indexes):
        if not initial_filename:
            out_filename = "{}_{}_{}_pathway.pdb".format(epoch, traj, snap)
        bt.main(int(traj), int(snap), epoch, outputPath, out_filename, top)


if __name__ == "__main__":
    crit, epoch, n, top, output_path, output_name, sort = parse_args()
    analyze_unbinding(crit, epoch=epoch, top=top, n=n, outputPath=output_path, out_filename=output_name, sort=sort)
