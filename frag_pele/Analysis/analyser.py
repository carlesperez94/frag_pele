import sys
import os
import re
import glob
import pandas as pd
import argparse


def parse_arguments():
    """
            Parse user arguments

            Output: list with all the user arguments
        """
    # All the docstrings are very provisional and some of them are old, they would be changed in further steps!!
    parser = argparse.ArgumentParser(description="""Computes the mean of the 25% lowest values of the sampling simulation
    for each fragment grown. """)
    required_named = parser.add_argument_group('required named arguments')
    # Growing related arguments
    required_named.add_argument("path_to_analyze",
                                help="""Path of the folder to be analyzed, where all simulations are stored.""")

    parser.add_argument("-r", "--rep_pref", default='report_',
                                help="""Prefix of the report file: usually 'report_'. """)
    parser.add_argument("-s", "--steps", default=False,
                        help="Used to filter by N PELE steps")
    parser.add_argument("-f", "--file", default=False,
                        help="Used to set the flag 'export file? to True'")
    parser.add_argument("-o", "--out", default="",
                        help="Name of the output file")
    parser.add_argument("-e", "--equil_folder", default="sampling_results*",
                        help="Prefix of the equilibration folder (where the results are)")
    parser.add_argument("-c", "--col", default="Binding Energy",
                        help="Name of the column that we will use to find the minimum value of the PELE report")
    parser.add_argument("-q", "--quart", default=0.25,
                        help="Quartile that will be used to compute the mean")
    parser.add_argument("-csv", "--export_csv", default=True,
                        help="By default a csv file with a summary of all reports will be stored. If you do not"
                        "want to save it, set the flag to False.")
    parser.add_argument("-l", "--limit_column", default=None, type=str,
                        help="Column used to subset the results.")
    parser.add_argument("-lu", "--limit_up", default=None, type=float,
                        help="Top value to create the new subset of samples.")
    parser.add_argument("-ld", "--limit_down", default=None, type=float,
                        help="Lowest value to create the new subset of samples.")

    args = parser.parse_args()

    return args.path_to_analyze, args.rep_pref, args.steps, args.file, args.out, args.equil_folder, args.col, \
           args.quart, args.export_csv, args.limit_column, args.limit_up, args.limit_down


def pele_report2pandas(path, export=True):
    """
        This function merge the content of different report for PELE simulations in a single file pandas Data Frame.
    """
    data = []
    report_list = glob.glob('{}*'.format(path))
    for report in report_list:
        tmp_data = pd.read_csv(report, sep='    ', engine='python')
        tmp_data = tmp_data.iloc[1:]  # We must discard the first row
        processor = re.findall('\d+$'.format(path), report)
        tmp_data['Processor'] = processor[0]
        data.append(tmp_data)
    result = pd.concat(data)
    if export:
        path_to_folder = "/".join(path.split("/")[0:-1])
        new_dir = path_to_folder + "/summary"
        if not os.path.exists(new_dir):
            os.mkdir(new_dir)
        while True:
            counter = 1
            try:
                result.to_csv("{}/summary_report_{}.csv".format(new_dir, counter))
                print("Data saved in {}/summary_report_{}.csv".format(new_dir, counter))
                break
            except FileExistsError:
                counter =+ 1
                result.to_csv("{}/summary_report_{}.csv".format(new_dir, counter))

    return result


def select_subset_by_steps(dataframe, steps):
    """
    Given a pandas dataframe from PELE's result it returns a subset of n PELE steps.
    :param dataframe: pandas object
    :param steps: int
    :return: dataframe
    """
    subset = dataframe[dataframe['Step'] <= int(steps)]
    return subset


def get_min_value(dataframe, column):
    minimum = dataframe.loc[dataframe[column] == dataframe[column].min()]
    if len(minimum) > 1:
        return minimum.iloc[0][column].min()
    else:
        return minimum[column].min()


def compute_mean_quantile(dataframe, column, quantile_value, limit_col=None, limit_up=None, limit_down=None):
    if limit_col:
        if not limit_up:
            raise ValueError("You must fill the argument '--limit_up'!")
        if not limit_down:
            raise ValueError("You must fill the argument '--limit_down'!")
        dataframe = dataframe[dataframe[limit_col] > float(limit_down)]
        dataframe = dataframe[dataframe[limit_col] < float(limit_up)]
    dataframe = dataframe[dataframe[column] < dataframe[column].quantile(quantile_value)]
    dataframe = dataframe[column]
    mean_subset = dataframe.mean()
    return mean_subset


def get_csv(dataframe, output_path):
    dataframe.to_csv(path_or_buf=output_path, sep='\t', header=True, index=False)


def compute_sterr(dataframe, column, quantile_value):
    dataframe = dataframe[dataframe[column] < dataframe[column].quantile(quantile_value)]
    dataframe = dataframe[column]
    err_subset = dataframe.sem()
    return err_subset


def get_score_for_folder(report_prefix, path_to_equilibration, steps=False,
                         column="Binding Energy", quantile_value=0.25, export=True,
                         limit_col=None, limit_up=None, limit_down=None):
    df = pele_report2pandas(os.path.join(path_to_equilibration, report_prefix), export)
    if steps:
        df = select_subset_by_steps(df, steps)
    mean_quartile = compute_mean_quantile(df, column, quantile_value, limit_col, limit_up, limit_down)
    results = [path_to_equilibration, mean_quartile]
    print("SCORING results: {}	{}".format(results[0], results[1]))
    return results


def analyse_at_epoch(report_prefix, path_to_equilibration, execution_dir, steps=False, column="Binding Energy", quantile_value=0.25):
    result = get_score_for_folder(report_prefix=report_prefix, path_to_equilibration=path_to_equilibration,
                                  steps=steps, column=column, quantile_value=quantile_value)
    out_file = os.path.join(execution_dir, "simulation_score_summary.tsv")
    if os.path.exists(out_file):
        #Pandas new version compatibility
        try:
            df = pd.read_csv(out_file, sep="\t", header=0, ignore_index=True)
        except TypeError:
            df = pd.read_csv(out_file, sep="\t", header=0, index_col=False)
        df = pd.concat([df, pd.DataFrame([result],  columns=["Fragment_Results_Folder", "Score"])])
    else:
        df = pd.DataFrame([result], columns=["Fragment_Results_Folder", "Score"])
    df.to_csv(out_file, sep="\t", index=False)


def main(report_prefix, path_to_equilibration, equil_pattern="equilibration*", steps=False, out_report=False,
         column="Binding Energy", quantile_value=0.25, export=True, limit_col=None, limit_up=None, limit_down=None):
    folder_list = glob.glob(os.path.join(path_to_equilibration, equil_pattern))
    for folder in list(folder_list):
        try:
            result = get_score_for_folder(report_prefix=report_prefix, path_to_equilibration=folder, steps=steps,
                                          column=column, quantile_value=quantile_value, export=export,
                                          limit_col=limit_col, limit_up=limit_up, limit_down=limit_down)
            line_of_report = "{}\t{:.2f}\n".format(result[0], float(result[1]))
            print(line_of_report)
            if out_report:
                if os.path.exists(out_report):
                    with open(out_report, "a") as report:
                        report.write(line_of_report)
                else:
                    with open(out_report, "w") as report:
                        report.write("# FragmentResultsFolder\tFrAG-score\n")
        except Exception as e:
            print("ERROR {} in {}".format(e, folder))


if __name__ == '__main__':
    path, rep_pref, steps, file, out, equil_folder, col, quart, export, limit_col, limit_up, limit_down= \
        parse_arguments()
    main(report_prefix=rep_pref, path_to_equilibration=path, equil_pattern=equil_folder, steps=steps, out_report=out,
         column=col, quantile_value=quart, export=export, limit_col=limit_col, limit_up=limit_up, limit_down=limit_down)



