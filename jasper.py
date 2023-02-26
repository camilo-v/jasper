"""
    Main Jasper script. BETA.
    Please do not use in production environments as this script is undergoing active development.
    If you would like to use Jasper in your research, please get in touch.
"""

import os, sys
import shutil
import argparse
import collections
import random
import time
import datetime
import math
import csv
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import shlex
import subprocess as sp

from utilities import  jasper_utilities as jutils


# ------------------------------------------------ Path to R Script ---------------------------------------------------
#
#   Get the path of the R script so that we can call it to draw the Hilbert curve.
#
r_script_path = os.path.join(os.path.dirname(__file__), 'R_scripts')


# ------------------------------------------------------ Main ---------------------------------------------------------
#
#
def main(args):
    """
    Main function of the app.
    Args:
        args: command line arguments.

    Returns:

    """
    jutils.print_jasper_header()

    # 	Pick up the command line arguments
    parser = argparse.ArgumentParser()
    parser.add_argument("--profile", required=True, type=str, help="Path to file with abundance profiles.")
    parser.add_argument("--type", required=True, default="matrix", const="matrix", nargs="?",
                        choices=["matrix", "flint", "kraken"], help="The type of profile (default: %(default)s)")
    parser.add_argument("--order", required=True, default="taxonomic", const="taxonomic", nargs="?",
                        choices=["taxonomic", "labeled"],
                        help="The ordering scheme for the Hilbert curve (default: %(default)s)")
    parser.add_argument("--min_num_levels", type=int, default=100,
                        help="Minimum number of child taxa to draw labels for.")
    parser.add_argument("--output", required=True, type=str,
                        help="Output directory path in which the images will be created.")
    parser.add_argument("--ignore_labels", action="store_true", required=False,
                        help="Ignores labels when '--type' is matrix.")
    parser.add_argument("--verbose", action="store_true", required=False, help="Wordy Terminal output.")

    args = parser.parse_args(args)


    #   Assign the arguments to variables we can use.
    path_to_profile = args.profile
    profile_type    = args.type
    order_scheme    = args.order
    output_path     = args.output
    ignore_labels   = args.ignore_labels
    min_num_levels  = args.min_num_levels
    verbose_output  = args.verbose

    print("[" + time.strftime('%d-%b-%Y %H:%M:%S', time.localtime()) + "]")
    print("[" + time.strftime('%d-%b-%Y %H:%M:%S', time.localtime()) + "] Starting Visualization...")
    print("[" + time.strftime('%d-%b-%Y %H:%M:%S', time.localtime()) + "] Ordering Scheme: " +
          str(order_scheme).upper())

    # -------------------------------------------- Output Directories -------------------------------------------------
    #
    #
    print("[" + time.strftime('%d-%b-%Y %H:%M:%S', time.localtime()) + "] Preparing Output Directory...")

    output_dir = output_path + "/jasper-images"

    if not os.path.exists(output_dir):
        print("[" + time.strftime('%d-%b-%Y %H:%M:%S', time.localtime()) + "]" +
              "  - Output directory does not exist. Creating...")
        os.makedirs(output_dir)

    print("[" + time.strftime('%d-%b-%Y %H:%M:%S', time.localtime()) + "]")


    # -------------------------------------------- Temp Directories ---------------------------------------------------
    #
    #
    temp_dir = output_dir + "/_tmp"
    if not os.path.exists(temp_dir):
        os.makedirs(temp_dir)


    # -----------------------------------------------------------------------------------------------------------------
    # ---------------------------------------------- Profile Type -----------------------------------------------------
    # -----------------------------------------------------------------------------------------------------------------
    #
    #   What type of profiles are we dealing with?
    #
    print("[" + time.strftime('%d-%b-%Y %H:%M:%S', time.localtime()) + "] Analyzing Profiles...")
    print("[" + time.strftime('%d-%b-%Y %H:%M:%S', time.localtime()) + "] Type: " + str(profile_type).upper())

    start_time = time.time()

    #
    #   'matrix' type contains profiles. Each 'row' is a sample, and it contains N-columns, with each column being
    #   a microbial strain, and the last column being a class label for the sample.
    #
    if profile_type == "matrix":

        if ignore_labels:
            print("[" + time.strftime('%d-%b-%Y %H:%M:%S', time.localtime()) + "]  - Ignoring labels in output.")

        print("[" + time.strftime('%d-%b-%Y %H:%M:%S', time.localtime()) + "] Loading...")

        #
        #   Profile's Data Frame. The main data structure in the program.
        #
        profile_df = pd.read_csv(path_to_profile, sep='\t', header=0, index_col=0, engine='python')

        feature_list   = list(profile_df.columns.values)
        sample_list    = list(profile_df.index.values)

        number_of_microbes = len(feature_list) - 1
        number_of_samples  = len(profile_df.index)
        loc_label = len(feature_list) - 1   # The column that contains the Sample's label.

        taxa_name_list = feature_list
        taxa_name_list.pop()

        print("[" + time.strftime('%d-%b-%Y %H:%M:%S', time.localtime()) + "] " + "Number of Microbes: " +
              '{:0,.0f}'.format(number_of_microbes))

        print("[" + time.strftime('%d-%b-%Y %H:%M:%S', time.localtime()) + "] " + "Number of Samples: " +
              '{:0,.0f}'.format(number_of_samples))

        print("[" + time.strftime('%d-%b-%Y %H:%M:%S', time.localtime()) + "]")

        print("[" + time.strftime('%d-%b-%Y %H:%M:%S', time.localtime()) + "] Creating Hilbert Images..." )


        # -------------------------------------------------------------------------------------------------------------
        #
        #   Ordering based on a Labeled Interpretation
        #
        if order_scheme == "labeled":
            print("[" + time.strftime('%d-%b-%Y %H:%M:%S', time.localtime()) + "] Reordering based on Labels...")

            #   Sort the dataframe based on the labels column
            df_sorted_by_labels = profile_df.sort_values("Label")

            #   New dataframe that will only have the average relative abundances of the cohorts
            df_mean_abundance = df_sorted_by_labels[0:0]

            #   Extract a unique list of the labels so that we can extract individual datasets.
            labels_list = list(set(list(df_sorted_by_labels["Label"])))
            labels_list.sort()
            print(labels_list)

            for label in labels_list:
                print("[" + time.strftime('%d-%b-%Y %H:%M:%S', time.localtime()) + "] Extracting " + str(label))

                #
                #   Extract the samples (rows) that have the same "Label" — or extract all the samples from a
                #   labeled condition identify by the label column.
                #
                cohort_df = df_sorted_by_labels[df_sorted_by_labels["Label"] == label]

                print("[" + time.strftime('%d-%b-%Y %H:%M:%S', time.localtime()) + "]   - Calculating mean...")

                cohort_df.loc["mean"] = cohort_df.mean()

                # Set the label of the last cell in the dataframe as the call to 'mean()' above leaves it blank.
                cohort_df.iloc[-1,-1] = label

                df_mean_abundance.loc[label] = cohort_df.loc["mean"]

            #   At this point, the dataframe "df_mean_abundance" should have as many rows as we have cohorts, and each
            #   row should contain the average relative abundance of the microbe in the cohort.

            #   We don't need the last column ('Label') anymore, so we discard it
            df_mean_abundance = df_mean_abundance.iloc[:, :-1]

            #   We'll transpose it so that we can easily sort it.
            df_mean_abundance_t = df_mean_abundance.transpose()
            df_col_names = list(df_mean_abundance_t.columns.values)
            df_mean_abundance_t = df_mean_abundance_t.sort_values(df_col_names, ascending=False)

            #   Boxplots
            boxplot = sns.boxplot(x="Sample Name",
                                  y="value",
                                  data=pd.melt(df_mean_abundance_t),
                                  palette="colorblind")

            boxplot.set_title("Raw Data")

            plot_file_path = temp_dir + "/boxplot-not_normalized.png"
            boxplot.figure.savefig(plot_file_path,
                                   format='png',
                                   dpi=120)


            #   Keep processing the Dataframe with more columns we need.
            df_mean_abundance_t["Max"] = df_mean_abundance_t[df_col_names].max(axis=1)
            df_mean_abundance_t["Max Label"] = df_mean_abundance_t[df_col_names].idxmax(axis=1)

            sample_tmp_file = temp_dir + "/jasper-tmp.txt"

            df_mean_abundance_t.to_csv(sample_tmp_file,
                                       sep="\t",
                                       encoding="utf-8",
                                       header=True,
                                       index=True,
                                       line_terminator="\n")



        # -------------------------------------------------------------------------------------------------------------
        #
        #   Ordering based on taxonomic tree
        #
        if order_scheme == "taxonomic":

            sample_count = 1

            for index, row in profile_df.iterrows():
                sample_id = str(index) # The row name
                sample_label = str(row[loc_label])
                sample_label = sample_label.replace("/","_")
                sample_label = sample_label.replace(" ", "_")

                print("[" + time.strftime('%d-%b-%Y %H:%M:%S', time.localtime()) + "]  - " + str(sample_count) +
                      " Processing " + sample_id + " (" + sample_label + ")", end=" ")

                # Each image will be placed in a subdirectory for the 'label' of the sample.
                sample_label_dir = output_dir + "/" + sample_label
                if not os.path.exists(sample_label_dir):
                    os.makedirs(sample_label_dir)

                abundance_list = list(row)
                abundance_list.pop() # Remove the last element which contains the label

                abundance_tuples = list(zip(taxa_name_list, abundance_list))
                sample_df = pd.DataFrame(abundance_tuples, columns=["Name", "Abundance"])

                sample_tmp_file = temp_dir + "/" + sample_id + ".txt"

                #   Store a row's profile (i.e., a Sample) into a temp file that will be fed to the R script for
                #   visualization.
                sample_df.to_csv(sample_tmp_file,
                                 sep="\t",
                                 encoding="utf-8",
                                 header=False,
                                 index=False,
                                 line_terminator="\n")

                #
                #   R
                #
                #   Prepare the necessary parameters for calling the "hilbert.R" script for drawing a Hilbert image.
                #
                r_work_dir   = temp_dir
                r_output_dir = sample_label_dir
                r_sample_id  = sample_id
                r_col_name   = 1
                r_col_abundance = 2
                r_min_num_labels = min_num_levels

                r_hilbert_cmd_str = r_script_path + '/hilbert.R --work_dir \"' + r_work_dir + '\"' + \
                                                        ' --output_dir \"' + r_output_dir + '\"' +\
                                                        ' --sample_id ' + r_sample_id + \
                                                        ' --col_name ' + str(r_col_name) + \
                                                        ' --col_abundance ' + str(r_col_abundance) + \
                                                        ' --min_num_labels ' + str(r_min_num_labels) + \
                                                        ' --labels'

                r_hilbert_cmd = shlex.split(r_hilbert_cmd_str)

                r_start_time = datetime.datetime.now()

                #   Dispatch an R subprocess that will call the "hilbert.R" script with the correct parameters
                #   for 'this' sample.
                try:
                    r_subprocess = sp.Popen(r_hilbert_cmd, stdin=sp.PIPE, stdout=sp.PIPE, stderr=sp.PIPE)

                    r_output, r_error = r_subprocess.communicate()

                    if verbose_output:
                        print("OUTPUT: " + str(r_output))
                        print("ERROR: " + str(r_error))

                except sp.CalledProcessError as err:
                    print("[Jasper - R ERROR] " + str(err))
                    sys.exit(-1)

                #   We'll print how long it took R to create the image for 'this' sample.
                r_end_time = datetime.datetime.now()
                r_runtime = r_end_time - r_start_time
                print("[" + str(datetime.timedelta(seconds=r_runtime.seconds)) + "]")

                sample_count += 1

                #   Once the subprocess has ended, we can remove the tmp profile file.
                try:
                    os.remove(sample_tmp_file)
                except OSError as e:
                    print("[ERROR] Removing TMP file: " + str(sample_tmp_file) )
                    print(e)

                # break


    #
    #   We are done so we can calculate the overall runtime, display it, and finish.
    #
    end_time = time.time()
    run_time = end_time - start_time

    print("[" + time.strftime('%d-%b-%Y %H:%M:%S', time.localtime()) + "]")
    print("[" + time.strftime('%d-%b-%Y %H:%M:%S', time.localtime()) + "] " + "Analysis Run Time: " +
          str(datetime.timedelta(seconds=run_time)))
    print("[" + time.strftime('%d-%b-%Y %H:%M:%S', time.localtime()) + "] Done.")
    print("[" + time.strftime('%d-%b-%Y %H:%M:%S', time.localtime()) + "]")



# ---------------------------------------------------------------------------------------------------------------------
# ------------------------------------------------------- Init --------------------------------------------------------
# ---------------------------------------------------------------------------------------------------------------------
#
#   App Initializer.
#
if __name__ == "__main__":
    main(sys.argv[1:])


# ---------------------------------------------------------------------------------------------------------------------
# ---------------------------------------------------- End of Line ----------------------------------------------------
# ---------------------------------------------------------------------------------------------------------------------
