import matplotlib
matplotlib.use('Agg')  # Non-interactive backend (for saving to file)

import seaborn as sns
import matplotlib.pyplot as plt
from matplotlib.ticker import ScalarFormatter
import numpy as np
import pandas as pd
from pathlib import Path
from datetime import datetime
from rpy2.robjects import r, globalenv
from rpy2.robjects.conversion import localconverter
from rpy2.robjects import default_converter



def write_progress_tsv(full_data: pd.DataFrame,output_folder: Path,iteration: int,modelname: str) -> pd.DataFrame:
    """
    Writes the classified scores of the most recent iteration to a tsv and updates the dataframe, used for plotting the confidence over time
    :param full_data: Dataframe  containing all classifier confidence scores of all iterations
    :param output_folder: Path to output directory
    :param iteration: Current iteration
    :param modelname: Name of model used for classification
    :return: Updated DataFrame with most recent iteration included
    """
    current_csv = Path(f"{output_folder}/iteration_{iteration}/merged_probes_methyl_calls_{modelname}_iteration_{iteration}.csv")
    current_results = pd.read_csv(current_csv,sep=',',header=None).T.iloc[1:] #Read csv, transpose, and drop number_probes row
    current_results.columns = ["class", "score"]
    #Add the time stamp
    now = datetime.now()
    TimeStamp = now.strftime("%H:%M")
    current_results.loc[len(current_results)] = ["TIME", TimeStamp]
    if (full_data.empty):
        current_results.rename(columns={"score": f"iteration_{iteration}"}, inplace=True)
        return current_results
    else:
        mgd = pd.merge(full_data, current_results, on="class")
        mgd.rename(columns={"score": f"iteration_{iteration}"}, inplace=True)
        return mgd

def plot_confidence_over_time(full_data: pd.DataFrame,output_file: str,color_translation: dict) -> None:
    """
    Creates confidence over time plots, in png and pdf format
    :param full_data: Dataframe containing all the classifier confidence scores of all iterations
    :param output_file: Path to output file
    :param color_translation: Dict containing per label a color hex code
    :return: None
    """
    #Filter out time from dataframe
    time_row = full_data[full_data["class"] == "TIME"]
    full_data = full_data[full_data["class"] != "TIME"]
    full_data.set_index("class",inplace=True)
    full_data = full_data.astype(float)


    #Initiate plot
    plt.figure(figsize=(12,6))
    sns.set_theme(style="whitegrid")
    #Plot each class
    for class_name, row in full_data.iterrows():

        y_values = row.values
        x_values = full_data.columns
        if "/" in class_name:
            classNameSubsitute = class_name.replace("/"," ")
            color = color_translation.get(classNameSubsitute,"gray")
        else:
            color = color_translation.get(class_name,"gray")



        plt.plot(x_values,y_values,label=class_name, color=color, marker="o",linestyle = "-")

        #Label if final confidence > 0.5
        last_confidence = float(y_values[-1])
        if last_confidence > 0.5:
            plt.text(len(y_values) - 1, last_confidence+0.05, class_name, fontsize=9)

    #Add threshold lines
    plt.axhline(0.95,color='red',linestyle='--', label="0.95 threshold")
    plt.axhline(0.80, color='orange', linestyle='--', label="0.80 threshold")

    #Set axis labels and title
    plt.xlabel("Iteration")
    plt.ylabel("Confidence")
    plt.title("Confidence over time")

    num_iterations = full_data.shape[1]
    x_upper = num_iterations + 5 if num_iterations > 5 else 10
    """Only show time labels every 4 iterations, for readability"""
    tick_labels = time_row.iloc[0,1:].tolist()
    padded_time_labels = tick_labels + [""] * (x_upper - len(tick_labels))
    tick_indices = list(range(x_upper))
    iteration_labels = [str(i + 1) for i in tick_indices]
    sparse_time_labels = [padded_time_labels[i] if (i + 1) % 4 == 0 else "" for i in tick_indices]
    combined_labels = [f"{iter_label}\n{time}" for iter_label, time in zip(iteration_labels, sparse_time_labels)]

    plt.xticks(ticks=tick_indices, labels=combined_labels, rotation=0)
    plt.xlim(-0.8,x_upper-0.5)
    plt.yticks(np.arange(0, 1.1, step=0.2))
    plt.ylim(-0.1,1.1)
    plt.gca().yaxis.set_major_formatter(ScalarFormatter())
    plt.tight_layout()
    # Save to files
    plt.savefig(f"{output_file}.png")
    plt.savefig(f"{output_file}.pdf")
    plt.close()

def plot_CNV_bam(input_bam: Path, output_file: Path, r_script: Path, utils: Path) -> None:
    """
    Runs the R script to create the CNV plot
    :param input_bam: Merged bam file used for generating CNV data
    :param output_file: path to output file
    :param r_script: path to r script used for plotting
    :return: None
    """
    with localconverter(default_converter) as cv:
        r(f'source("{r_script}")')
        r_plot_cnv = globalenv['plot_cnv_from_bam_DNAcopy']
        r_plot_cnv(input_bam, output_file, utils)


