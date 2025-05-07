from unittest.mock import inplace
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
import rpy2.robjects as robjects


def write_progress_tsv(full_data,output_folder,iteration,modelname):
    current_csv = Path(f"{output_folder}/iteration_{iteration}/merged_probes_methyl_calls_{modelname}_iteration_{iteration}.csv")
    current_results = pd.read_csv(current_csv,sep=',',header=None).T.iloc[1:] #Read csv, transpose, and drop number_probes row
    current_results.columns = ["class", "score"]
    #Add the color configuration
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

def plot_confidence_over_time(full_data,output_file,color_translation):
    """
    Write plot confidence over time
    y-axis: confidence
    x-axis: iteration
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
        color = color_translation.get(class_name,"gray")



        plt.plot(x_values,y_values,label=class_name, color=color, marker="o",linestyle = "")

        #Label if final confidence > 0.5
        for i,value in enumerate(y_values):
            confidence = float(value)
            if confidence > 0.5:
                plt.text(i, confidence+0.05, class_name, fontsize=9)

    #Add threshold lines
    plt.axhline(0.95,color='red',linestyle='--', label="0.95 threshold")
    plt.axhline(0.80, color='orange', linestyle='--', label="0.80 threshold")

    #Set axis labels and title
    plt.xlabel("Iteration")
    plt.ylabel("Confidence")
    plt.title("Confidence over time")

    num_iterations = full_data.shape[1]
    x_upper = num_iterations + 5 if num_iterations > 5 else 10
    tick_labels = time_row.iloc[0,1:].tolist()
    iteration_labels = [str(i + 1) for i in range(x_upper)]


    padded_time_labels = tick_labels + [""] * (x_upper - len(tick_labels))
    combined_labels = [f"{iter_label}\n{time}" for iter_label, time in zip(iteration_labels, padded_time_labels)]
    plt.xticks(ticks=range(x_upper), labels=combined_labels, rotation=0)
    plt.xlim(-0.8,x_upper-0.5)
    plt.yticks(np.arange(0, 1.1, step=0.2))
    plt.ylim(-0.1,1.1)
    plt.gca().yaxis.set_major_formatter(ScalarFormatter())
    plt.tight_layout()
    # Save to files
    plt.savefig(f"{output_file}.png")
    plt.savefig(f"{output_file}.pdf")
    plt.close()

def plot_CNV_bam(input_bam, output_file, r_script):
    with localconverter(default_converter) as cv:
        r(f'source("{r_script}")')
        r_plot_cnv = globalenv['plot_cnv_from_bam_DNAcopy']
        r_plot_cnv(input_bam, output_file)


