from pathlib import Path
import yaml
import zipfile
import numpy as np
import pandas as pd
import click
from matplotlib import pyplot as plt

from live_prediction_wrapper import SturgeonLogging as SL
app_log = SL._get_app_logger()

def bar_plot(input, model,output):
    with zipfile.ZipFile(model, 'r') as zipf:
        app_log.info("Loading colors dict")
        try:
            color_dict = yaml.safe_load(zipf.open('cns-v2/classification_system.yaml'))
        except FileNotFoundError:
            app_log.info("No colors dict found in zip file")
            color_dict = None

    prediction_df = pd.read_csv(input)
    output_prefix = input.stem


    SMALL_SIZE = 8
    MEDIUM_SIZE = 12
    BIGGER_SIZE = 18

    plt.rc('font', size=SMALL_SIZE)  # controls default text sizes
    plt.rc('axes', titlesize=MEDIUM_SIZE)  # fontsize of the axes title
    plt.rc('axes', labelsize=MEDIUM_SIZE)  # fontsize of the x and y labels
    plt.rc('xtick', labelsize=MEDIUM_SIZE)  # fontsize of the tick labels
    plt.rc('ytick', labelsize=MEDIUM_SIZE)  # fontsize of the tick labels
    plt.rc('legend', fontsize=SMALL_SIZE)  # legend fontsize
    plt.rc('figure', titlesize=BIGGER_SIZE)  # fontsize of the figure title

    non_label_columns = ['probes', 'timestamp']
    label_columns = []
    for c in prediction_df.columns:
        if c not in non_label_columns:
            label_columns.append(c)

    title = Path(output_prefix).stem + '\n' + 'Measured probes: {}'.format(prediction_df['probes'].item())

    plt.figure(figsize=(len(label_columns) // 4, 7))
    for i, c in enumerate(label_columns):
        v = prediction_df[c].item()

        if color_dict is not None:
            color = color_dict['colors']['type'][c]
        else:
            color = 'grey'

        plt.axhline(0.8, color='grey', linestyle='--', zorder=1)
        plt.axhline(0.95, color='grey', linestyle='--', zorder=1)
        plt.bar(
            x=i,
            height=v,
            edgecolor='black',
            color=color,
            zorder=2
        )

    plt.xlim(-1, len(label_columns))
    plt.ylim(0, 1.05)
    plt.xticks(np.arange(len(label_columns)), label_columns, rotation=90)
    plt.ylabel('Score')
    plt.title(title)
    plt.tight_layout()
    plt.savefig(f"{output}/{output_prefix}.pdf", bbox_inches="tight", dpi=300)
    plt.savefig(f"{output}/{output_prefix}.png", bbox_inches="tight", dpi=300)
    plt.close()

def click_command(func):
    @click.command()
    @click.option(
        "-i", "--input", type=click.Path(path_type=Path, exists=False, dir_okay=True), default=None, help="CSV file with prediction scores or current iteration"
    )
    @click.option(
        "-m", "--model", type=str, default=None, help="Location of model used for sturgeon prediction"
    )
    @click.option(
        "-o", "--output",type=click.Path(path_type=Path, exists=False, dir_okay=True), default=None, help="Directory where plots will be saved"
    )
    def wrapper(*args, **kwargs):
        return func(*args, **kwargs)

    return wrapper

@click_command
def main(input: Path, model: Path, output:Path) -> None:
    bar_plot(input,model,output)

if __name__ == "__main__":
    main()