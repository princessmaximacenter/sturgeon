# Sturgeon wrapper for intra-operative prediction

This repository contains the Sturgeon CNS classifier and a collection of wrapper scripts to enable live prediction.
During live prediction, the wrapper script watches the input directory for any bam files that are already present and will continue scanning for newly created bam files. 
With modkit, the methylation data is extracted from the bam files and the Sturgeon Classifier is run on the resulting bed file.
```
Per iteration the following are generated:
  - A csv file with confidence score of all subtypes in the Sturgeon model
  - A bar plot with the confidence score of all subtypes in the Sturgeon model
  - A confidence over time plot showing the progression of the confidence scores per sturgeon iteration

Depending on the set frequency, a CNV plot is also created after every n iterations
```


# Sturgeon

Sturgeon is a CNS neural network classifier based on the reference dataset published by [Capper et al., 2018](https://doi.org/10.1038/nature26000).

For more information on the classifier please refer to our [paper](https://www.nature.com/articles/s41586-023-06615-2).

![image](https://github.com/marcpaga/sturgeon/assets/35224321/e512824b-689a-4437-8664-48407db744f3)

If you use Sturgeon, please cite our work and consider starring this repository 🌟

```
Vermeulen, C., Pagès-Gallego, M., Kester, L. et al. Ultra-fast deep-learned CNS tumour classification during surgery. Nature (2023). https://doi.org/10.1038/s41586-023-06615-2
```
## System requirements
Software:
```
Python version =>3.9 <3.10
modkit (or Rust/Cargo for installation): https://github.com/nanoporetech/modkit

```

Hardware:
Testing and validation for the wrapper scripts have been performed on the PromethION 2 Integrated (https://nanoporetech.com/products/sequence/promethion-2)
As long as the output file structure is the same, the scripts should work on any ONT sequencing device. 
```
An Oxford Nanopore Technology Sequencing device. 
```
## Installation

Get the repository.

```
git clone https://github.com/princessmaximacenter/sturgeon
```

### Optional 
To include the models directly to the installation path, download the desired models, see below for available models and links.

And the move the models:
```
cd sturgeon
mv DOWNLOADED_MODEL.zip sturgeon/include/models/DOWNLOADED_MODEL.zip
```

Otherwise, during prediction, you can just pass the path to the zip file.

Install Sturgeon.

```
cd sturgeon # if you haven't

python3 -m venv venv
source venv/bin/activate
python3 -m pip install --upgrade pip
pip3 install -e . --no-cache-dir
```

If you use a pre-compiled binary then there's no necessity for installation.

## Available models

Model files are zip files with the an `onnx` model file and several other files with information about the classes and model calibration. 

The following models are available:

### `General`

- This model contains a similar classification scheme (some classes are merged, see below) as [Capper et al., 2018](https://doi.org/10.1038/nature26000). This means that there are a total of 87 classes (9 control and 78 tumor classes). 
- We recommend to always use this model in parallel, even when using other models.
- Score recommendations:
    
    - `score < 0.8`: inconclusive result. We recommend to wait for additional data; if not possible, consider the top3 highest scoring classes or their class families as whole.
    - `0.8 <= score < 0.95`: confident result that the class is correct. If used during live sequencing, we recommend to wait for additional sequencing data so that the score gets higher, or the score is stable within this range (e.g several predictions like 0.82, 0.87, 0.84, ...).
    - `score >= 0.95`: high confident result that the class is correct.

Download link: https://www.dropbox.com/s/yzca4exl40x9ukw/general.zip?dl=0

Merged classes (these are merged since the subtype differences are based on tumour location or age):

- Embryonal - MB SHH - CHL AD INF merges: "Embryonal - MB SHH - CHL AD" and "Embryonal - MB SHH - INF"

- Ependymal - SUBEPN - ALL merges: "Ependymal - SUBEPN - PF", "Ependymal - SUBEPN - SPINE" and "Ependymal - SUBEPN - ST"

- Other glioma - LGG PA - PA merges: "Other glioma - LGG PA - PA MID" and "Other glioma - LGG PA - PA PF"


### `Brainstem`

- This model constains a reduced number of classes (30) that can only occur in the brainstem region. Included in these classes are also the control classes and a `Other - Non brainstem` class.
- We recommend to always use the general model in parallel with this one.
- This model is only meant to be used when the tumor is located in the brainstem.
- Score recommendations:
    - `score < 0.95`: inconclusive result. We recommend to wait for additional data; if not possible, consider the top3 highest scoring classes or their class families as whole.
    - `score >= 0.95`: high confident result that the class is correct, but should be treated as inconclusive is the predicted class is `Other - Non brainstem`.

Download link: https://www.dropbox.com/s/55hypw7i8tidr0a/brainstem.zip?dl=0

## Quickstart

This program has four main utilities:

- `inputtobed`: convert input files from Guppy (bam) or from Megalodon (txt), to bed files suitable to be used for prediction using this tool.
- `predict`: predict the CNS type from an input file(s) in bed format.
- `live`: watch over an output folder where output files are written live. Then convert and predict them as they come. This is meant to be used during live sequencing.
- `models`: list, add and delete models. Not strictly necessary, as models can be passed to the previous utilities by path. Ignore for binary programs.

Once launched, a log file will be created in the `logs` folder. The location of that folder will be dependent from the current working directory (from where the tool is executed).
Please refer to each utility `--help` for additional info.

### CRITICAL NOTES:
- Guppy can output bam index files (.bai), do NOT use those, they are not compatible. Either provide them via `samtools` or do not provide them. This program comes with pysam and will create the index files itself.
- Prediction relies on proper mapping of the methylation calls. This program expects that the reads have aligned to the [T2T reference genome v2.0](https://s3-us-west-2.amazonaws.com/human-pangenomics/T2T/CHM13/assemblies/analysis_set/chm13v2.0.fa.gz). We also support alignment to `hg38`, but have NOT throughly tested it.
- While this tool can be run with the output from Guppy or Megalodon, we found that the methylation calls from Megalodon work best and we therefore recommend it for optimal results.


## Preparing data into the right format: `inputtobed`

Convert input files (bam or txt) to bed files that can be used as input to predict the CNS type.

It is critical that the data is correctly aligned to a reference genome, otherwise the results will be wrong. We have used the Telomere-to-telomere reference genome (CHM13v2) for all of our work, and it is the default setting for this program. We encourage you to also use it, as we have not thoroughly tested other reference genomes. If you have aligned your data to `hg38`, we encourage you to re-align (or "lift it" with [CrossMap](https://crossmap.sourceforge.net/#)) it to the T2T reference genome. If you still prefer to use `hg38`, then you can pass `--reference-genome hg38` to use the correct probe coordinates. The `hg38` coordinates have been generated via liftover, this lead to the loss of 26 out of 427680 probes.

### Extract modification calls: modkit (recommended)

If you have bam files that contain modifications basecalled with Guppy or Dorado, this is the recommended option. You will have to install [modkit](https://github.com/nanoporetech/modkit), the official ONT tool to extract modifications. Once installed, given a bam file you should do the following:

1. (optional) If, besides 5mC you also have 5hmC modification calls, you should add them up via:
```
modkit adjust-mods --convert h m INPUT.bam OUTPUT.bam
```
The reasoning behind this is that both 5mC and 5hmC poorly react to bisulfite treatment [Huang et al., 2010](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC2811190/), therefore 5hmC sites would appear as 5mC in a methylation microarray experiment. We trained on those kind of experiments, we therefore think that adding the two scores up reflects the most the training data. 

2. You can then extract the 5mC scores using:
```
modkit extract OUTPUT.bam OUTPUT.txt
```

This `OUTPUT.txt` can then be directly used in sturgeon using:
```
sturgeon inputtobed -i MODKIT_OUTPUT_DIR -o OUTPUT_DIR -s modkit
```

Please note that Sturgeon will only use 5mC scores, other modification scores will be filtered out.
Please refer to the [modkit wiki](https://nanoporetech.github.io/modkit/quick_start.html) in case these commands change.

### Alignment bam files: guppy (not recommended)

Please strongly consider using [modkit](https://github.com/nanoporetech/modkit), to extract methylation calls from bam files from guppy we rely on [modbampy](https://github.com/epi2me-labs/modbam2bed/tree/master), which is currently deprecated. Updates to guppy might break modbampy compatibility and could give WRONG results.

Convert a bam that contains methylation calls into the adequate format (.bed) so that it can be used for prediction.

To obtain bam files from Guppy you can add the following arguments:
```
.guppy_basecaller \
--config dna_r9.4.1_e8.1_modbases_5mc_cg_hac.cfg \
--bam_out \
--align_ref chm13v2.0.fa
```

Do NOT pass the `--index` flag to get index files, these index files are not compatible.
Also remember to pass the T2T reference genome v2.0 as reference for alignment.

Convert bam to bed example with demo data:

```
sturgeon inputtobed -i demo/bam -o demo/bam/out -s guppy
```

In `demo/bam/out` you should have some output files, the most important one is `merged_probes_methyl_calls.bed` as it is the input for the prediction tool. This input contains the aggregation of all methylation calls across ALL bam files. If each bam file is a different sample, then process them independently in different folders.

```
sturgeon inputtobed -i demo/bam/example_1.bam -o demo/bam/out_1 -s guppy
sturgeon inputtobed -i demo/bam/example_2.bam -o demo/bam/out_2 -s guppy
sturgeon inputtobed -i demo/bam/example_3.bam -o demo/bam/out_3 -s guppy
```

### Per read methylation txt files: megalodon (not recommended)

Please strongly consider using [modkit](https://github.com/nanoporetech/modkit), to extract methylation calls since megalodon is deprecated by ONT.

Convert a txt file that contains per read methylation calls into the adequate format (.bed) so that it can be used for prediction.

To obtain the txt files you can run megalodon with the following arguments:
```
./megalodon YOUR_FAST5_PATH \
--outputs mods basecalls mappings \
--mappings-format bam \
--reference chm13v2.0.fa \
--write-mods-text \
--mod-motif m CG 0 \
--processes 10 \
--guppy-config res_dna_r941_prom_modbases_5mC_CpG_v001.cfg \
--devices cuda:0 \
--output-directory YOUR_OUTPUT_PATH
```

For further information on how to run megalodon please refer to their [github](https://github.com/nanoporetech/megalodon) page.

Convert txt files from demo data:
```
sturgeon inputtobed -i demo/mega -o demo/mega/out -s megalodon
```

In `demo/mega/out` you should have some output files, the most important one is `merged_probes_methyl_calls.bed` as it is the input for the prediction tool. This input contains the aggregation of all methylation calls across ALL bam files. If each bam file is a different sample, then process them independently in different folders.

## CNS type prediction: `predict`

This mode predicts a set of samples in bed file format given a set of models.
Each bed file is considered to be a single sample and is treated independently.

Example usage with the demo data:
```
sturgeon predict \
-i demo/bed \
-o demo/bed/results/ \
--model-files PATH_TO_MODEL_DIR/general.zip \
--plot-results
```

In `demo/results` there should be a `.csv` file for each sample with the scores for each CNS class. There should also be a `.pdf` with a barplot representing the predicted scores.

Values indicate the score that the model gave to each class. Higher scores indicate higher confidence in the prediction. 

## CNS type prediction while sequencing: `live or post-sequencing`

This program can be used during live basecalling and alignment. It watches over a folder and waits for bam files to be written there. Then it processes them as they come. This program expects that all bam files in that folder come from the same sample, therefore the amount of sequencing for that sample increases over time. In this line, each bam file will not be treated independently, but instead they will be added in a cumulative manner. 
It is assumed that basecalling and alignment are performed directly by the ONT sequencing device. 


Usage:
```commandline
SturgeonLivePrediction
Usage: SturgeonLivePrediction [OPTIONS]

Options:
  -i, --input PATH            Directory of sequencing run for sturgeon
                              analysis
  -o, --output PATH           Directory where results are written.
  -l, --lock PATH             Name of lock file.
  -s, --sturgeon_script PATH  Path to the script that will be called for
                              processing.
  -b, --barcode TEXT          Barcode used in library preparation.
  -f, --freq INTEGER          Number of iterations before merging BAMs and
                              plotting CNV.
  -m, --model PATH            Location of model used for sturgeon prediction
  -u, --utils PATH            Location of utils directory
  -r, --r_script PATH         Location of R script for plotting CNV
  -g, --gridion BOOLEAN       If run is a gridion verification run, some
                              parameters are changed
  -sf, --shutdown_file PATH   Location of shutdown flag
  --gui_activated             Flag to indicate script is run through GUI
  -lr, --live_run             Flag to indicate whether sequencing and sturgeon analysis is live. Default: FALSE
  --help                      Show this message and exit.
```
### **Important notes about usage**
The input directory can be given in the following ways:
```commandline
/location/to/sequencing/run/
/location/to/sequencing/run/bam_pass/
/location/to/sequencing/run/bam_pass/barcodeXX/
```
The code assumes that the files of interest are in a "barcodeXX" directory, so an alias through the use of a samplesheet **cannot** currently be used. 

The output directory cannot exist yet, it will be made by the wrapper script.

A config.yaml file can be found in the python_scripts directory, which is intended for the use with the Docker Container. 

If the --gridion flag is set to True, previous analysis runs that were performed with sturgeon V1.0.0 can be re-analyzed for validation purposes. 
### **The config.yaml found in python_scripts/ assumes the script is run through docker. Change the parameters in config.yaml to fit your situation**

## Example live run
```commandline
SturgeonLivePrediction --input /location/to/sequencing/run/ --output /location/to/output/live_sturgeon_run --barcode 5 \
--freq 10 --live_run
```
## Example post-sequencing run
```commandline
SturgeonLivePrediction --input /location/to/sequencing/run/bam_pass/barcode05 --output /location/to/output/live_sturgeon_run --barcode 5 \
--freq 10
```
The post-sequencing run will automatically shutdown after all bam files in the input directory have been processed. 

# Docker usage
A docker container has also been created for the live prediction. \
This can be installed with:
```commandline
docker pull princessmaximacenter/sturgeon:v2.0.0
```
### **The current docker image does not include the --live_run parameter, it assumes every run is "live" <br>**
### **See ```https://github.com/princessmaximacenter/sturgeon/tree/v2.0.0-docker``` for code-base for docker image**

Usage example: \
The model file still needs to be installed seperately
```commandline
docker run --rm -i -v /location/to/sequencing/run/:/home/docker/sturgeon/input \
-v /location/where/output/dir/is_created/:/home/docker/sturgeon/output \
-v /location/where/modelfile/is.zip:/opt/sturgeon/sturgeon/include/models/general.zip \
-v $PWD:/home/docker/sturgeon \
princessmaximacenter/sturgeon:v2.0.0 SturgeonLivePrediction \
--input /home/docker/sturgeon/input \
--output /home/docker/sturgeon/output/prediction_output/ \
--barcode {barcode}



```