#!/bin/bash
<<com
This script runs the sturgeon classifier on the bam files output from the P2 basecalling and alignment.
First modkit is used to extract the methylation signal from the bam file
Secondly, sturgeon inputtobed is run to create the input file for the sturgeon predictor
Lastly, sturgeon predict is run for prediction, using model version 1 or 2.
com
INPUT_BAM="$1"
OUTPUT_DIR="$2"
MODEL_FILE="$3"
ITERATION="$4"
Sturgeon_V2=${5:-'false'}
CONF_PLOT="$6"

SAMPLE_NAME=$(basename "$INPUT_BAM" .bam)
OUTPUT_SUFFIX=("pdf" "png" "csv")
ITERATION_DIR="${OUTPUT_DIR}/iteration_${ITERATION}"

function extract_methylation_calls() {
  echo "Extract modifications from bam file with modkit..."
  echo "modkit adjust-mods --convert h m ${INPUT_BAM} ${OUTPUT_DIR}/modkit/${SAMPLE_NAME}_modkit.bam"
  modkit adjust-mods --convert h m "${INPUT_BAM}" "${OUTPUT_DIR}/modkit/${SAMPLE_NAME}_modkit.bam"
  echo "modkit extract full ${OUTPUT_DIR}/modkit/${SAMPLE_NAME}_modkit.bam ${OUTPUT_DIR}/modkit/${SAMPLE_NAME}_modkit.txt"
  modkit extract full "${OUTPUT_DIR}/modkit/${SAMPLE_NAME}_modkit.bam" "${OUTPUT_DIR}/modkit/${SAMPLE_NAME}_modkit_${ITERATION}.txt"

  echo "sturgeon inputtobed --margin 50 -i ${OUTPUT_DIR}/modkit/ -o ${OUTPUT_DIR}/modkit -s modkit"
  sturgeon inputtobed --margin 50 -i "${OUTPUT_DIR}/modkit/" -o "${OUTPUT_DIR}/modkit" -s modkit
}

function sturgeon_prediction_v1() {
  echo "Sturgeon Model Version 1 will be used for the prediction"
  echo "sturgeon predict -p --i ${OUTPUT_DIR}/modkit/ -o ${OUTPUT_DIR} -m ${MODEL_FILE}"
  sturgeon predict -p --i "${OUTPUT_DIR}/modkit/" -o ${OUTPUT_DIR} -m ${MODEL_FILE}
}

function sturgeon_prediction_v2() {
  echo "Sturgeon Model Version 2 will be used for the prediction"
  echo "sturgeon-v2 -i "${OUTPUT_DIR}/modkit/merged_probes_methyl_calls.bed" --model ${MODEL_FILE} -o "${OUTPUT_DIR}/merged_probes_methyl_calls_v2.csv" -f bed"
  sturgeon-v2 -i "${OUTPUT_DIR}/modkit/merged_probes_methyl_calls.bed" --model ${MODEL_FILE} -o "${OUTPUT_DIR}/merged_probes_methyl_calls_cns-v2.csv" -f bed
  # Create the prediction over time plot
  python ${CONF_PLOT} -i "${OUTPUT_DIR}/merged_probes_methyl_calls_cns-v2.csv" -m ${MODEL_FILE} -o ${OUTPUT_DIR}


  }

function move_results() {
  echo "Moving files..."

  mkdir ${ITERATION_DIR}
  for EXT in "${OUTPUT_SUFFIX[@]}"; do
    for FILE in "${OUTPUT_DIR}"/*."${EXT}"; do
      # Skip already renamed results
      if [[ "$(basename "$FILE")" == *iteration* ]]; then
        continue
      fi
      FILENAME=$(basename "$FILE")
      mv "$FILE" "${ITERATION_DIR}/${FILENAME%.*}_iteration_${ITERATION}.${FILENAME##*.}"
    done
  done

  echo "Results moved to ${ITERATION_DIR}"
}


echo "FLAG: starting sturgeon for iteration_${ITERATION}"

extract_methylation_calls

if [[ "$Sturgeon_V2" = "true" || "$Sturgeon_V2" = "True" ]] ; then
  sturgeon_prediction_v2
else
  sturgeon_prediction_v1

fi

move_results

