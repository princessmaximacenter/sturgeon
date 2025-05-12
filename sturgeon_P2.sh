#!/bin/bash
INPUT_BAM="$1"
OUTPUT_DIR="$2"
MODEL_FILE="$3"
ITERATION="$4"


SAMPLE_NAME=$(basename "$INPUT_BAM" .bam)
OUTPUT_SUFFIX=("pdf" "png" "csv")

echo "Extract modifications from bam file with modkit..."
echo "modkit adjust-mods --convert h m ${INPUT_BAM} ${OUTPUT_DIR}/modkit/${SAMPLE_NAME}_modkit.bam"
modkit adjust-mods --convert h m "${INPUT_BAM}" "${OUTPUT_DIR}/modkit/${SAMPLE_NAME}_modkit.bam"
echo "modkit extract full ${OUTPUT_DIR}/modkit/${SAMPLE_NAME}_modkit.bam ${OUTPUT_DIR}/modkit/${SAMPLE_NAME}_modkit.txt"
modkit extract full "${OUTPUT_DIR}/modkit/${SAMPLE_NAME}_modkit.bam" "${OUTPUT_DIR}/modkit/${SAMPLE_NAME}_modkit.txt"

echo "Running sturgeon..."
echo "sturgeon inputtobed --margin 50 -i ${OUTPUT_DIR}/modkit/${SAMPLE_NAME}_modkit.txt -o ${OUTPUT_DIR}/modkit -s modkit"
sturgeon inputtobed --margin 50 -i "${OUTPUT_DIR}/modkit/${SAMPLE_NAME}_modkit.txt" -o "${OUTPUT_DIR}/modkit" -s modkit
echo "sturgeon predict -p --i ${OUTPUT_DIR}/modkit/ -o ${OUTPUT_DIR} -m ${MODEL_FILE}"
sturgeon predict -p --i "${OUTPUT_DIR}/modkit/" -o ${OUTPUT_DIR} -m ${MODEL_FILE}

echo "Moving files..."
ITERATION_DIR="${OUTPUT_DIR}/iteration_${ITERATION}"
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