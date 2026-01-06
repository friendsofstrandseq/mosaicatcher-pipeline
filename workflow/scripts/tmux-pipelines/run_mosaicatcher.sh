#!/bin/bash

# Check if all arguments are provided
if [ "$#" -ne 2 ]; then
    echo "Usage: $0 <ARG1> <ARG2>"
    exit 1
fi

ARG1=$1
ARG2=$2

# Set the session name
SESSION_NAME=${ARG1}

# Start tmux session if it doesn't already exist
tmux has-session -t $SESSION_NAME 2>/dev/null
if [ $? != 0 ]; then
    tmux new-session -d -s $SESSION_NAME 
    
    # Separate samples into two batches
    FIRST_BATCH=()
    SECOND_BATCH=()

    for SAMPLE in ${ARG2}/*; do
        SAMPLE_NAME=$(basename "$SAMPLE")
        if [[ $SAMPLE_NAME == GM* || $SAMPLE_NAME == HG* || $SAMPLE_NAME == *pool* || $SAMPLE_NAME == *POOL* ]]; then
            FIRST_BATCH+=("$SAMPLE_NAME")
        else
            SECOND_BATCH+=("$SAMPLE_NAME")
        fi
    done

    # If no first batch, use a single batch
    if [ ${#FIRST_BATCH[@]} -eq 0 ]; then
        SECOND_BATCH=($(basename -a ${ARG2}/*))
    fi

    # Create temporary directories for each batch
    TMP_DIR_BASE="/g/korbel/STOCKS_WF/mosaicatcher-pipeline-TMP/${ARG1}"
    FIRST_BATCH_DIR="${TMP_DIR_BASE}/FIRST_BATCH"
    SECOND_BATCH_DIR="${TMP_DIR_BASE}/SECOND_BATCH"

    mkdir -p "$FIRST_BATCH_DIR"
    mkdir -p "$SECOND_BATCH_DIR"

    # Symlink samples into their respective temporary directories
    for SAMPLE in "${FIRST_BATCH[@]}"; do
        ln -fs "${ARG2}/$SAMPLE" "$FIRST_BATCH_DIR/$SAMPLE"
    done

    for SAMPLE in "${SECOND_BATCH[@]}"; do
        ln -fs "${ARG2}/$SAMPLE" "$SECOND_BATCH_DIR/$SAMPLE"
    done

    # Create a string for samples_to_process
    FIRST_BATCH_STRING="\"[$(IFS=,; echo "${FIRST_BATCH[*]}" | sed 's/ /,/g')]\""
    SECOND_BATCH_STRING="\"[$(IFS=,; echo "${SECOND_BATCH[*]}" | sed 's/ /,/g')]\""
    echo "First batch string: ${FIRST_BATCH_STRING}"
    echo "Second batch string: ${SECOND_BATCH_STRING}"

    

    # Create tmux window for the first batch if it exists
    if [ ${#FIRST_BATCH[@]} -gt 0 ]; then
        tmux new-window -t $SESSION_NAME -n "first_batch"
        tmux send-keys -t $SESSION_NAME:first_batch "zsh" C-m
        tmux send-keys -t $SESSION_NAME:first_batch "cd /g/korbel2/weber/workspace/StrandSeq_workspace/DEV/mosaicatcher-pipeline-friendsofstrandseq" C-m
        tmux send-keys -t $SESSION_NAME:first_batch "conda activate snakemake-v7-19-1" C-m
        tmux send-keys -t $SESSION_NAME:first_batch "snakemake \
            --config \
                data_location=${FIRST_BATCH_DIR} \
                samples_to_process=${FIRST_BATCH_STRING} \
                multistep_normalisation=False \
                MultiQC=False \
                genome_browsing_files_generation=False \
                ashleys_pipeline=True \
                bypass_ashleys=False \
            --rerun-triggers mtime \
            --rerun-incomplete \
            --profile workflow/snakemake_profiles/mosaicatcher-pipeline/v7/HPC/slurm_EMBL/ \
            --set-resources \
                ashleys_mark_duplicates:constraint="rome" \
                ashleys_generate_features:time="10:00:00" \
                ashleys_bwa_strandseq_to_reference_alignment:time="24:00:00" \
            --jobs 500 \
            -k" C-m    
    fi

    # Create tmux window for the second batch
    if [ ${#SECOND_BATCH[@]} -gt 0 ]; then
        tmux new-window -t $SESSION_NAME -n "second_batch"
        tmux send-keys -t $SESSION_NAME:second_batch "zsh" C-m
        tmux send-keys -t $SESSION_NAME:second_batch "cd /g/korbel2/weber/workspace/StrandSeq_workspace/DEV/mosaicatcher-pipeline-friendsofstrandseq" C-m
        tmux send-keys -t $SESSION_NAME:second_batch "conda activate snakemake-v7-19-1" C-m
        tmux send-keys -t $SESSION_NAME:second_batch "snakemake \
            --config \
                data_location=${SECOND_BATCH_DIR} \
                samples_to_process=${SECOND_BATCH_STRING} \
                multistep_normalisation=False \
                MultiQC=False \
                genome_browsing_files_generation=False \
                ashleys_pipeline=True \
                bypass_ashleys=True \
            --rerun-triggers mtime \
            --rerun-incomplete \
            --profile workflow/snakemake_profiles/mosaicatcher-pipeline/v7/HPC/slurm_EMBL/ \
            --set-resources \
                ashleys_mark_duplicates:constraint="rome" \
                ashleys_generate_features:time="10:00:00" \
                ashleys_bwa_strandseq_to_reference_alignment:time="24:00:00" \
            --jobs 500 \
            -k" C-m    
    fi

    # Attach to the tmux session after sending commands
    tmux attach -t $SESSION_NAME
fi

