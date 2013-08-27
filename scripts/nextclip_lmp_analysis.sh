#!/bin/bash

# Script:  nextclip_lmp_analysis.sh
# Purpose: NextClip Long Mate Pair analysis pipeline
# Author:  Richard Leggett
# Contact: richard.leggett@tgac.ac.uk
#          http://www.tgac.ac.uk/richard-leggett/

# Include any lines here to add paths to BWA, TexLive or R - at TGAC, we use source command
# source bwa-0.6.1
# source texlive-2012
# source R-2.15.1

# Location of NextClip tool - change as appropriate
nextclip=/usr/users/ga002/leggettr/programs/nextclip/bin/nextclip

# Location of NextClip scripts - change as appropriate
scriptdir=/data/workarea/leggettr/matepairs/scripts

# Scheduler currently supported options are "LSF" or "none"
scheduler="LSF"

# Threads to use for BWA
numthreads=8

# Minimum mapping quality
minmapq=10

# ==================================================================================================================================
# Clip transposons
# ==================================================================================================================================
function run_clipping
{
    command="${nextclip} --input_one ${readsdir}/${read_one} --input_two ${readsdir}/${read_two} --output_prefix ${readsdir}/${lib} --log ${logdir}/nextclip_${lib}_log.txt --min_length 25 --number_of_reads 100000000 --trim_ends 0"
    logfile="${logdir}/nextclip_${lib}_lsf.txt"

    if [ "${scheduler}" == "LSF" ] ; then
        bsub -J ${lib}clip -oo ${logfile} -R "rusage[mem=4000]" "${command}"
    elif [ "${scheduler}" == "none" ] ; then
        eval ${command} | tee ${logfile}
    fi
}

# ==================================================================================================================================
# Generate alignments with BWA
# ==================================================================================================================================
function run_alignment
{
    for type in A B C D
    do
        for read in 1 2
        do
            fastqfile=${readsdir}/${lib}_${type}_R${read}.fastq
            saifile=${lib}/bwa/${lib}_${type}_R${read}.sai
            samfile=${lib}/bwa/${lib}_${type}_R${read}.sam

            if [ ${read_length} -gt 101 ] ; then
                command="bwa bwasw -M -t ${numthreads} ${reference} ${fastqfile} > ${samfile}"
                logfile="${logdir}/sam_${type}_R${read}_lsf.txt"
                if [ "${scheduler}" == "LSF" ] ; then
                    bsub -J ${lib}_sam_${type}_R${read} -w "ended(${lib}clip)" -oo ${logfile} -R "span[hosts=1]" -n ${numthreads} "${command}"
                elif [ "${scheduler}" == "none" ] ; then
                    eval ${command} 2>&1 | tee ${logfile}
                fi
            else
                command="bwa aln -t ${numthreads} ${reference} ${fastqfile} > ${saifile}"
                logfile="${logdir}/sai_${type}_R${read}_lsf.txt"
                if [ "${scheduler}" == "LSF" ] ; then
                    bsub -J ${lib}_sai_${type}_R${read} -w "ended(${lib}clip)" -oo ${logfile} -R "span[hosts=1]" -n ${numthreads} "${command}"
                elif [ "${scheduler}" == "none" ] ; then
                    eval ${command} 2>&1 | tee ${logfile}
                fi

                command="bwa samse ${reference} ${saifile} ${fastqfile} > ${samfile}"
                logfile="${logdir}/sam_${type}_R${read}_lsf.txt"
                if [ "${scheduler}" == "LSF" ] ; then
                    bsub -J ${lib}_sam_${type}_R${read} -w "ended(${lib}_sai_${type}_R${read})" -oo ${logfile} "${command}"
                elif [ "${scheduler}" == "none" ] ; then
                    eval ${command} 2>&1 | tee ${logfile}
                fi
            fi
        done
    done
}

# ==================================================================================================================================
# Parse SAM files
# ==================================================================================================================================
function parse_alignment
{
    for type in A B C D
    do
        sam_one=${lib}/bwa/${lib}_${type}_R1.sam
        sam_two=${lib}/bwa/${lib}_${type}_R2.sam
        out_prefix=${lib}/analysis/${lib}_${type}

        command="perl ${scriptdir}/nextclip_sam_parse.pl -one ${sam_one} -two ${sam_two} -out ${out_prefix} -reference ${reference} -refminsize ${ref_min_size} -minmapq ${minmapq}"
        logfile="${logdir}/parse_${type}_lsf.txt"
        if [ "${scheduler}" == "LSF" ] ; then
            bsub -J ${lib}_parse_${type} -w "ended(${lib}_sam_${type}*)" -oo ${logfile} "${command}"
        elif [ "${scheduler}" == "none" ] ; then
            eval ${command} | tee ${logfile}
        fi
    done
}

# ==================================================================================================================================
# Plot graphs
# ==================================================================================================================================
function plot_graphs
{
    command="${scriptdir}/nextclip_plot_graphs.sh ${lib}"
    logfile="${logdir}/plot_graphs_${lib}_lsf.txt"

    if [ "${scheduler}" == "LSF" ] ; then
        # TGAC only - use specific hosts for R - comment out if not at TGAC
        # bsub -J ${lib}_plot_graphs -w "ended(${lib}_parse*)" -oo ${logfile} -R "select[hname=='cn-128-23' || hname=='cn-128-03']" "${command}"
        # General line - use outside TGAC
        bsub -J ${lib}_plot_graphs -w "ended(${lib}_parse*)" -oo ${logdir}/plot_graphs_${lib}_lsf.txt "${command}"
    elif [ "${scheduler}" == "none" ] ; then
        eval ${command} | tee ${logfile}
    fi
}

# ==================================================================================================================================
# Latex generation
# ==================================================================================================================================
function make_report
{
    command="${scriptdir}/nextclip_make_report.sh -d ${lib} -l ${lib} -r \"${reference}\" -o \"${organism}\""
    logfile="${logdir}/${lib}_latex_lsf.txt"

    if [ "${scheduler}" == "LSF" ] ; then
        bsub -J ${lib}_latex -w "ended(${lib}_plot_graphs)" -oo ${logfile} "${command}"
    elif [ "${scheduler}" == "none" ] ; then
        eval "${command}" | tee ${logfile}
    fi
}

# ==================================================================================================================================
# Main script
# ==================================================================================================================================

echo ""
echo "---------- Launching LMP analysis... ----------"
echo "  Library: ${lib}"
echo "   Read 1: ${read_one}"
echo "   Read 2: ${read_two}"
echo " Organism: ${organism}"
echo "Reference: ${reference}"
echo ""

logdir=${lib}/logs
readsdir=${lib}/reads

read_one_clipped=${read_one}.clipped
read_two_clipped=${read_two}.clipped
read_one_noadapt=${read_one}.noadaptor
read_two_noadapt=${read_two}.noadaptor
summaryfile=${lib}/summary.txt

if [ ! -d ${lib} ] ; then
    echo "Making directory ${lib}"
    mkdir ${lib}
fi

for dir in bwa reads logs graphs analysis latex
do
    if [ ! -d ${lib}/${dir} ] ; then
        mkdir ${lib}/${dir}
    fi
done

read_length=`head -n 4 ${lib}/reads/${read_one} | wc -L`
machine=`head -n 1 ${lib}/reads/${read_one} | awk -F':' '{print $1}'`

echo "Read length is ${read_length}"
echo "Machine ID is ${machine}"

run_clipping
run_alignment
parse_alignment
plot_graphs
make_report
