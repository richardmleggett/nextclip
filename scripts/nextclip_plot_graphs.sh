#!/bin/bash

# Script:  nextclip_plot_graphs.sh
# Purpose: Plot insert length and read length graphs for NextClip
# Author:  Richard Leggett
# Contact: richard.leggett@tgac.ac.uk
#          http://www.tgac.ac.uk/richard-leggett/

# Location of NextClip scripts - change as appropriate
scriptdir=/Users/leggettr/Downloads/nextclip-master/scripts

libdir=$1
lib=$2

readsdir=${libdir}/reads
graphsdir=${libdir}/graphs
analysisdir=${libdir}/analysis

inputfile=${readsdir}/${lib}_duplicates.txt
outputfile=${graphsdir}/${lib}_duplicates.pdf
echo "Plotting duplication rate graph"
echo " Input: ${inputfile}"
echo "Output: ${outputfile}"

if [ -f ${inputfile} ] ; then
    Rscript ${scriptdir}/nextclip_plot_duplication.R ${inputfile} ${outputfile}
else
    echo "ERROR: Can't find ${inputfile}"
fi

for read in 1 2
do
    inputfile=${readsdir}/${lib}_R${read}_gc.txt
    outputfile=${graphsdir}/${lib}_R${read}_gc.pdf
    echo "Plotting GC content for read ${read}"
    echo " Input: ${inputfile}"
    echo "Output: ${outputfile}"

    if [ -f ${inputfile} ] ; then
        Rscript ${scriptdir}/nextclip_plot_gc.R ${inputfile} ${outputfile}
    else
        echo "ERROR: Can't find ${inputfile}"
    fi
done

for type in A B C D
do
    inputfile=${readsdir}/${lib}_${type}_pair_hist.txt
    outputfile=${graphsdir}/${lib}_cumulative_pairs_${type}.pdf
    echo "Plotting cumulative pair lengths for category ${type}"
    echo " Input: ${inputfile}"
    echo "Output: ${outputfile}"

    if [ -f ${inputfile} ] ; then
        Rscript ${scriptdir}/nextclip_plot_pair_lengths.R ${inputfile} ${outputfile}
    else
        echo "ERROR: Can't find ${inputfile}"
    fi

    for read in 1 2
    do
        inputfile=${readsdir}/${lib}_${type}_R${read}_hist.txt
        outputfile=${graphsdir}/${lib}_lengths_${type}_R${read}.pdf
        echo "Plotting lengths for category ${type} read ${read}"
        echo " Input: ${inputfile}"
        echo "Output: ${outputfile}"

        if [ -f ${inputfile} ] ; then
            Rscript ${scriptdir}/nextclip_plot_lengths.R ${inputfile} ${outputfile}
        else
            echo "ERROR: Can't find ${inputfile}"
        fi
    done

    for kind in mp pe tandem
    do
        inputfile=${analysisdir}/${lib}_${type}_${kind}.txt
        outputfile=${graphsdir}/${lib}_${kind}_${type}.pdf
        echo "Plotting insert size lengths for category ${type} kind ${kind}"
        echo " Input: ${inputfile}"
        echo "Output: ${outputfile}"

        if [ -f ${inputfile} ] ; then
            Rscript ${scriptdir}/nextclip_plot_inserts.R ${inputfile} ${outputfile}
        else
            echo "ERROR: Can't find ${inputfile}"
        fi
    done
done

