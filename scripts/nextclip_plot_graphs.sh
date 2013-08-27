#!/bin/bash

# Script:  nextclip_plot_graphs.sh
# Purpose: Plot insert length and read length graphs for NextClip
# Author:  Richard Leggett
# Contact: richard.leggett@tgac.ac.uk
#          http://www.tgac.ac.uk/richard-leggett/

# Location of NextClip scripts - change as appropriate
scriptdir=/data/workarea/leggettr/matepairs/scripts

lib=$1

echo "Plotting duplication rate graph"
inputfile=${lib}/reads/${lib}_duplicates.txt
outputfile=${lib}/graphs/${lib}_duplicates.pdf

if [ -f ${inputfile} ] ; then
    Rscript ${scriptdir}/nextclip_plot_duplication.R ${inputfile} ${outputfile}
else
    echo "ERROR: Can't find ${inputfile}"
fi

for read in 1 2
do
    echo "Plotting GC content for read ${read}"
    inputfile=${lib}/reads/${lib}_R${read}_gc.txt
    outputfile=${lib}/graphs/${lib}_R${read}_gc.pdf

    if [ -f ${inputfile} ] ; then
        Rscript ${scriptdir}/nextclip_plot_gc.R ${inputfile} ${outputfile}
    else
        echo "ERROR: Can't find ${inputfile}"
    fi
done

for type in A B C D
do
    echo "Plotting cumulative pair lengths for category ${type}"
    inputfile=${lib}/reads/${lib}_${type}_pair_hist.txt
    outputfile=${lib}/graphs/${lib}_cumulative_pairs_${type}.pdf

    if [ -f ${inputfile} ] ; then
        Rscript ${scriptdir}/nextclip_plot_pair_lengths.R ${inputfile} ${outputfile}
    else
        echo "ERROR: Can't find ${inputfile}"
    fi

    for read in 1 2
    do
        echo "Plotting lengths for category ${type} read ${read}"
        inputfile=${lib}/reads/${lib}_${type}_R${read}_hist.txt
        outputfile=${lib}/graphs/${lib}_lengths_${type}_R${read}.pdf

        if [ -f ${inputfile} ] ; then
            Rscript ${scriptdir}/nextclip_plot_lengths.R ${inputfile} ${outputfile}
        else
            echo "ERROR: Can't find ${inputfile}"
        fi
    done

    for kind in mp pe tandem
    do
        echo "Plotting insert size lengths for category ${type} kind ${kind}"
        inputfile=${lib}/analysis/${lib}_${type}_${kind}.txt
        outputfile=${lib}/graphs/${lib}_${kind}_${type}.pdf

        if [ -f ${inputfile} ] ; then
            Rscript ${scriptdir}/nextclip_plot_inserts.R ${inputfile} ${outputfile}
        else
            echo "ERROR: Can't find ${inputfile}"
        fi
    done
done

