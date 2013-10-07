#!/bin/bash

# Script:  nextclip_make_report.sh
# Purpose: Produce PDF report of NextClip analysis
# Author:  Richard Leggett
# Contact: richard.leggett@tgac.ac.uk
#          http://www.tgac.ac.uk/richard-leggett/

reference="Unknown"
organism="Unknown"

# ----------------------------------------------------------------------
# Function: output_latex_header
# Purpose:  Ouput preamble to LaTeX file
# ----------------------------------------------------------------------
function output_latex_header
{
    echo "\\documentclass[a4paper,11pt,oneside]{article}" > ${latex_file}
    echo "\\usepackage{graphicx}" >> ${latex_file}
    echo "\\usepackage{url}" >> ${latex_file}
    echo "\\usepackage{subcaption}" >> ${latex_file}
    echo "\\usepackage{rotating}" >> ${latex_file}
    echo "\\usepackage{color}" >> ${latex_file}
    echo "\\usepackage[compact]{titlesec}" >> ${latex_file}
    echo "\\usepackage[portrait,top=1cm, bottom=2cm, left=1cm, right=1cm]{geometry}" >> ${latex_file}
    echo "\\begin{document}" >> ${latex_file}
    echo "\\renewcommand*{\familydefault}{\sfdefault}" >> ${latex_file}
}

# ----------------------------------------------------------------------
# Function: output_latex_footer
# Purpose:  End LaTeX document
# ----------------------------------------------------------------------
function output_latex_footer
{
    echo "\\end{document}" >> ${latex_file}
}

# ----------------------------------------------------------------------
# Function: make_pdf
# Purpose:  Run pdflatex
# ----------------------------------------------------------------------
function make_pdf
{
    echo "--- Making PDF ---"
    pdflatex --file-line-error --shell-escape --output-directory ${output_dir} --interaction=batchmode ${latex_file}
    pdf_file=`echo ${latex_file} | sed -e s/[.]tex/\.pdf/`
    echo "Written PDF file ${pdf_file} from ${latex_file}"
}

# ----------------------------------------------------------------------
# Function: output_category_report
# Purpose:  Output LaTeX for each category (A, B, C etc.)
# ----------------------------------------------------------------------
function output_category_report
{
    sumfile=$1
    suffix=$2

    mp_limit=`cat ${sumfile} | grep "MP limit" | perl -nae '$_ =~ /MP limit:\t(\d+)/ ; print $1'`
    pe_limit=`cat ${sumfile} | grep "PE limit" | perl -nae '$_ =~ /PE limit:\t(\d+)/ ; print $1'`
    tandem_limit=`cat ${sumfile} | grep "Tandem limit" | perl -nae '$_ =~ /Tandem limit:\t(\d+)/ ; print $1'`
    mapq_threshold=`cat ${sumfile} | grep "MAPQ threshold" | perl -nae '$_ =~ /MAPQ threshold:\t(\d+)/ ; print $1'`
    ref_min=`cat ${sumfile} | grep "Reference minimum contig size" | perl -nae '$_ =~ /Reference minimum contig size:\t(\d+)/ ; print $1'`
    num_alignments=`cat ${sumfile} | grep "Number of alignments" | perl -nae '$_ =~ /Number of alignments:\t(\d+)/ ; print $1'`

    n_mp=`cat ${sumfile} | grep "Number of mate pair" | head -n 1 | perl -nae '$_ =~ /Number of mate pair:\t(\d+)\t(\S+)\t(\S+)\t(\S+)/ ; print $1'`
    pc_mp=`cat ${sumfile} | grep "Number of mate pair" | head -n 1 | perl -nae '$_ =~ /Number of mate pair:\t(\d+)\t(\S+)\t(\S+)\t(\S+)/ ; print $2'`
    mp_r1_av=`cat ${sumfile} | grep "Number of mate pair" | head -n 1 | perl -nae '$_ =~ /Number of mate pair:\t(\d+)\t(\S+)\t(\S+)\t(\S+)/ ; print $3'`
    mp_r2_av=`cat ${sumfile} | grep "Number of mate pair" | head -n 1 | perl -nae '$_ =~ /Number of mate pair:\t(\d+)\t(\S+)\t(\S+)\t(\S+)/ ; print $4'`
    n_mp_in_range=`cat ${sumfile} | grep "Number of mate pair in limit" | perl -nae '$_ =~ /Number of mate pair in limit:\t(\d+)\t(\S+)/ ; print $1'`
    pc_mp_in_range=`cat ${sumfile} | grep "Number of mate pair in limit" | perl -nae '$_ =~ /Number of mate pair in limit:\t(\d+)\t(\S+)/ ; print $2'`
    n_mp_out_of_range=`cat ${sumfile} | grep "Number of mate pair out of limit" | perl -nae '$_ =~ /Number of mate pair out of limit:\t(\d+)\t(\S+)/ ; print $1'`
    pc_mp_out_of_range=`cat ${sumfile} | grep "Number of mate pair out of limit" | perl -nae '$_ =~ /Number of mate pair out of limit:\t(\d+)\t(\S+)/ ; print $2'`

    n_pe=`cat ${sumfile} | grep "Number of pair end" | head -n 1 | perl -nae '$_ =~ /Number of pair end:\t(\d+)\t(\S+)\t(\S+)\t(\S+)/ ; print $1'`
    pc_pe=`cat ${sumfile} | grep "Number of pair end" | head -n 1 | perl -nae '$_ =~ /Number of pair end:\t(\d+)\t(\S+)\t(\S+)\t(\S+)/ ; print $2'`
    pe_r1_av=`cat ${sumfile} | grep "Number of pair end" | head -n 1 | perl -nae '$_ =~ /Number of pair end:\t(\d+)\t(\S+)\t(\S+)\t(\S+)/ ; print $3'`
    pe_r2_av=`cat ${sumfile} | grep "Number of pair end" | head -n 1 | perl -nae '$_ =~ /Number of pair end:\t(\d+)\t(\S+)\t(\S+)\t(\S+)/ ; print $4'`
    n_pe_in_range=`cat ${sumfile} | grep "Number of pair end in limit" | perl -nae '$_ =~ /Number of pair end in limit:\t(\d+)\t(\S+)/ ; print $1'`
    pc_pe_in_range=`cat ${sumfile} | grep "Number of pair end in limit" | perl -nae '$_ =~ /Number of pair end in limit:\t(\d+)\t(\S+)/ ; print $2'`
    n_pe_out_of_range=`cat ${sumfile} | grep "Number of pair end out of limit" | perl -nae '$_ =~ /Number of pair end out of limit:\t(\d+)\t(\S+)/ ; print $1'`
    pc_pe_out_of_range=`cat ${sumfile} | grep "Number of pair end out of limit" | perl -nae '$_ =~ /Number of pair end out of limit:\t(\d+)\t(\S+)/ ; print $2'`

    n_tandem=`cat ${sumfile} | grep "Number of tandem" | head -n 1 | perl -nae '$_ =~ /Number of tandem:\t(\d+)\t(\S+)\t(\S+)\t(\S+)/ ; print $1'`
    pc_tandem=`cat ${sumfile} | grep "Number of tandem" | head -n 1 | perl -nae '$_ =~ /Number of tandem:\t(\d+)\t(\S+)\t(\S+)\t(\S+)/ ; print $2'`
    tandem_r1_av=`cat ${sumfile} | grep "Number of tandem" | head -n 1 | perl -nae '$_ =~ /Number of tandem:\t(\d+)\t(\S+)\t(\S+)\t(\S+)/ ; print $3'`
    tandem_r2_av=`cat ${sumfile} | grep "Number of tandem" | head -n 1 | perl -nae '$_ =~ /Number of tandem:\t(\d+)\t(\S+)\t(\S+)\t(\S+)/ ; print $4'`
    n_tandem_in_range=`cat ${sumfile} | grep "Number of tandem in limit" | perl -nae '$_ =~ /Number of tandem in limit:\t(\d+)\t(\S+)/ ; print $1'`
    pc_tandem_in_range=`cat ${sumfile} | grep "Number of tandem in limit" | perl -nae '$_ =~ /Number of tandem in limit:\t(\d+)\t(\S+)/ ; print $2'`
    n_tandem_out_of_range=`cat ${sumfile} | grep "Number of tandem out of limit" | perl -nae '$_ =~ /Number of tandem out of limit:\t(\d+)\t(\S+)/ ; print $1'`
    pc_tandem_out_of_range=`cat ${sumfile} | grep "Number of tandem out of limit" | perl -nae '$_ =~ /Number of tandem out of limit:\t(\d+)\t(\S+)/ ; print $2'`

    n_unmapped=`cat ${sumfile} | grep "Number with both unmapped" | perl -nae '$_ =~ /Number with both unmapped:\t(\d+)\t(\S+)/ ; print $1'`
    pc_unmapped=`cat ${sumfile} | grep "Number with both unmapped" | perl -nae '$_ =~ /Number with both unmapped:\t(\d+)\t(\S+)/ ; print $2'`
    n_singlemapped=`cat ${sumfile} | grep "Number with one read unmapped" | perl -nae '$_ =~ /Number with one read unmapped:\t(\d+)\t(\S+)/ ; print $1'`
    pc_singlemapped=`cat ${sumfile} | grep "Number with one read unmapped" | perl -nae '$_ =~ /Number with one read unmapped:\t(\d+)\t(\S+)/ ; print $2'`
    n_differentmapped=`cat ${sumfile} | grep "Number that map differently" | perl -nae '$_ =~ /Number that map differently:\t(\d+)\t(\S+)/ ; print $1'`
    pc_differentmapped=`cat ${sumfile} | grep "Number that map differently" | perl -nae '$_ =~ /Number that map differently:\t(\d+)\t(\S+)/ ; print $2'`
    n_poormapping=`cat ${sumfile} | grep "Number with poor mapping" | perl -nae '$_ =~ /Number with poor mapping:\t(\d+)\t(\S+)/ ; print $1'`
    pc_poormapping=`cat ${sumfile} | grep "Number with poor mapping" | perl -nae '$_ =~ /Number with poor mapping:\t(\d+)\t(\S+)/ ; print $2'`

    good_maps=`cat ${sumfile} | grep "Total number with good mapping" | perl -nae '$_ =~ /Total number with good mapping:\t(\d+)\t(\S+)/ ; print $1, " (", $2, "\\\%)"'`
    good_maps_in_range=`cat ${sumfile} | grep "Total good maps in limit" | perl -nae '$_ =~ /Total good maps in limit:\t(\d+)\t(\S+)/ ; print $1, " (", $2, "\\\%)"'`
    good_maps_out_of_range=`cat ${sumfile} | grep "Total good maps out of limit" | perl -nae '$_ =~ /Total good maps out of limit:\t(\d+)\t(\S+)/ ; print $1, " (", $2, "\\\%)"'`
    bad_maps=`cat ${sumfile} | grep "Total number with bad mapping" | perl -nae '$_ =~ /Total number with bad mapping:\t(\d+)\t(\S+)/ ; print $1, " (", $2, "\\\%)"'`

    echo "\\begin{table}[h!]" >> ${latex_file}
    echo "{\\footnotesize" >> ${latex_file}
    echo "\\fontsize{9pt}{11pt}\selectfont" >> ${latex_file}
    echo "\\begin{tabular}{l c c c c c}" >> ${latex_file}
    echo "& {\bf Total} & {\bf In range} & {\bf Out of range} & {\bf R1 mean} & {\bf R2 mean} \\\\" >> ${latex_file}
    echo "{\bf Pairs producing good mappings} & ${good_maps} & ${good_maps_in_range} & ${good_maps_out_of_range} & & \\\\" >> ${latex_file}
    echo "\hspace{5mm}In Mate Pair orientation & ${n_mp} (${pc_mp}\%) & ${n_mp_in_range} (${pc_mp_in_range}\%) & ${n_mp_out_of_range} (${pc_mp_out_of_range}\%) & ${mp_r1_av} & ${mp_r2_av} \\\\" >> ${latex_file}
    echo "\hspace{5mm}In Pair End orientation & ${n_pe} (${pc_pe}\%) & ${n_pe_in_range} (${pc_pe_in_range}\%) & ${n_pe_out_of_range} (${pc_pe_out_of_range}\%) & ${pe_r1_av} & ${pe_r2_av} \\\\" >> ${latex_file}
    echo "\hspace{5mm}In Tandem orientation & ${n_tandem} (${pc_tandem}\%) & ${n_tandem_in_range} (${pc_tandem_in_range}\%) & ${n_tandem_out_of_range} (${pc_tandem_out_of_range}\%) & ${tandem_r1_av} & ${tandem_r2_av} \\\\" >> ${latex_file}
    echo "{\bf Pairs producing bad mappings} & ${bad_maps} & & & \\\\" >> ${latex_file}
    echo "\hspace{5mm}Both reads unmapped & ${n_unmapped} (${pc_unmapped}\%) & & & \\\\" >> ${latex_file}
    echo "\hspace{5mm}One read unmapped & ${n_singlemapped} (${pc_singlemapped}\%) & & & \\\\" >> ${latex_file}
    echo "\hspace{5mm}Reads map to different ID & ${n_differentmapped} (${pc_differentmapped}\%) & & & \\\\" >> ${latex_file}
    echo "\hspace{5mm}Reads with MAPQ {\textless} ${mapq_threshold} & ${n_poormapping} (${pc_poormapping}\%) & & & \\\\" >> ${latex_file}
    echo "\\end{tabular}" >> ${latex_file}
    echo "}" >> ${latex_file}
    echo "\\end{table}" >> ${latex_file}
    echo "\\vspace{-5mm}" >> ${latex_file}

    echo "\\begin{figure}[h!]" >> ${latex_file}
    echo "\\centering" >> ${latex_file}
    echo "\\begin{subfigure}{.3\textwidth}" >> ${latex_file}
    echo "\\centering" >> ${latex_file}
    echo "\\includegraphics[width=.9\linewidth]{${librarydir}/graphs/${libname}_mp_${suffix}.pdf}" >> ${latex_file}
    echo "\\caption*{\\footnotesize Mate pair}" >> ${latex_file}
    echo "\\end{subfigure}" >> ${latex_file}
    echo "\\begin{subfigure}{.3\textwidth}" >> ${latex_file}
    echo "\\centering" >> ${latex_file}
    echo "\\includegraphics[width=.9\linewidth]{${librarydir}/graphs/${libname}_pe_${suffix}.pdf}" >> ${latex_file}
    echo "\\caption*{\\footnotesize Pair end}" >> ${latex_file}
    echo "\\end{subfigure}" >> ${latex_file}
    echo "\\begin{subfigure}{.3\textwidth}" >> ${latex_file}
    echo "\\centering" >> ${latex_file}
    echo "\\includegraphics[width=.9\linewidth]{${librarydir}/graphs/${libname}_tandem_${suffix}.pdf}" >> ${latex_file}
    echo "\\caption*{\\footnotesize Tandem}" >> ${latex_file}
    echo "\\end{subfigure}" >> ${latex_file}
    echo "\\end{figure}" >> ${latex_file}
}

# ----------------------------------------------------------------------
# Function: get_main_stats
# Purpose:  Get main NextClip stats
# ----------------------------------------------------------------------
function get_main_stats
{
    nextclip_output=$librarydir/logs/nextclip_${libname}_lsf.txt
    minimum_read_size=`grep 'Minimum read size' ${nextclip_output} | perl -nae '$_ =~ /Minimum read size: (\d+)/ ; print $1'`
    number_of_pairs=`grep 'Number of read pairs' ${nextclip_output} | perl -nae '$_ =~ /Number of read pairs: (\d+)/ ; print $1'`

    r1_with_adaptor=`grep 'R1 Num reads with adaptor' ${nextclip_output} | perl -nae '$_ =~ /R1 Num reads with adaptor: (\d+)\t(\S+)/ ; print $1, " (", $2, "\\\%)"'`
    r1_with_adaptor_and_external=`grep 'R1 Num with external also' ${nextclip_output} | perl -nae '$_ =~ /R1 Num with external also: (\d+)\t(\S+)/ ; print $1, " (", $2, "\\\%)"'`
    r1_with_adaptor_long=`grep 'R1 long adaptor reads' ${nextclip_output} | perl -nae '$_ =~ /R1 long adaptor reads: (\d+)\t(\S+)/ ; print $1, " (", $2, "\\\%)"'`
    r1_with_adaptor_short=`grep 'R1 reads too short' ${nextclip_output} | perl -nae '$_ =~ /R1 reads too short: (\d+)\t(\S+)/ ; print $1, " (", $2, "\\\%)"'`
    r1_without_adaptor=`grep 'R1 Num reads no adaptor' ${nextclip_output} | perl -nae '$_ =~ /R1 Num reads no adaptor: (\d+)\t(\S+)/ ; print $1, " (", $2, "\\\%)"'`
    r1_without_adaptor_but_external=`grep 'R1 no adaptor but external' ${nextclip_output} | perl -nae '$_ =~ /R1 no adaptor but external: (\d+)\t(\S+)/ ; print $1, " (", $2, "\\\%)"'`

    r2_with_adaptor=`grep 'R2 Num reads with adaptor' ${nextclip_output} | perl -nae '$_ =~ /R2 Num reads with adaptor: (\d+)\t(\S+)/ ; print $1, " (", $2, "\\\%)"'`
    r2_with_adaptor_and_external=`grep 'R2 Num with external also' ${nextclip_output} | perl -nae '$_ =~ /R2 Num with external also: (\d+)\t(\S+)/ ; print $1, " (", $2, "\\\%)"'`
    r2_with_adaptor_long=`grep 'R2 long adaptor reads' ${nextclip_output} | perl -nae '$_ =~ /R2 long adaptor reads: (\d+)\t(\S+)/ ; print $1, " (", $2, "\\\%)"'`
    r2_with_adaptor_short=`grep 'R2 reads too short' ${nextclip_output} | perl -nae '$_ =~ /R2 reads too short: (\d+)\t(\S+)/ ; print $1, " (", $2, "\\\%)"'`
    r2_without_adaptor=`grep 'R2 Num reads no adaptor' ${nextclip_output} | perl -nae '$_ =~ /R2 Num reads no adaptor: (\d+)\t(\S+)/ ; print $1, " (", $2, "\\\%)"'`
    r2_without_adaptor_but_external=`grep 'R2 no adaptor but external' ${nextclip_output} | perl -nae '$_ =~ /R2 no adaptor but external: (\d+)\t(\S+)/ ; print $1, " (", $2, "\\\%)"'`

    cat_a_total=`grep 'Total pairs in category A' ${nextclip_output} | perl -nae '$_ =~ /Total pairs in category A: (\d+)\t(\S+)/ ; print $1, " (", $2, "\\\%)"'`
    cat_a_long=`grep 'A pairs long enough' ${nextclip_output} | perl -nae '$_ =~ /A pairs long enough: (\d+)\t(\S+)/ ; print $1, " (", $2, "\\\%)"'`
    cat_a_short=`grep 'A pairs too short' ${nextclip_output} | perl -nae '$_ =~ /A pairs too short: (\d+)\t(\S+)/ ; print $1, " (", $2, "\\\%)"'`
    cat_a_external=`grep 'A external clip in 1 or both' ${nextclip_output} | perl -nae '$_ =~ /A external clip in 1 or both: (\d+)\t(\S+)/ ; print $1, " (", $2, "\\\%)"'`
    cat_b_total=`grep 'Total pairs in category B' ${nextclip_output} | perl -nae '$_ =~ /Total pairs in category B: (\d+)\t(\S+)/ ; print $1, " (", $2, "\\\%)"'`
    cat_b_long=`grep 'B pairs long enough' ${nextclip_output} | perl -nae '$_ =~ /B pairs long enough: (\d+)\t(\S+)/ ; print $1, " (", $2, "\\\%)"'`
    cat_b_short=`grep 'B pairs too short' ${nextclip_output} | perl -nae '$_ =~ /B pairs too short: (\d+)\t(\S+)/ ; print $1, " (", $2, "\\\%)"'`
    cat_b_external=`grep 'B external clip in 1 or both' ${nextclip_output} | perl -nae '$_ =~ /B external clip in 1 or both: (\d+)\t(\S+)/ ; print $1, " (", $2, "\\\%)"'`
    cat_c_total=`grep 'Total pairs in category C' ${nextclip_output} | perl -nae '$_ =~ /Total pairs in category C: (\d+)\t(\S+)/ ; print $1, " (", $2, "\\\%)"'`
    cat_c_long=`grep 'C pairs long enough' ${nextclip_output} | perl -nae '$_ =~ /C pairs long enough: (\d+)\t(\S+)/ ; print $1, " (", $2, "\\\%)"'`
    cat_c_short=`grep 'C pairs too short' ${nextclip_output} | perl -nae '$_ =~ /C pairs too short: (\d+)\t(\S+)/ ; print $1, " (", $2, "\\\%)"'`
    cat_c_external=`grep 'C external clip in 1 or both' ${nextclip_output} | perl -nae '$_ =~ /C external clip in 1 or both: (\d+)\t(\S+)/ ; print $1, " (", $2, "\\\%)"'`
    cat_d_total=`grep 'Total pairs in category D' ${nextclip_output} | perl -nae '$_ =~ /Total pairs in category D: (\d+)\t(\S+)/ ; print $1, " (", $2, "\\\%)"'`
    cat_d_long=`grep 'D pairs long enough' ${nextclip_output} | perl -nae '$_ =~ /D pairs long enough: (\d+)\t(\S+)/ ; print $1, " (", $2, "\\\%)"'`
    cat_d_short=`grep 'D pairs too short' ${nextclip_output} | perl -nae '$_ =~ /D pairs too short: (\d+)\t(\S+)/ ; print $1, " (", $2, "\\\%)"'`
    cat_d_external=`grep 'D external clip in 1 or both' ${nextclip_output} | perl -nae '$_ =~ /D external clip in 1 or both: (\d+)\t(\S+)/ ; print $1, " (", $2, "\\\%)"'`

    all_too_short=`grep 'All categories too short' ${nextclip_output} | perl -nae '$_ =~ /All categories too short: (\d+)\t(\S+)/ ; print $1, " (", $2, "\\\%)"'`
    all_long_enough=`grep 'All long enough' ${nextclip_output} | perl -nae '$_ =~ /All long enough: (\d+)\t(\S+)/ ; print $1, " (", $2, "\\\%)"'`

    total_usable=`grep 'Total usable pairs' ${nextclip_output} | perl -nae '$_ =~ /Total usable pairs: (\d+)\t(\S+)/ ; print $1, " (", $2, "\\\%)"'`

    number_of_duplicates=`grep 'Number of duplicate pairs' ${nextclip_output} | perl -nae '$_ =~ /Number of duplicate pairs: (\d+)\t(\S+)/ ; print $1, " (", $2, "\\\%)"'`
    number_of_pairs_containing_n=`grep 'Number of pairs containing N' ${nextclip_output} | perl -nae '$_ =~ /Number of pairs containing N: (\d+)\t(\S+)/ ; print $1, " (", $2, "\\\%)"'`

    gc_content=`grep 'Overall GC content' ${nextclip_output} | perl -nae '$_ =~ /Overall GC content: (\S+)/ ; print $1, , "\\\%"'`
}

# ----------------------------------------------------------------------
# Function: write_overall_section
# Purpose:  Output overall stats
# ----------------------------------------------------------------------
function write_overall_section
{
    echo "\\subsection*{Overall}" >> ${latex_file};
    echo "\\vspace{-3mm}" >> ${latex_file}

    echo "\\begin{table}[h!]" >> ${latex_file}
    echo "{\\footnotesize" >> ${latex_file}
    echo "\\fontsize{9pt}{11pt}\selectfont" >> ${latex_file}
    echo "\\begin{tabular}{c c c}" >> ${latex_file}
    echo "{\bf R1} & & {\bf R2} \\\\" >> ${latex_file}
    echo "${number_of_pairs} & Number of reads & ${number_of_pairs} \\\\" >> ${latex_file}
    echo "${r1_with_adaptor} & With junction adaptor & ${r2_with_adaptor} \\\\" >> ${latex_file}
    echo "${r1_with_adaptor_long} & of which long enough ($\ge$ ${minimum_read_size}) & ${r2_with_adaptor_long} \\\\" >> ${latex_file}
    echo "${r1_with_adaptor_short} & and too short ($<$ ${minimum_read_size}) & ${r2_with_adaptor_short} \\\\" >> ${latex_file}
    echo "${r1_without_adaptor} & Without junction adaptor & ${r2_without_adaptor} \\\\" >> ${latex_file}
    echo "\\end{tabular}" >> ${latex_file}
    echo "}" >> ${latex_file}
    echo "\\end{table}" >> ${latex_file}
    echo "\\vspace{-5mm}" >> ${latex_file}

    echo "\\begin{table}[h!]" >> ${latex_file}
    echo "{\\footnotesize" >> ${latex_file}
    echo "\\fontsize{9pt}{11pt}\selectfont" >> ${latex_file}
    echo "\\begin{tabular}{l c c c}" >> ${latex_file}
    echo "{\bf Category} & {\bf Number of pairs} & {\bf Too short ($<$ ${minimum_read_size})} & {\bf Long enough ($\ge$ ${minimum_read_size})} \\\\" >> ${latex_file}
    echo "Adaptor in R1 and R2 (A) & ${cat_a_total} & ${cat_a_short} & ${cat_a_long} \\\\" >> ${latex_file}
    echo "Adaptor in R2 only (B) & ${cat_b_total} & ${cat_b_short} & ${cat_b_long} \\\\" >> ${latex_file}
    echo "Adaptor in R1 only (C) & ${cat_c_total} & ${cat_c_short} & ${cat_c_long} \\\\" >> ${latex_file}
    echo "Adaptor in neither (D) & ${cat_d_total} & ${cat_d_short} & ${cat_d_long} \\\\" >> ${latex_file}
    echo "All categories & ${number_of_pairs} (100\%) & ${all_too_short} & ${all_long_enough} \\\\" >> ${latex_file}
    echo "Total usable (A,B,C) & & & ${total_usable} \\\\" >> ${latex_file}
    echo "\\end{tabular}" >> ${latex_file}
    echo "}" >> ${latex_file}
    echo "\\end{table}" >> ${latex_file}
}

# ----------------------------------------------------------------------
# Function: write_gc_content_section
# Purpose:  Output GC stats
# ----------------------------------------------------------------------
function write_gc_content_section
{
    echo "\\subsection*{GC content}" >> ${latex_file};

    echo "\\begin{table}[h!]" >> ${latex_file}
    echo "{\\footnotesize" >> ${latex_file}
    echo "\\fontsize{9pt}{11pt}\selectfont" >> ${latex_file}
    echo "\\begin{tabular}{l c}" >> ${latex_file}
    echo "Overall GC content & ${gc_content} \\\\" >> ${latex_file}
    echo "\\end{tabular}" >> ${latex_file}
    echo "}" >> ${latex_file}
    echo "\\end{table}" >> ${latex_file}

    echo "\\vspace{-3mm}" >> ${latex_file}
    echo "\\begin{figure}[h!]" >> ${latex_file}
    echo "\\centering" >> ${latex_file}
    echo "\\begin{subfigure}{.33\textwidth}" >> ${latex_file}
    echo "\\centering" >> ${latex_file}
    echo "\\includegraphics[width=.9\linewidth]{${librarydir}/graphs/${libname}_R1_gc.pdf}" >> ${latex_file}
    echo "\\caption*{\\footnotesize R1}" >> ${latex_file}
    echo "\\end{subfigure}" >> ${latex_file}
    echo "\\begin{subfigure}{.33\textwidth}" >> ${latex_file}
    echo "\\centering" >> ${latex_file}
    echo "\\includegraphics[width=.9\linewidth]{${librarydir}/graphs/${libname}_R2_gc.pdf}" >> ${latex_file}
    echo "\\caption*{\\footnotesize R2}" >> ${latex_file}
    echo "\\end{subfigure}" >> ${latex_file}
    echo "\\end{figure}" >> ${latex_file}
}

# ----------------------------------------------------------------------
# Function: write_shortest_pair_section
# Purpose:  Output shortest pair stats
# ----------------------------------------------------------------------
function write_shortest_pair_section
{
    echo "\\subsection*{Shortest pair length}" >> ${latex_file};
    echo "\\vspace{-3mm}" >> ${latex_file}
    echo "\\begin{figure}[h!]" >> ${latex_file}
    echo "\\centering" >> ${latex_file}
    echo "\\begin{subfigure}{.22\textwidth}" >> ${latex_file}
    echo "\\centering" >> ${latex_file}
    echo "\\includegraphics[width=.9\linewidth]{${librarydir}/graphs/${libname}_cumulative_pairs_A.pdf}" >> ${latex_file}
    echo "\\caption*{\\footnotesize Category A pairs}" >> ${latex_file}
    echo "\\end{subfigure}" >> ${latex_file}
    echo "\\begin{subfigure}{.22\textwidth}" >> ${latex_file}
    echo "\\centering" >> ${latex_file}
    echo "\\includegraphics[width=.9\linewidth]{${librarydir}/graphs/${libname}_cumulative_pairs_B.pdf}" >> ${latex_file}
    echo "\\caption*{\\footnotesize Category B pairs}" >> ${latex_file}
    echo "\\end{subfigure}" >> ${latex_file}
    echo "\\begin{subfigure}{.22\textwidth}" >> ${latex_file}
    echo "\\includegraphics[width=.9\linewidth]{${librarydir}/graphs/${libname}_cumulative_pairs_C.pdf}" >> ${latex_file}
    echo "\\caption*{\\footnotesize Category C pairs}" >> ${latex_file}
    echo "\\end{subfigure}" >> ${latex_file}
    echo "\\end{figure}" >> ${latex_file}
}

# ----------------------------------------------------------------------
# Function: write_clipped_read_length_section
# Purpose:  Output clipped read length stats
# ----------------------------------------------------------------------
function write_clipped_read_length_section
{
    echo "\\subsection*{Clipped read lengths}" >> ${latex_file};
    echo "\\vspace{-3mm}" >> ${latex_file}
    echo "\\begin{figure}[h!]" >> ${latex_file}
    echo "\\centering" >> ${latex_file}
    echo "\\begin{subfigure}{.22\textwidth}" >> ${latex_file}
    echo "\\centering" >> ${latex_file}
    echo "\\includegraphics[width=.9\linewidth]{${librarydir}/graphs/${libname}_lengths_A_R1.pdf}" >> ${latex_file}
    echo "\\caption*{\\footnotesize Category A R1 lengths}" >> ${latex_file}
    echo "\\end{subfigure}" >> ${latex_file}
    echo "\\begin{subfigure}{.22\textwidth}" >> ${latex_file}
    echo "\\centering" >> ${latex_file}
    echo "\\includegraphics[width=.9\linewidth]{${librarydir}/graphs/${libname}_lengths_A_R2.pdf}" >> ${latex_file}
    echo "\\caption*{\\footnotesize Category A R2 lengths}" >> ${latex_file}
    echo "\\end{subfigure}" >> ${latex_file}
    echo "\\begin{subfigure}{.22\textwidth}" >> ${latex_file}
    echo "\\includegraphics[width=.9\linewidth]{${librarydir}/graphs/${libname}_lengths_B_R2.pdf}" >> ${latex_file}
    echo "\\caption*{\\footnotesize Category B R2 lengths}" >> ${latex_file}
    echo "\\end{subfigure}" >> ${latex_file}
    echo "\\begin{subfigure}{.22\textwidth}" >> ${latex_file}
    echo "\\centering" >> ${latex_file}
    echo "\\includegraphics[width=.9\linewidth]{${librarydir}/graphs/${libname}_lengths_C_R1.pdf}" >> ${latex_file}
    echo "\\caption*{\\footnotesize Category C R1 lengths}" >> ${latex_file}
    echo "\\end{subfigure}" >> ${latex_file}
    echo "\\end{figure}" >> ${latex_file}
}

# ----------------------------------------------------------------------
# Function: write_duplication_section
# Purpose:  Output PCR duplication stats
# ----------------------------------------------------------------------
function write_duplication_section
{
    echo "\\subsection*{Duplication}" >> ${latex_file};
    echo "\\vspace{-3mm}" >> ${latex_file}

    echo "\\begin{table}[h!]" >> ${latex_file}
    echo "{\\footnotesize" >> ${latex_file}
    echo "\\fontsize{9pt}{11pt}\selectfont" >> ${latex_file}
    echo "\\begin{tabular}{l c}" >> ${latex_file}
    echo "Number of read pairs in library & ${number_of_pairs} \\\\" >> ${latex_file}
    echo "PCR duplicates & ${number_of_duplicates} \\\\" >> ${latex_file}
    echo "\\end{tabular}" >> ${latex_file}
    echo "}" >> ${latex_file}
    echo "\\end{table}" >> ${latex_file}

    echo "\\vspace{-5mm}" >> ${latex_file}
    echo "\\begin{figure}[h!]" >> ${latex_file}
    echo "\\centering" >> ${latex_file}
    echo "\\includegraphics[width=.4\linewidth]{${librarydir}/graphs/${libname}_duplicates.pdf}" >> ${latex_file}
    echo "\\end{figure}" >> ${latex_file}
}

# ----------------------------------------------------------------------
# Function: write_external_adaptor_section
# Purpose:  Output external adaptor stats
# ----------------------------------------------------------------------
function write_external_adaptor_section
{
    echo "\\subsection*{External adaptor}" >> ${latex_file};
    echo "\\vspace{-3mm}" >> ${latex_file}
    echo "\\begin{table}[h!]" >> ${latex_file}
    echo "{\\footnotesize" >> ${latex_file}
    echo "\\fontsize{9pt}{11pt}\selectfont" >> ${latex_file}
    echo "\\begin{tabular}{c c c}" >> ${latex_file}
    echo "{\bf R1} & & {\bf R2} \\\\" >> ${latex_file}
    echo "${number_of_pairs} & Number of reads & ${number_of_pairs} \\\\" >> ${latex_file}
    echo "${r1_with_adaptor_and_external} & With junction adaptor and external adaptor & ${r2_with_adaptor_and_external} \\\\" >> ${latex_file}
    echo "${r1_without_adaptor_but_external} & Without junction adaptor, but with external adaptor & ${r2_without_adaptor_but_external} \\\\" >> ${latex_file}
    echo "\\end{tabular}" >> ${latex_file}
    echo "}" >> ${latex_file}
    echo "\\end{table}" >> ${latex_file}
    echo "\\vspace{-5mm}" >> ${latex_file}
    echo "\\begin{table}[h!]" >> ${latex_file}
    echo "{\\footnotesize" >> ${latex_file}
    echo "\\fontsize{9pt}{11pt}\selectfont" >> ${latex_file}
    echo "\\begin{tabular}{l c}" >> ${latex_file}
    echo "Category A pairs also trimmed for external adaptor & ${cat_a_external} \\\\" >> ${latex_file}
    echo "Category B pairs also trimmed for external adaptor & ${cat_b_external} \\\\" >> ${latex_file}
    echo "Category C pairs also trimmed for external adaptor & ${cat_c_external} \\\\" >> ${latex_file}
    echo "Category D pairs also trimmed for external adaptor & ${cat_d_external} \\\\" >> ${latex_file}
    echo "\\end{tabular}" >> ${latex_file}
    echo "}" >> ${latex_file}
    echo "\\end{table}" >> ${latex_file}
}

# ----------------------------------------------------------------------
# Function: write_ambiguous_bases_section
# Purpose:  Output ambiguous base stats
# ----------------------------------------------------------------------
function write_ambiguous_bases_section
{
    echo "\\subsection*{Ambiguous bases}" >> ${latex_file};
    echo "\\vspace{-3mm}" >> ${latex_file}
    echo "\\begin{table}[h!]" >> ${latex_file}
    echo "{\\footnotesize" >> ${latex_file}
    echo "\\fontsize{9pt}{11pt}\selectfont" >> ${latex_file}
    echo "\\begin{tabular}{l c}" >> ${latex_file}
    echo "Number of pairs containing Ns & ${number_of_pairs_containing_n} \\\\" >> ${latex_file}
    echo "\\end{tabular}" >> ${latex_file}
    echo "}" >> ${latex_file}
    echo "\\end{table}" >> ${latex_file}
}

# ----------------------------------------------------------------------
# Function: write_notes_section
# Purpose:  Output notes
# ----------------------------------------------------------------------
function write_notes_section
{
    echo "\\subsection*{Notes}" >> ${latex_file};
    echo "\\vspace{-3mm}" >> ${latex_file}
    echo "\\begin{table}[h!]" >> ${latex_file}
    echo "{\\footnotesize" >> ${latex_file}
    echo "\\fontsize{9pt}{11pt}\selectfont" >> ${latex_file}
    echo "\\begin{tabular}{l c}" >> ${latex_file}
    echo "Minimum contig size for alignment to reference & ${ref_min} \\\\" >> ${latex_file}
    echo "Maximum allowed MP insert & ${mp_limit} \\\\" >> ${latex_file}
    echo "Maximum allowed PE insert & ${pe_limit} \\\\" >> ${latex_file}
    echo "Maximum allowed tandem insert & ${tandem_limit} \\\\" >> ${latex_file}
    echo "\\end{tabular}" >> ${latex_file}
    echo "}" >> ${latex_file}
    echo "\\end{table}" >> ${latex_file}
}

# ----------------------------------------------------------------------
# Main program
# ----------------------------------------------------------------------

# Loop through command line options
while getopts d:l:r:o: OPTION
do
    case $OPTION in
        d) librarydir=$OPTARG;;
        l) libname=$OPTARG;;
        o) organism=$OPTARG;;
        r) reference=`basename $OPTARG | sed 's/_/\\\\_/g'`;;
    esac
done

title=`echo ${libname} | sed 's/_/\\\\_/g'`
if [ "${organism}" != "Unknown" ] ; then
    title="NextClip report for ${title} (${organism})"
fi
output_dir=${librarydir}/latex
latex_file=${librarydir}/latex/${libname}.tex

echo ""
echo "Report builder"
echo ""
echo "Library name is ${libname}"
echo "Output dir is ${output_dir}"
echo "Organism is ${organism}"
echo "Reference is ${reference}"
echo "LaTeX file is ${latex_file}"

get_main_stats
output_latex_header

echo "{\\section*{\large{${title}}}" >> ${latex_file}

write_overall_section

echo "\\subsection*{Fragments with junction adaptor in R1 and R2 (A)}" >> ${latex_file};
echo "\\vspace{-3mm}" >> ${latex_file}
output_category_report ${librarydir}/logs/parse_A_lsf.txt A
echo "\\subsection*{Fragments with junction adaptor in R2 only (B)}" >> ${latex_file};
echo "\\vspace{-3mm}" >> ${latex_file}
output_category_report ${librarydir}/logs/parse_B_lsf.txt B

echo "\\clearpage" >> ${latex_file}

echo "\\subsection*{Fragments with junction adaptor in R1 only (C)}" >> ${latex_file};
echo "\\vspace{-3mm}" >> ${latex_file}
output_category_report ${librarydir}/logs/parse_C_lsf.txt C
echo "\\subsection*{Fragments where neither read contain the junction adaptor (D)}" >> ${latex_file};
echo "\\vspace{-3mm}" >> ${latex_file}
output_category_report ${librarydir}/logs/parse_D_lsf.txt D

write_gc_content_section

echo "\\clearpage" >> ${latex_file}

write_shortest_pair_section
write_clipped_read_length_section
write_duplication_section
write_ambiguous_bases_section
write_external_adaptor_section
write_notes_section
output_latex_footer
make_pdf
