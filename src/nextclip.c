#include <stdlib.h>
#include <stdio.h>
#include <stdint.h>
#include <string.h>
#include <getopt.h>
#include <math.h>
#include "global.h"
#include "binary_kmer.h"
#include "element.h"
#include "hash_table.h"

//#define USE_MULTIPLE_HASHES

/*----------------------------------------------------------------------*
 * Constants
 *----------------------------------------------------------------------*/
#define NEXTCLIP_VERSION "0.4"
#define MAX_PATH_LENGTH 1024
#define NUMBER_OF_CATEGORIES 5
#define SEPARATE_KMER_SIZE 11
#define TOTAL_KMER_SIZE (4 * SEPARATE_KMER_SIZE)
#define MAX_DUPLICATES 100
#define NUMBER_OF_HASHES 3

/*----------------------------------------------------------------------*
 * Structures
 *----------------------------------------------------------------------*/
typedef struct {
    char read_header[1024];
    char read[MAX_READ_LENGTH];
    char quality_header[1024];
    char qualities[MAX_READ_LENGTH];
    int read_size;
} FastQRead;

typedef struct {
    int read_size;
    int score;
    int total_matches;
    double total_identity;
    int total_alignment_length;
    int position;
    int matches[2];
    int mismatches[2];
    int alignment_length[2];
    double identity[2];
    int read_start;
    int read_end;
    int transposon_start;
    int transposon_end;
    int accepted;
} AlignmentResult;

typedef struct {
    int read_length;
    FILE* input_fp[2];
    FILE* output_fp[NUMBER_OF_CATEGORIES][2];
    FILE* log_fp;
    FILE* duplicates_fp;
    char input_filenames[2][MAX_PATH_LENGTH];
    char output_prefix[MAX_PATH_LENGTH];
    char output_filenames[NUMBER_OF_CATEGORIES][MAX_PATH_LENGTH];
    char log_filename[MAX_PATH_LENGTH];
    int num_read_pairs;
    int count_adaptor_found[2];
    int count_no_adaptor[2];
    int count_long_enough[2];
    int count_too_short[2];
    double percent_adaptor_found[2];
    double percent_no_adaptor[2];
    double percent_long_enough[2];
    double percent_too_short[2];
    int count_by_category[NUMBER_OF_CATEGORIES];
    int count_by_category_long_enough[NUMBER_OF_CATEGORIES];
    int count_by_category_too_short[NUMBER_OF_CATEGORIES];
    int count_by_category_relaxed_hit[NUMBER_OF_CATEGORIES];
    double percent_by_category[NUMBER_OF_CATEGORIES];
    double percent_by_category_long_enough[NUMBER_OF_CATEGORIES];
    double percent_by_category_too_short[NUMBER_OF_CATEGORIES];
    double percent_by_category_relaxed_hit[NUMBER_OF_CATEGORIES];
    int total_too_short;
    double percent_total_too_short;
    int total_long_enough;
    int total_usable;
    double percent_usable;
    double percent_total_long_enough;
    int read_length_counts[NUMBER_OF_CATEGORIES][2][MAX_READ_LENGTH];
    int read_pair_length_counts[NUMBER_OF_CATEGORIES][MAX_READ_LENGTH];
    int n_duplicates;
    int n_invalid_for_duplicate;
    double percent_duplicates;
    int duplicates_not_written;
    double percent_duplicates_not_written;
    long int at_bases;
    long int gc_bases;
    int gc_content[2][101];
    double percent_gc;
    long int pairs_containing_n;
    double percent_pairs_containing_n;
} MPStats;

/*----------------------------------------------------------------------*
 * Globals
 *----------------------------------------------------------------------*/
char adaptor[] = "CTGTCTCTTATACACATCT";
char transposon[] = "CTGTCTCTTATACACATCTAGATGTGTATAAGAGACAG";
int minimum_read_size = 25;
int strict_double_match = 34;
int strict_single_match = 18;
int relaxed_double_match = 32;
int relaxed_single_match = 17;
int discarded = 0;
int use_category_e = 0;
int remove_duplicates = 0;
int num_categories = NUMBER_OF_CATEGORIES;
int trim_ends = 19;
int approximate_reads = 20000000;

/*
 * Single hash option algorithm
 *
 * We take a kmer from the start of each read and one from the middle, eg.
 *
 *   R1: TGACTGACTGACTGACTGACTGACTGACTG     R2: TGACTGACTGACTGACTGACTGACTGACTG
 * kmer: AAA           BBB                      aaa           bbb
 *
 * Then we make a hash AAABBBaaabbb.
 *
 *
 * Multiple hash option algorithm
 *
 * We take a kmer from each read, from a third of the way and from two thirds of the way through, eg.
 *
 *   R1: TGACTGACTGACTGACTGACTGACTGACTG     R2: TGACTGACTGACTGACTGACTGACTGACTG
 * kmer: AAA       BBB       CCC                aaa       bbb       ccc
 *
 * Then we make three kmer hashes:
 * kmer_hashes[0] = AAABBBaaabbb
 * kmer_hashes[1] = AAACCCaaaccc
 * kmer_hashes[2] = BBBCCCbbbccc
 *
 * Then, for something to be declared a duplicate, we ask for a match to any of the three hash tables.
 * This allows one of the three kmers to not match - eg. because of miscalled base.
 */

#ifdef USE_MULTIPLE_HASHES
HashTable* kmer_hashes[NUMBER_OF_HASHES];
int kmer_offsets[NUMBER_OF_HASHES][2];
#else
HashTable* duplicate_hash = NULL;
#endif

/*----------------------------------------------------------------------*
 * Function:   initialise_stats
 * Purpose:    Initialise MPStats structure
 * Parameters: stats -> an MPStats structure
 * Returns:    None
 *----------------------------------------------------------------------*/
void initialise_stats(MPStats* stats)
{
    int i, j;
    
    for (i=0; i<2; i++) {
        stats->input_filenames[i][0] = 0;
        stats->count_adaptor_found[i] = 0;
        stats->count_no_adaptor[i] = 0;
        stats->count_too_short[i] = 0;
        stats->count_long_enough[i] = 0;
        stats->input_fp[i] = 0;
    }

    for (i=0; i<NUMBER_OF_CATEGORIES; i++) {
        stats->output_filenames[i][0] = 0;
        stats->output_fp[i][0] = 0;
        stats->output_fp[i][1] = 0;
        stats->count_by_category[i] = 0;
        stats->count_by_category_long_enough[i] = 0;
        stats->count_by_category_too_short[i] = 0;
        stats->count_by_category_relaxed_hit[i] = 0;

        for (j=0; j<MAX_READ_LENGTH; j++) {
            stats->read_length_counts[i][0][j] = 0;
            stats->read_length_counts[i][1][j] = 0;
            stats->read_pair_length_counts[i][j] = 0;
        }
    }

    for (i=0; i<=100; i++) {
        stats->gc_content[0][i] = 0;
        stats->gc_content[1][i] = 0;
    }

    stats->read_length = 0;
    stats->log_filename[0] = 0;
    stats->output_prefix[0] = 0;
    stats->log_fp = 0;
    stats->duplicates_fp = 0;
    stats->num_read_pairs = 0;
    stats->n_duplicates = 0;
    stats->n_invalid_for_duplicate = 0;
    stats->duplicates_not_written = 0;
    stats->total_usable = 0;
    stats->gc_bases = 0;
    stats->at_bases = 0;
    stats->pairs_containing_n = 0;
}

/*----------------------------------------------------------------------*
 * Function:   initialise_result
 * Purpose:    Initialise AlignmentResult structure
 * Parameters: result -> an AlignmentResult structure
 * Returns:    None
 *----------------------------------------------------------------------*/
void initialise_result(AlignmentResult* result)
{
    result->score = -1;
    result->position = -1;
    result->total_matches = 0;
    result->total_alignment_length = 0;
    result->total_identity = 0;
    result->mismatches[0] = 0;
    result->matches[0] = 0;
    result->mismatches[1] = 0;
    result->matches[1] = 0;
    result->read_start= 0;
    result->read_end = 0;
    result->transposon_start = 0;
    result->transposon_end = 0;
    result->accepted = 0;
}

/*----------------------------------------------------------------------*
 * Function:   reverse_compliment
 * Purpose:    Make reverse compliment of sequence.
 * Parameters: fwd -> sequence to produce reverse compliment of,
 *             rev -> resulting reverse compliment.
 * Returns:    Pointer to reverse compliment
 *----------------------------------------------------------------------*/
char* reverse_compliment(char* fwd, char* rev)
{
    int i;
    
    for (i=0; i<strlen(fwd); i++) {
        switch(fwd[strlen(fwd) - i - 1]) {
            case 'A': rev[i] = 'T'; break;
            case 'C': rev[i] = 'G'; break;
            case 'G': rev[i] = 'C'; break;
            case 'T': rev[i] = 'A'; break;
            default: printf("ERROR: Bad nucleotide in %s.\n", fwd); exit(2); break; 
        }
    }
    
    rev[strlen(fwd)] = 0;
    
    return rev;
}

/*----------------------------------------------------------------------*
 * Function:   usage
 * Purpose:    Report program usage.
 * Parameters: None
 * Returns:    None
 *----------------------------------------------------------------------*/
void usage(void)
{
    printf("Clip and analyse Illumina Nextera Long Mate Pair reads\n" \
           "\nSyntax: nextclip [-r | -p] [-i r1.fastq] [-j r2.fastq]\n" \
           "\nOptions:\n" \
           "    [-d | --remove_duplicates] Remove PCR duplicates\n"
           "    [-e | --category_e] Use category E\n"
           "    [-h | --help] This help screen\n" \
           "    [-i | --input_one] Input FASTQ R1 file\n" \
           "    [-j | --input_two] Input FASTQ R2 file\n" \
           "    [-l | --log] Log filename\n" \
           "    [-m | --min_length] Minimum usable read length\n" \
           "    [-n | --number_of_reads] Approximate number of reads (default 20,000,000)\n" \
           "    [-s | --adaptor_sequence] Adaptor sequence - default CTGTCTCTTATACACATCT\n" \
           "    [-o | --output_prefix] Prefix for output files\n" \
           "    [-t | --trim_ends] Trim ends of non-matching category B and C reads by amount (default 19)\n" \
           "    [-x | --strict_match] Strict alignment matches (default '34,18')\n" \
           "    [-y | --relaxed_match] Relaxed alignment matches (default '32,17')\n" \
           "\nComments/suggestions to richard.leggett@tgac.ac.uk\n" \
           "\n");
}

/*----------------------------------------------------------------------*
 * Function:   chomp
 * Purpose:    Remove hidden characters from end of line
 * Parameters: str -> String to change
 * Returns:    None
 *----------------------------------------------------------------------*/
void chomp(char* str)
{
    int i = strlen(str) - 1;
    
    while ((i > 0) && (str[i] < ' ')) {
        str[i--] = 0;
    }    
}

/*----------------------------------------------------------------------*
 * Function:   parse_csv_params
 * Purpose:    
 * Parameters: 
 * Returns:    None
 *----------------------------------------------------------------------*/
void parse_csv_params(char* s, int* a, int* b)
{
    char* comma = strchr(s, ',');
    *comma = 0;
    *a = atoi(s);
    *b = atoi(comma+1);    
}

/*----------------------------------------------------------------------*
 * Function:   parse_command_line
 * Purpose:    Parse command line options
 * Parameters: argc = number of arguments
 *             argv -> array of arguments
 * Returns:    None
 *----------------------------------------------------------------------*/
void parse_command_line(int argc, char* argv[], MPStats* stats)
{
    static struct option long_options[] = {
        {"remove_duplicates", no_argument, NULL, 'd'},
        {"use_category_e", no_argument, NULL, 'e'},
        {"help", no_argument, NULL, 'h'},
        {"input_one", required_argument, NULL, 'i'},
        {"input_two", required_argument, NULL, 'j'},
        {"log", required_argument, NULL, 'l'},
        {"min_length", required_argument, NULL, 'm'},
        {"number_of_reads", required_argument, NULL, 'n'},
        {"output_prefix", required_argument, NULL, 'o'},
        {"adaptor_sequence", required_argument, NULL, 's'},
        {"trim_ends", required_argument, NULL, 't'},
        {"strict_match", required_argument, NULL, 'x'},
        {"relaxed_match", required_argument, NULL, 'y'},
        {0, 0, 0, 0}
    };
    int opt;
    int longopt_index;
    int i;
    
    if (argc == 1) {
        usage();
        exit(0);
    }
    
    while ((opt = getopt_long(argc, argv, "dehi:j:l:m:n:o:s:t:x:y:z:", long_options, &longopt_index)) > 0)
    {
        switch(opt) {
            case 'd':
                remove_duplicates=1;
                break;
            case 'h':
                usage();
                exit(0);
                break;
            case 'i':
                if (optarg==NULL) {
                    printf("Error: [-i | --input_one] option requires an argument.\n");
                    exit(1);
                }
                strcpy(stats->input_filenames[0], optarg);
                break;
            case 'j':
                if (optarg==NULL) {
                    printf("Error: [-j | --input_two] option requires an argument.\n");
                    exit(1);
                }
                strcpy(stats->input_filenames[1], optarg);
                break;
            case 'l':
                if (optarg==NULL) {
                    printf("Error: [-l | --log] option requires an argument.\n");
                    exit(1);
                }
                strcpy(stats->log_filename, optarg);
                break;                
            case 'm':
                if (optarg==NULL) {
                    printf("Error: [-m | --minium_length] option requires an argument.\n");
                    exit(1);
                }
                minimum_read_size = atoi(optarg);
                break;
            case 'n':
                if (optarg==NULL) {
                    printf("Error: [-n | --number_of_reads] option requires an argument.\n");
                    exit(1);
                }
                approximate_reads = atoi(optarg);
                break;
            case 'o':
               if (optarg==NULL) {
                    printf("Error: [-o | --output_prefix] option requires an argument.\n");
                    exit(1);
                }
                strcpy(stats->output_prefix, optarg);
                break;                
            case 'p':
                use_category_e = 1;
                break;
            case 's':
                if (optarg==NULL) {
                    printf("Error: [-s | --adaptor_sequence] option requires an argument.\n");
                    exit(1);
                }
                strcpy(adaptor, optarg);
                break;
            case 't':
                if (optarg==NULL) {
                    printf("Error: [-t | --trim_ends] option requires an argument.\n");
                    exit(1);
                }
                trim_ends = atoi(optarg);
                break;                
            case 'x':
                if (optarg==NULL) {
                    printf("Error: [-x | --strict_match] option requires an argument.\n");
                    exit(1);
                }
                if (strchr(optarg, ',') == 0) {
                    printf("Error: [-x | --strict_match] option is of the format 'a,b'.\n");
                    exit(1);                    
                }
                parse_csv_params(optarg, &strict_double_match, &strict_single_match);
                break;
            case 'y':
                if (optarg==NULL) {
                    printf("Error: [-x | --relaxed_match] option requires an argument.\n");
                    exit(1);
                }
                if (strchr(optarg, ',') == 0) {
                    printf("Error: [-x | --relaxed_match] option is of the format 'a,b'.\n");
                    exit(1);                    
                }
                parse_csv_params(optarg, &relaxed_double_match, &relaxed_single_match);
                break;
            default:
                printf("Error: Unknown option %c\n", opt);
                exit(1);
                break;
        }
    }

    if (use_category_e == 1) {
        num_categories = 5;
    } else {
        num_categories = 4;
    }

    if (stats->output_prefix[0] != 0) {
        for (i=0; i<num_categories; i++) {
            sprintf(stats->output_filenames[i], "%s_%c", stats->output_prefix, i+'A');
        }
    }
    
    for (i=0; i<2; i++) {
        if (stats->input_filenames[i][0] == 0) {
            printf("Error: you must specify two input filenames\n");
            exit(2);
        }
    }
    
    for (i=0; i<num_categories; i++) {
        if (stats->output_filenames[i][0] == 0) {
            printf("Error: you must specify four output filenames\n");
            exit(2);
        }
    }
}

/*----------------------------------------------------------------------*
 * Function:   strict_check
 * Purpose:    
 * Parameters: 
 * Returns:    None
 *----------------------------------------------------------------------*/
void strict_check(AlignmentResult* result)
{
    if (result->total_matches >= strict_double_match) {
        result->accepted = 1;
        return;
    } 
   
    if ((result->matches[0] >= strict_single_match) || (result->matches[1] >= strict_single_match)) {
        result->accepted = 1;
        return;
    }
}


/*----------------------------------------------------------------------*
 * Function:   relaxed_check
 * Purpose:    
 * Parameters: 
 * Returns:    None
 *----------------------------------------------------------------------*/
void relaxed_check(AlignmentResult* result)
{
    if (result->total_matches >= relaxed_double_match) {
        result->accepted = 1;
        return;
    } 
    
    if ((result->matches[0] >= relaxed_single_match) || (result->matches[1] >= relaxed_single_match)) {
        result->accepted = 1;
        return;
    }    
}

/*----------------------------------------------------------------------*
 * Function:   find_transposons
 * Purpose:    
 * Parameters: 
 * Returns:    None
 *----------------------------------------------------------------------*/
void find_transposons(FastQRead* read, AlignmentResult* result)
{
    int x, p;
    int t_length = strlen(transposon);

    // Initialise a result structure to store information
    initialise_result(result);
    result->read_size = read->read_size;
    
    // Start searching for the transposon... x is the position in the read where we start to compare the transposon
    for (x=-t_length+5; x<read->read_size-5; x++) {
        int matches[2] = {0, 0};
        int mismatches[2] = {0, 0};
        int score = 0;
        int r_start = -1;
        int r_end = -1;
        int t_start = -1;
        int t_end = -1;
 
        // Go through each base of transposon, count matches and store start and end of match
        // For interest, we store the 19nt sequence and it's reverse as a part 1 and part 2!
        for (p=0; p<t_length; p++) {
            if (((x+p) >= 0) && ((x+p) < read->read_size)) {
                if (transposon[p] == read->read[x+p]) {
                    matches[p < 19 ? 0:1]++;
                    if (r_start == -1) {
                        r_start = x+p;
                        t_start = p;
                    } else {
                        r_end = x+p;
                        t_end = p;
                    }
                } else {
                    mismatches[p < 19 ? 0:1]++;
                }
            }
        }
        
        // Score is simply matches for part 1 and 2
        score = matches[0] + matches[1];
                    
        // Is this the best score yet?
        if (score > result->score) {
            result->score = score;
            result->matches[0] = matches[0];
            result->mismatches[0] = mismatches[0];
            result->matches[1] = matches[1];
            result->mismatches[1] = mismatches[1];
            result->position = x;
            result->read_start = r_start;
            result->read_end = r_end;
            result->transposon_start = t_start;
            result->transposon_end = t_end;
            result->alignment_length[0] = t_start < 19 ? 19 - t_start : 0;
            result->alignment_length[1] = t_end >= 19 ? t_end-18 : 0;
            result->total_matches = matches[0] + matches[1];
            result->total_alignment_length = 1 + (t_end - t_start);
            result->total_identity = 100.0 * result->total_matches / result->total_alignment_length;
            result->identity[0] = 100.0 * result->matches[0] / result->alignment_length[0];
            result->identity[1] = 100.0 * result->matches[1] / result->alignment_length[1];
            strict_check(result);
        }
    }
}

/*----------------------------------------------------------------------*
 * Function:   output_alignment
 * Purpose:    
 * Parameters: 
 * Returns:    None
 *----------------------------------------------------------------------*/
void output_alignment(MPStats* stats, FastQRead* read, AlignmentResult* result)
{
    char s1[MAX_READ_LENGTH];
    char s2[MAX_READ_LENGTH];
    int i;

    fprintf(stats->log_fp, "\n    Id: %s\n", read->read_header);
    fprintf(stats->log_fp, "  Read: %s\n", read->read);    

    if (result->score < 0) {
        fprintf(stats->log_fp, "        No alignment\n");
    } else {
        for (i=0; i<read->read_size; i++) {
            int p = i - result->position;
            
            s1[i] = ' ';
            s2[i] = ' ';
            
            if ((p >= 0) && (p < strlen(transposon))) {
                if (transposon[p] == read->read[i]) {
                    s1[i] = '|';
                }
                s2[i] = transposon[p];
            }                
        }
   
        s1[read->read_size] = 0;
        s2[read->read_size] = 0;
        
        fprintf(stats->log_fp, "        %s\n", s1);
        fprintf(stats->log_fp, "        %s\n", s2);
    }
    
    fprintf(stats->log_fp, " Match: read base %d to %d, transposon base %d to %d\n", result->read_start, result->read_end, result->transposon_start, result->transposon_end);
    fprintf(stats->log_fp, 
            " Score: %d Id %.2f length %d (breakdown %d,%d matches %d,%d mismatches %d,%d lengths %.2f,%.2f identity)\n",
           result->score, 
           result->total_identity,
           result->total_alignment_length,
           result->matches[0], result->matches[1],
           result->mismatches[0], result->mismatches[1],
           result->alignment_length[0], result->alignment_length[1],
           result->identity[0], result->identity[1]);
    fprintf(stats->log_fp, "Result: %s\n", result->accepted == 1 ? "GOOD ALIGNMENT":"BAD");
}

/*----------------------------------------------------------------------*
 * Function:   get_read
 * Purpose:    
 * Parameters: 
 * Returns:    None
 *----------------------------------------------------------------------*/
int get_read(FILE* fp, FastQRead* read)
{
    int got_read = 1;
    
    
    if (!fgets(read->read_header, 1024, fp)) {
        got_read = 0;
    }

    if (!fgets(read->read, MAX_READ_LENGTH, fp)) {
        got_read = 0;
    }
    
    if (!fgets(read->quality_header, 1024, fp)) {
        got_read = 0;
    }
    
    if (!fgets(read->qualities, MAX_READ_LENGTH, fp)) {
        got_read = 0;
    }

    if (got_read == 1) {
        chomp(read->read_header);
        chomp(read->read);
        chomp(read->quality_header);
        chomp(read->qualities);
        read->read_size=strlen(read->read);
    }    
    
    return got_read;
}

/*----------------------------------------------------------------------*
 * Function:   write_read
 * Purpose:    
 * Parameters: 
 * Returns:    None
 *----------------------------------------------------------------------*/
void write_read(FastQRead* read, FILE *fp)
{
    fprintf(fp, "%s\n", read->read_header);
    fprintf(fp, "%s\n", read->read);
    fprintf(fp, "%s\n", read->quality_header);
    fprintf(fp, "%s\n", read->qualities);
}

/*----------------------------------------------------------------------*
 * Function:   decide_category
 * Purpose:    
 * Parameters: 
 * Returns:    None
 *----------------------------------------------------------------------*/
void decide_category(MPStats* stats, FastQRead* read_one, AlignmentResult* result_one, FastQRead* read_two, AlignmentResult* result_two)
{
    int category = -1;
    int trim_one = 0;
    int trim_two = 0;
    int l_one;
    int l_two;
    int smallest;
    
    if ((result_one->accepted == 1) && (result_two->accepted == 1)) {
        // Category A
        category = 0;
    } else if ((result_one->accepted == 0) && (result_two->accepted == 1)) {
        // Category B
        category = 1;
        // If using category E, check if we can get a relaxed hit on R1
        if (use_category_e == 1) {
            relaxed_check(result_one);
            if (result_one->accepted == 1) {
                read_one->read[result_one->read_start] = 0;
                read_one->qualities[result_one->read_start] = 0;
                stats->count_by_category_relaxed_hit[category]++;
                category=4;
            } else {
                trim_one = 1;
            }
        } else {
            trim_one = 1;
        }
    } else if ((result_one->accepted == 1) && (result_two->accepted == 0)) {
        // Category C
        category = 2;
        // If using category E, check if we can get a relaxed hit on R2
        if (use_category_e == 1) {
            relaxed_check(result_two);
            if (result_two->accepted == 1) {
                read_two->read[result_two->read_start] = 0;
                read_two->qualities[result_two->read_start] = 0;
                stats->count_by_category_relaxed_hit[category]++;
                category=4;
            } else {
                trim_two = 1;
            }
        } else {
            trim_two = 1;
        }
    } else {
        // Category D
        category = 3;
        trim_one = 1;
        trim_two = 1;
    }

    // Trim ends?
    if (trim_ends > 0) {
        if (trim_one) {
            int new_length = strlen(read_one->read) - trim_ends;
            read_one->read[new_length] = 0;
            read_one->qualities[new_length] = 0;
        }
        if (trim_two) {
            int new_length = strlen(read_two->read) - trim_ends;
            read_two->read[new_length] = 0;
            read_two->qualities[new_length] = 0;
        }
    }
    
    // Keep count
    stats->count_by_category[category]++;
    
    if (stats->log_fp != 0) {
        fprintf(stats->log_fp, "\n-------------------- Category %c --------------------\n", 'A' + category);   
    }
    
    l_one = strlen(read_one->read);
    l_two = strlen(read_two->read);
    smallest = l_one < l_two ? l_one:l_two;
    
    stats->read_length_counts[category][0][l_one]++;
    stats->read_length_counts[category][1][l_two]++;
    
    stats->read_pair_length_counts[category][smallest]++;
    
    if ((strlen(read_one->read) < (minimum_read_size-1)) ||
        (strlen(read_two->read) < (minimum_read_size-1))) {
        stats->count_by_category_too_short[category]++;
    } else {
        stats->count_by_category_long_enough[category]++;
        // Write reads
        write_read(read_one, stats->output_fp[category][0]);
        write_read(read_two, stats->output_fp[category][1]);
    }

}

/*----------------------------------------------------------------------*
 * Function:   process_file
 * Purpose:    
 * Parameters: 
 * Returns:    None
 *----------------------------------------------------------------------*/
void check_read_ids(MPStats* stats, FastQRead* read_one, FastQRead* read_two)
{
    int p=0;
    
    while ((p < strlen(read_one->read_header)) &&
           (p < strlen(read_two->read_header))) {
        if (read_one->read_header[p] != read_two->read_header[p]) {
            printf("Error: headers don't match up: %s and %s\n", read_one->read_header, read_two->read_header);
            exit(2);
        }
        
        if (read_one->read_header[p] <= ' ') {
            break;
        } else {
            p++;
        }
    }
    
    //if (stats->log_fp != 0) {
    //    char read_id[1024];
    //    strncpy(read_id, read_one->read_header+1, p-1);
    //    read_id[p-1] = 0;
    //    fprintf(stats->log_fp, "==================== %s ====================\n", read_id); 
    //}
}

/*----------------------------------------------------------------------*
 * Function:   valid_bases
 * Purpose:
 * Parameters:
 * Returns:    None
 *----------------------------------------------------------------------*/
boolean valid_bases(char *string)
{
    int i;
    
    for (i=0; i<strlen(string); i++) {
        if ((string[i] != 'A') && (string[i] != 'C') && (string[i] != 'G') && (string[i] != 'T')) {
            return false;
        }
    }
    
    return true;
}

/*----------------------------------------------------------------------*
 * Function:   check_valid_bases_and_gc_content
 * Purpose:
 * Parameters:
 * Returns:    None
 *----------------------------------------------------------------------*/
boolean check_valid_bases_and_gc_content(char *string, int* gc)
{
    int i; 

    *gc = 0;

    for (i=0; i<strlen(string); i++) {
        if ((string[i] == 'G') || (string[i] == 'C')) {
            *gc = (*gc) + 1;
        } else if  ((string[i] != 'A') && (string[i] != 'T')) {
            return false;
        }
    }

    return true;
}

#ifdef USE_MULTIPLE_HASHES
/*----------------------------------------------------------------------*
 * Function:   check_pcr_duplicates
 * Purpose:
 * Parameters:
 * Returns:    None
 *----------------------------------------------------------------------*/
boolean check_pcr_duplicates(FastQRead* read_one, FastQRead* read_two, MPStats* stats)
{
    char kmer_string[TOTAL_KMER_SIZE+1];
    BinaryKmer kmer;
    BinaryKmer tmp_kmer;
    Element* e;
    boolean is_duplicate = 0;
    boolean found = false;
    int i;

    if ((!valid_bases(read_one->read)) || (!valid_bases(read_two->read))) {
        stats->pairs_containing_n++;
        return false;
    }
    
    if (kmer_offsets[0][0] == -1) {
        kmer_offsets[0][0] = 0;
        kmer_offsets[0][1] = read_one->read_size / 3;
        kmer_offsets[1][0] = 0;
        kmer_offsets[1][1] = 2 * read_one->read_size / 3;
        kmer_offsets[2][0] = read_one->read_size / 3;
        kmer_offsets[2][1] = 2 * read_one->read_size / 3;
    }
    
    for (i=0; i<NUMBER_OF_HASHES; i++) {
        strncpy(kmer_string,                        read_one->read + kmer_offsets[i][0], SEPARATE_KMER_SIZE);
        strncpy(kmer_string+(1*SEPARATE_KMER_SIZE), read_one->read + kmer_offsets[i][1], SEPARATE_KMER_SIZE);
        strncpy(kmer_string+(2*SEPARATE_KMER_SIZE), read_two->read + kmer_offsets[i][0], SEPARATE_KMER_SIZE);
        strncpy(kmer_string+(3*SEPARATE_KMER_SIZE), read_two->read + kmer_offsets[i][1], SEPARATE_KMER_SIZE);
        kmer_string[TOTAL_KMER_SIZE]=0;

        seq_to_binary_kmer(kmer_string, TOTAL_KMER_SIZE, &kmer);
        Key key = element_get_key(&kmer, TOTAL_KMER_SIZE, &tmp_kmer);
        e = hash_table_find_or_insert(key, &found, kmer_hashes[i]);
        if (found) {
            is_duplicate = true;
            e->count++;
        } else {
            e->flags = 1;
            e->count = 1;
        }
    }
    
    if (is_duplicate) {
        stats->n_duplicates++;
        fprintf(stats->duplicates_fp, "Match: %s\n", kmer_string);
        fprintf(stats->duplicates_fp, "   R1: %s\n", read_one->read);
        fprintf(stats->duplicates_fp, "   R2: %s\n\n", read_two->read);
    }
    
    return is_duplicate;
}
#else
/*----------------------------------------------------------------------*
 * Function:   check_pcr_duplicates
 * Purpose:
 * Parameters:
 * Returns:    None
 *----------------------------------------------------------------------*/
boolean check_pcr_duplicates(FastQRead* read_one, FastQRead* read_two, MPStats* stats)
{
    char kmer_string[TOTAL_KMER_SIZE+1];
    BinaryKmer kmer;
    BinaryKmer tmp_kmer;
    Element* e;
    boolean is_duplicate = false;
    boolean found = false;
    int gc_one;
    int gc_two;
    
    if ((!check_valid_bases_and_gc_content(read_one->read, &gc_one)) || (!check_valid_bases_and_gc_content(read_two->read, &gc_two))) {
        stats->pairs_containing_n++;
        stats->n_invalid_for_duplicate++;
        return false;
    }

    stats->gc_bases+=gc_one;
    stats->gc_bases+=gc_two;
    stats->at_bases+=(read_one->read_size - gc_one);
    stats->at_bases+=(read_two->read_size - gc_two);
    
    gc_one=(gc_one*100)/read_one->read_size;
    gc_two=(gc_two*100)/read_two->read_size;
    
    stats->gc_content[0][gc_one]++;
    stats->gc_content[1][gc_two]++;

    strncpy(kmer_string, read_one->read, SEPARATE_KMER_SIZE);
    strncpy(kmer_string+(1*SEPARATE_KMER_SIZE), (read_one->read) + (read_one->read_size / 2), SEPARATE_KMER_SIZE);
    strncpy(kmer_string+(2*SEPARATE_KMER_SIZE), read_two->read, SEPARATE_KMER_SIZE);
    strncpy(kmer_string+(3*SEPARATE_KMER_SIZE), (read_two->read) + (read_two->read_size / 2), SEPARATE_KMER_SIZE);
    //strncpy(kmer_string, read_one->read, TOTAL_KMER_SIZE);
    kmer_string[TOTAL_KMER_SIZE]=0;
    
    seq_to_binary_kmer(kmer_string, TOTAL_KMER_SIZE, &kmer);
    Key key = element_get_key(&kmer, TOTAL_KMER_SIZE, &tmp_kmer);
    e = hash_table_find_or_insert(key, &found, duplicate_hash);
     
    if (found) {
        e->count++;
        is_duplicate = true;
            
        fprintf(stats->duplicates_fp, "Match: %s\n", kmer_string);
        fprintf(stats->duplicates_fp, "   R1: %s\n", read_one->read);
        fprintf(stats->duplicates_fp, "   R2: %s\n\n", read_two->read);
    } else {
        e->flags = 1;
        e->count = 1;
    }

    return is_duplicate;
}
#endif

/*----------------------------------------------------------------------*
 * Function:   process_file
 * Purpose:    
 * Parameters: 
 * Returns:    None
 *----------------------------------------------------------------------*/
void process_files(MPStats* stats)
{
    FastQRead reads[2];
    AlignmentResult alignments[2];
    int i, j;
    int is_duplicate;
    
    if (stats->log_filename[0] != 0) {
        char duplicates_filename[1024];

        stats->log_fp = fopen(stats->log_filename, "w");
        sprintf(duplicates_filename, "%s.pcr.txt", stats->log_filename);
        stats->duplicates_fp = fopen(duplicates_filename, "w");
    }
    
    
    // Open input files
    for (i=0; i<2; i++) {
        printf("Opening input filename %s\n", stats->input_filenames[i]);
        stats->input_fp[i] = fopen(stats->input_filenames[i], "r");
        if (!stats->input_fp[i]) {
            printf("Error: can't open file %s\n", stats->input_filenames[i]);
            exit(2);
        }
    }

    // Open output files
    for (i=0; i<num_categories; i++) {
        for (j=0; j<2; j++) {
            char filename[MAX_PATH_LENGTH];
            sprintf(filename, "%s_R%d.fastq", stats->output_filenames[i], j+1);
            printf("Opening output file %s\n", filename);
            stats->output_fp[i][j] = fopen(filename, "w");
            if (!stats->output_fp[i][j]) {
                printf("Error: can't open file %s\n", filename);
                exit(2);
            }
        }
    }
    
    // Read each entry in FASTQ files
    while (!feof(stats->input_fp[0])) {
        // Go next pair of reads
        int n_reads = 0;
        for (i=0; i<2; i++) {
            if (get_read(stats->input_fp[i], &reads[i]) == 1) {
                n_reads++;
            }

            if (stats->read_length == 0) {
                stats->read_length = strlen(reads[i].read);
            }
        }
        
        // Process pair
        if (n_reads == 2) {
            // Check read IDs match up
            check_read_ids(stats, &reads[0], &reads[1]);

            // Count pairs
            stats->num_read_pairs++;
            
            // Handle PCR duplicates
            is_duplicate = check_pcr_duplicates(&reads[0], &reads[1], stats);
            
            if ((remove_duplicates == 0) ||
                ((remove_duplicates == 1) && (is_duplicate == 0))) {
            
                for (i=0; i<2; i++) {
                    // Find transposons
                    find_transposons(&reads[i], &alignments[i]);

                    // Display log information
                    if (stats->log_fp != 0) {
                        output_alignment(stats, &reads[i], &alignments[i]);                        
                    }
                                    
                    // If adaptor found...
                    if (alignments[i].accepted == 1) {
                        // Count
                        stats->count_adaptor_found[i]++;
                        // Trim
                        reads[i].read[alignments[i].read_start] = 0;
                        reads[i].qualities[alignments[i].read_start] = 0;
                        // Check for too short
                        if (alignments[i].read_start < (minimum_read_size-1)) {
                            stats->count_too_short[i]++;
                        } else {
                            stats->count_long_enough[i]++;
                        }
                    } else {
                        stats->count_no_adaptor[i]++;
                    }
                }
                
                // Decide category (A, B, C, D)
                decide_category(stats, &reads[0], &alignments[0], &reads[1], &alignments[1]);
            } else {
                stats->duplicates_not_written++;
            }
        } else if (n_reads == 1) {
            printf("Error: Only managed to get one read!\n");
            exit(2);            
        }
    }
    
    // Close files
    for (i=0; i<2; i++) {        
        fclose(stats->input_fp[i]);
    }
    
    for (i=0; i<num_categories; i++) {
        for (j=0; j<2; j++) {
            fclose(stats->output_fp[i][j]);
        }
    }
    
    if (stats->log_fp != 0) {
        fprintf(stats->log_fp, "\nDONE\n");
        fclose(stats->log_fp);
    }
    
    if (stats->duplicates_fp != 0) {
        fclose(stats->duplicates_fp);
    }
}

/*----------------------------------------------------------------------*
 * Function:   process_adaptor
 * Purpose:    
 * Parameters: 
 * Returns:    None
 *----------------------------------------------------------------------*/
void process_adaptor(void)
{
    char reverse[1024];

    reverse_compliment(adaptor, reverse);
    strcpy(transposon, adaptor);
    strcat(transposon, reverse);

    printf("Adaptor: %s\n\n", transposon);
}

/*----------------------------------------------------------------------*
 * Function:   calculate_stats
 * Purpose:    
 * Parameters: 
 * Returns:    None
 *----------------------------------------------------------------------*/
void calculate_stats(MPStats* stats)
{
    int i;
    
    if (stats->num_read_pairs < 1) {
        printf("Error: Number of read pairs < 1!\n");
        exit(2);
    }
    
    stats->percent_pairs_containing_n = (100.0 * (double)stats->pairs_containing_n) / (double)stats->num_read_pairs;
        
    for (i=0; i<2; i++) {
        if (stats->count_adaptor_found[i] > 0) {
            stats->percent_adaptor_found[i] = (100.0 * (double)stats->count_adaptor_found[i]) / (double)stats->num_read_pairs;
        } else {
            stats->percent_adaptor_found[i] = 0;
        }
        
        if (stats->count_no_adaptor[i] > 0) {
            stats->percent_no_adaptor[i] = (100.0 * (double)stats->count_no_adaptor[i]) / (double)stats->num_read_pairs;
        } else {
            stats->percent_no_adaptor[i] = 0;
        }
        
        if (stats->count_long_enough[i] > 0) {
            stats->percent_long_enough[i] = (100.0 * (double)stats->count_long_enough[i]) / (double)stats->num_read_pairs;
        } else {
            stats->percent_long_enough[i] = 0;
        }
        
        if (stats->count_too_short[i] > 0) {            
            stats->percent_too_short[i] = (100.0 * (double)stats->count_too_short[i]) / (double)stats->num_read_pairs;
        } else {
            stats->percent_too_short[i] = 0;
        }
    }
    
    stats->total_too_short = 0;
    stats->total_long_enough = 0;   
 
    for (i=0; i<num_categories; i++) {
        if (stats->count_by_category[i] > 0) {
            stats->percent_by_category[i] = (100.0 * (double)stats->count_by_category[i]) / (double)stats->num_read_pairs;
        } else {
            stats->percent_by_category[i] = 0;
        }
        
        if (stats->count_by_category_long_enough[i] > 0) {
            stats->percent_by_category_long_enough[i] = (100.0 * (double)stats->count_by_category_long_enough[i]) / (double)stats->num_read_pairs;
        } else {
            stats->percent_by_category_long_enough[i] = 0;
        }
        
        if (stats->count_by_category_too_short[i] > 0) {
            stats->percent_by_category_too_short[i] = (100.0 * (double)stats->count_by_category_too_short[i]) / (double)stats->num_read_pairs;
        } else {
            stats->percent_by_category_too_short[i] = 0;
        }
        
        if (stats->count_by_category_relaxed_hit[i] > 0) {
            stats->percent_by_category_relaxed_hit[i] = (100.0 * (double)stats->count_by_category_relaxed_hit[i]) / (double)stats->num_read_pairs;            
        } else {
            stats->percent_by_category_relaxed_hit[i] = 0;
        }
        
        stats->total_too_short += stats->count_by_category_too_short[i];
        stats->total_long_enough += stats->count_by_category_long_enough[i];
    }
    
    // Total usable is categories A, B, C, E which are long enough
    stats->total_usable = stats->count_by_category_long_enough[0] + stats->count_by_category_long_enough[1] + stats->count_by_category_long_enough[2] + stats->count_by_category_long_enough[4];
    stats->percent_usable = (100.0 * (double)stats->total_usable) / (double)stats->num_read_pairs;
    
    stats->percent_total_too_short = (100.0 * (double)stats->total_too_short) / (double)stats->num_read_pairs;
    stats->percent_total_long_enough = (100.0 * (double)stats->total_long_enough) / (double)stats->num_read_pairs;
    
    if (stats->duplicates_not_written > 0) {
        stats->percent_duplicates_not_written = (100.0 * (double)stats->duplicates_not_written) / (double)stats->num_read_pairs;
    } else {
        stats->percent_duplicates_not_written = 0;
    }
    
    // GC content
    printf("GC bases: %ld  AT bases: %ld\n", stats->gc_bases, stats->at_bases);
    stats->percent_gc = (100.0 * stats->gc_bases) / (stats->gc_bases + stats->at_bases);
}

/*----------------------------------------------------------------------*
 * Function:   report_stats
 * Purpose:    
 * Parameters: 
 * Returns:    None
 *----------------------------------------------------------------------*/
void report_stats(MPStats* stats)
{
    int i;
    
    printf("\nSUMMARY\n\n");
    printf("     Strict match parameters: %d, %d\n", strict_double_match, strict_single_match);

    if (use_category_e == 1) {
        printf("    Relaxed match parameters: %d, %d\n", relaxed_double_match, relaxed_single_match);
    }
    
    printf("           Minimum read size: %d\n", minimum_read_size);
    printf("                   Trim ends: %d\n", trim_ends);
    printf("\n");
    printf("        Number of read pairs: %d\n", stats->num_read_pairs);
    printf("   Number of duplicate pairs: %d\t%.2f %%\n", stats->n_duplicates, stats->percent_duplicates);
    printf("Number of pairs containing N: %ld\t%.2f %%\n", stats->pairs_containing_n, stats->percent_pairs_containing_n);
    
    for (i=0; i<2; i++) {
        printf("   R%d Num reads with adaptor: %d\t%.2f %%\n", i+1, stats->count_adaptor_found[i], stats->percent_adaptor_found[i]);
        printf("       R%d long adaptor reads: %d\t%.2f %%\n", i+1, stats->count_long_enough[i], stats->percent_long_enough[i]);
        printf("          R%d reads too short: %d\t%.2f %%\n", i+1, stats->count_too_short[i], stats->percent_too_short[i]);
        printf("     R%d Num reads no adaptor: %d\t%.2f %%\n", i+1, stats->count_no_adaptor[i], stats->percent_no_adaptor[i]);
    }
        
    for (i=0; i<num_categories; i++) {
        printf("   Total pairs in category %c: %d\t%.2f %%\n", 'A'+i, stats->count_by_category[i], stats->percent_by_category[i]);
        printf("         %c pairs long enough: %d\t%.2f %%\n", 'A'+i, stats->count_by_category_long_enough[i], stats->percent_by_category_long_enough[i]);
        printf("           %c pairs too short: %d\t%.2f %%\n", 'A'+i, stats->count_by_category_too_short[i], stats->percent_by_category_too_short[i]);
    }
 
    printf("          Total usable pairs: %d\t%.2f %%\n", stats->total_usable, stats->percent_usable);
    printf("             All long enough: %d\t%.2f %%\n", stats->total_long_enough, stats->percent_total_long_enough);
    printf("    All categories too short: %d\t%.2f %%\n", stats->total_too_short, stats->percent_total_too_short);
    printf("      Duplicates not written: %d\t%.2f %%\n", stats->duplicates_not_written, stats->percent_duplicates_not_written);

    if (use_category_e == 1) {
        printf("         Category B became E: %d\t%.2f %%\n", stats->count_by_category_relaxed_hit[1], stats->percent_by_category_relaxed_hit[1]);
        printf("         Category C became E: %d\t%.2f %%\n", stats->count_by_category_relaxed_hit[2], stats->percent_by_category_relaxed_hit[2]);
    }

    printf("          Overall GC content: %.2f %%\n", stats->percent_gc);
    
    printf("\n");
}

/*----------------------------------------------------------------------*
 * Function:   report_stats
 * Purpose:    
 * Parameters: 
 * Returns:    None
 *----------------------------------------------------------------------*/
void output_histograms(MPStats* stats)
{
    char filename[MAX_PATH_LENGTH];
    FILE* fp;
    int c, r, i;
    
    for (c=0; c<num_categories; c++) {
        for (r=0; r<2; r++) {            
            sprintf(filename, "%s_%c_R%d_hist.txt", stats->output_prefix, 'A'+c, r+1);
            fp = fopen(filename, "w");
            if (fp) {
                for (i=1; i<=stats->read_length; i++) {
                    fprintf(fp, "%d\t%d\n", i, stats->read_length_counts[c][r][i]);
                }
                fclose(fp);
            } else {
                printf("Error: can't open histogram file %s\n", filename);
            }
        }
        
        sprintf(filename, "%s_%c_pair_hist.txt", stats->output_prefix, 'A'+c);
        fp = fopen(filename, "w");
        if (fp) {
            long int cumulative[stats->read_length+1];
            
            // Make cumulative totals
            for (i=stats->read_length; i>=1; i--) {
                if (i==stats->read_length) {
                    cumulative[i]=0;
                } else {
                    cumulative[i]=cumulative[i+1];
                }
                cumulative[i]+=stats->read_pair_length_counts[c][i];
            }
            
            for (i=1; i<=stats->read_length; i++) {
                fprintf(fp, "%d\t%d\t%ld\n", i, stats->read_pair_length_counts[c][i], cumulative[i]);
            }
            fclose(fp);
        } else {
            printf("Error: can't open histogram file %s\n", filename);
        }
    }

    for (r=0; r<2; r++) {
        sprintf(filename, "%s_R%d_gc.txt", stats->output_prefix, r+1);
        fp = fopen(filename, "w");
        if (fp) {
            for (i=0; i<=99; i+=2) {
                fprintf(fp, "%d\t%d\n", i, stats->gc_content[r][i]+stats->gc_content[r][i+1]);
            }
            fclose(fp);
        } else {
            printf("Error: can't open GC histogram file %s\n", filename);
        }
    }
}

/*----------------------------------------------------------------------*
 * Function:   calculate_pcr_duplicate_stats
 * Purpose:
 * Parameters:
 * Returns:    None
 *----------------------------------------------------------------------*/
#ifdef USE_MULTIPLE_HASHES
void calculate_pcr_duplicate_stats(MPStats* stats)
{
    int n_duplicates;
    double pc_duplicates;
    int t;

    void store_duplicates(Element * node) {
        if (node->count > 1) {
            n_duplicates+=(node->count-1);
        } 
    }

    printf("\n");
    for (t=0; t<3; t++) {
        printf("Counting duplicates for table %d\n", t);
        hash_table_print_stats(kmer_hashes[t]);
        n_duplicates=0;
        hash_table_traverse(&store_duplicates, kmer_hashes[t]);
        pc_duplicates=(double)((100.0*n_duplicates)/(double)stats->num_read_pairs);
        printf("\nDuplicates: %d (%.2f %%)\n", n_duplicates, pc_duplicates);
    }
 
    stats->percent_duplicates=(double)((100.0*stats->n_duplicates)/(double)stats->num_read_pairs);
}
#else
void calculate_pcr_duplicate_stats(MPStats* stats)
{
    FILE* fp;
    char filename[MAX_PATH_LENGTH];
    int duplicate_counts[MAX_DUPLICATES];
    int largest_count = 0;
    int i;
    
    void store_duplicates(Element * node) {
        int count = node->count;

        if (count > largest_count) {
            largest_count = count;
        }
        
        if (count >= MAX_DUPLICATES) {
            count = MAX_DUPLICATES - 1;
            printf("Warning: count (%d) exceeds maximum - treated as %d\n", count, MAX_DUPLICATES-1);
        }
        
        duplicate_counts[count]++;
        
        if (count > 1) {
            stats->n_duplicates+=(count-1);
        }
    }
    
    printf("\n");
    hash_table_print_stats(duplicate_hash);
    printf("\nCounting duplicates...\n");
    
    for (i=0; i<MAX_DUPLICATES; i++) {
        duplicate_counts[i] = 0;
    }
    
    stats->n_duplicates=0;
    
    hash_table_traverse(&store_duplicates, duplicate_hash);
    
    stats->percent_duplicates=(double)((100.0*stats->n_duplicates)/(double)stats->num_read_pairs);

    sprintf(filename, "%s_duplicates.txt", stats->output_prefix);
    
    fp = fopen(filename, "w");
    if (fp) {
        fprintf(fp, "n\tCount\tPercent\n");
        for (i=1; ((i<MAX_DUPLICATES) && (i<=largest_count)); i++) {
            int count = (i*duplicate_counts[i]);
            double percent = (100.0*count)/(double)stats->num_read_pairs;
            
            if (i == 1) {
                count += stats->n_invalid_for_duplicate;
            }
            
            fprintf(fp, "%d\t%d\t%.2f\n", i, count, percent);
        }
        fclose(fp);
    } else {
        printf("Error: can't open duplicates file %s\n", filename);
    }
}
#endif 

/*----------------------------------------------------------------------*
 * Function:   create_hash_table
 * Purpose:
 * Parameters:
 * Returns:    None
 *----------------------------------------------------------------------*/
void create_hash_table(void)
{
    double utilisation = 0.8;
    double required_entries = approximate_reads / utilisation;
    int b = 100;
    double c = log(required_entries/b)/log(2);
    int n = ceil(c);
    int entries = pow(2.0, (double)n)*b;
    long int memory = entries*sizeof(Element);
#ifdef USE_MULTIPLE_HASHES
    int i;

    printf("Creating hash tables for duplicate storage...\n");
    printf("                n: %d\n", n);
    printf("                b: %d\n", b);
    printf("          Entries: %d\n", entries);
    printf("       Entry size: %ld\n", sizeof(Element));
    printf(" Number of tables: %d\n", NUMBER_OF_HASHES);
    printf("  Memory required: %ld MB\n\n", (memory/(1024*1024))*NUMBER_OF_HASHES);
    
    for (i=0; i<NUMBER_OF_HASHES; i++) {
        kmer_hashes[i] = hash_table_new(n, b, 25, TOTAL_KMER_SIZE);
        hash_table_print_stats(kmer_hashes[i]);
        kmer_offsets[i][0] = -1;
        kmer_offsets[i][1] = -1;
    }
#else
    printf("Creating hash table for duplicate storage...\n");
    printf("                n: %d\n", n);
    printf("                b: %d\n", b);
    printf("          Entries: %d\n", entries);
    printf("       Entry size: %ld\n", sizeof(Element));
    printf("  Memory required: %ld MB\n\n", memory/(1024*1024));
    
    duplicate_hash = hash_table_new(n, b, 25, TOTAL_KMER_SIZE);
    hash_table_print_stats(duplicate_hash);
#endif
}


/*----------------------------------------------------------------------*
 * Function:   main
 *----------------------------------------------------------------------*/
int main(int argc, char* argv[])
{
    MPStats stats;
    
    printf("\nNextClip v%s\n\n", NEXTCLIP_VERSION);
    
    initialise_stats(&stats);
    parse_command_line(argc, argv, &stats);
    create_hash_table();

    process_adaptor();
    process_files(&stats);
    calculate_stats(&stats);
    calculate_pcr_duplicate_stats(&stats);
    report_stats(&stats);
    output_histograms(&stats);
    
    printf("\nDone.\n");
    
    return 0;
}
