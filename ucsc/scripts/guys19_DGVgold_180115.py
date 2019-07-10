#!/usr/bin/env python2
'''
guys19_DGVgold_180115.py

Usage:
guys19_DGVgold_180115.py -i <input_file> [ -o <output_prefix> ]

Takes the DGV Gold Standard Variant input file (.GFF3) and returns two BED detail files for use as a
UCSC annotation track. One output file contains all DGV Gold 'Gain' variants, while the second
contains DGV Gold 'Loss' variants. Default output prefix is guys19_DGVgold_YYMMDD_gain/loss.bed

This script was tested with the DGV Gold Standard Variant file for GRch37(hg19) [2016-05-15],
available at: http://dgv.tcag.ca/dgv/docs/DGV.GS.March2016.50percent.GainLossSep.Final.hg19.gff3
md5sum : b2e8022d1073ceb33f3250d5ea87ebb1  DGV.GS.March2016.50percent.GainLossSep.Final.hg19.gff3

Author : Nana Mensah
Created : 15th Jan 2018
(Updated: Andy Bond, 5th April 2018)
'''
import re
import argparse
import datetime


# Create command-line argument parser for script input and outputs
parser = argparse.ArgumentParser()
# Set input file
parser.add_argument('-i', metavar='input_file', help='Input DGV Gold GFF3 file', required=True)
# Set output file prefix. This argument has a default attribute set: 'guys19_DGVgold_YYMMDD', where
# YYMMDD is formatted to the present date using datetime.today().
parser.add_argument('-o', metavar='output_prefix', help='Prefix for output *_gain.bed and' +
                    '*_loss.bed files. Default: guys19_DGVgold_YYMMDD',
                    default='guys19_DGVgold_' + datetime.datetime.today().strftime('%y%m%d'))
# Read command line arguments into the object args for use
args = parser.parse_args()

# Compile a regular expression for capturing track data from the final column (#9) of the input file:
regex = re.compile(r'ID=(?P<ID>\w+);.*' +  # Get the ID
                   r'variant_sub_type=(?P<subtype>\w+);.*' +  # Get the subtype (Gain/Loss)
                   r'inner_start=(?P<istart>\d+);' +  # Get the common start position of the track
                   r'inner_end=(?P<iend>\d+);.*' +  # Get the common stop position of the track
                   r'num_samples=(?P<positive_samples>\d+);.*' + # Number of samples positive for CNV
                   r'Frequency=(?P<frequency>.+?)%;.*' +  # Get the frequency (float/percentage)
                   r'num_unique_samples_tested=(?P<total_samples>\d+)') # Number of samples tested

# Set output file streams using command line argument (-o) or default
gain_out = open(args.o + '_gain.bed', 'w')
loss_out = open(args.o + '_loss.bed', 'w')
# Write UCSC BED detail header to appropriate output file stream.
gain_out.write('browser position chr9:140513444-140678845\n' +
               'track name="DGV Gold Gain" type=bedDetail db=hg19 visibility=3 itemRgb=On\n')
loss_out.write('browser position chr9:140513444-140678845\n' +
               'track name="DGV Gold Loss" type=bedDetail db=hg19 visibility=3 itemRgb=On\n')

# Intialise a set for storing DGV Gold SV IDs. Each ID is seen in triplicate in the GFF3 input file,
# as the inner- and outer-most probe locations are represented as 1bp calls. The output BED detail
# file requires the inner_start and inner_end annotations from any triplicates. This list is used to
# store previously written IDs and skip them in subsequent calls to process().
ID_list = set()


def bed_html(frequency, positive_samples, total_samples, ID):
    '''Return the HTML string for use in the last column of the BED output file. The arguments
    frequency, total_samples and ID are written into the html string using str.format()'''
    return ("<br /><h3>Frequency {}% ({} out of {} unique samples tested)</h3>" +
            "<a href=http://dgv.tcag.ca/gb2/gbrowse_details/dgv2_hg19?name={ID}>" +
            "DGV Gold call ID {ID}</a> (click to go to DGV)<br /><br />" +
            "DGV Gold Standard Variants Release Date 2016-05-15.<br />").format(frequency, positive_samples, total_samples, ID=ID)


def parse_data(line):
    '''Parse a line from the input GFF3 file by extracting tab-separated fields and reading the
    annotation field using the compiled regular expression. Returns a dictionary of the parsed data'''

    # Split line on tabs. Creates a list of 9 strings in the following format:
    # Chromosome, CNV, Text, Start,  Stop, '.', '.', '.', Annotation
    data = line.split('\t')
    # Set chromsome name from first entry of line to be used downstream in bed_string
    chromosome = data[0]
    # Extract data from the [Annotation] field using the compiled regular expression.
    # Produces a regular expression match object, where the group() method can be used to access
    # the named capture groups. E.g. annot.group('ID').
    annot = regex.search(data[8])

    # Set output file name and RGB colour based on DGV SV subtype
    if annot.group('subtype') == 'Gain':
        outfile = gain_out
        RGB = "0,0,255"  # Blue for DGV Gold Gain
    elif annot.group('subtype') == 'Loss':
        outfile = loss_out
        RGB = "255,0,0"  # Red for DGV Gold Loss
    # Raise error if subtype other than 'Gain' or 'Loss' found
    else:
        raise ValueError('Unkown subtype {} in GFF3 annotation'.format(annot.group('subtype')))

    return {'chromosome': chromosome, 'ID': annot.group('ID'), 'start': annot.group('istart'),
            'end': annot.group('iend'), 'frequency': annot.group('frequency'), 'outfile': outfile,
            'positive_samples': annot.group('positive_samples'), 'total_samples': annot.group('total_samples'), 
            'RGB': RGB}


def process(line):
    # Get dictionary of cleaned data from the input line
    data = parse_data(line)

    # If the ID for this DGV Gold call has been processed before, return None. When process() is
    # called in main(), this behaviour skips to the next line if the condition is met.
    if data['ID'] in ID_list:
        return None
    # Else, add the ID to ID_list and continue
    ID_list.add(data['ID'])

    # Set output file HTML string
    html = bed_html(data['frequency'], data['positive_samples'], data['total_samples'], data['ID'])

    # Create a list of the BED detail fields for this entry using the data dictionary.
    # Format: [chromosome], [start], [stop], [frequency(num_unique_samples_tested)], 0, .,
    #         [start], [stop], [RGB], [frequency(num_unique_samples_tested)], [html]
    bed_list = [data['chromosome'], data['start'], data['end'], data['frequency'] +
                '%(' + data['positive_samples'] + '/' + data['total_samples'] + ')', '0', '.', data['start'], data['end'], data['RGB'],
                data['frequency'] + '%(' + data['positive_samples'] + '/' + data['total_samples'] + ')', html]
    # Create BED detail string by joining fields in bed_list with tab delimiters
    bed_string = '\t'.join(bed_list)

    # Append the bed_string to the destination output file stream
    data['outfile'].write(bed_string + "\n")


def main():
    '''Convert input DGV Gold GFF3 file to BED detail format.'''

    # Open input file stream
    with open(args.i) as infile:
        # Loop through lines of input file and write bed detail file for UCSC tracks
        for line in infile:
            process(line)

    # Close output file streams
    gain_out.close()
    loss_out.close()


if __name__ == '__main__':
    main()
