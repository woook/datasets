#!/usr/env/bin python3
'''
Usage: guys_ClinGenHaplo.py ClinGen_haploinsufficiency_gene.bed

Convert ClinGen Haploinsufficiency BED file to BED-detail format for use in UCSC
annotation track. Outputs the file 'guys_ClinGenHaplo_YYMMDD.bed', where YYMMDD is formatted with
the date the script was called.

The input file (ClinGen_haploinsufficiency_gene.bed) was downloaded from:
ftp://ftp.ncbi.nlm.nih.gov/pub/dbVar/clingen/ClinGen_haploinsufficiency_gene_GRCh37.bed

Author: Nana Mensah
Created: 12th Feb 2018
'''
import sys
import datetime

# Create dictionary using the haploinsufficiency dosage scores as keys. Each key should contain a
# nested dictionary for the RGB color code (for BED-detail HTML) and ClinGen description for
# each score, where:
# 0 = No evidence available (Grey)
# 1 = Little evidence (Pink Light),
# 2 = Some evidence (Pink Dark),
# 3 = Sufficient evidience (Red),
# 30 = Associated with autosomal recessive phenotype (Green Light),
# 40 = Dosage sensitivity unlikely (Black)
RGB = {0: {"color": "220,220,220", "description": "No evidence available"},
       1: {"color": "255,210,210", "description": "Little evidence for dosage pathogenicity"},
       2: {"color": "255,155,155", "description": "Some evidence for dosage pathogenicity"},
       3: {"color": "255,0,0", "description": "Sufficient evidence for dosage pathogenicity"},
       30: {"color": "170,255,170", "description": "Gene associated with autosomal recessive phenotype"},
       40: {"color": "0,0,0", "description": "Dosage sensitivity unlikely"}}


def html(ID, score):
    '''Returns HTML string for annotation field of output BED-detail file. The arguments ID and
    score correspond to the Gene ID and dosage score parsed from a line of the input file. These
    values are written into the HTML string using str.format().'''
    return (("<br /><h2>{ID} - {description} (score:{score})</h2>" +
             "<a href=https://www.ncbi.nlm.nih.gov/projects/dbvar/clingen/clingen_gene.cgi?sym={ID}>" +
             "<strong>Click here to go to the ClinGen evidence page for this gene</strong></a><br /><br />")
            .format(ID=ID, description=RGB[score]["description"], score=score))


def parse(line):
    '''Return the corresponding BED-detail string for a line in the input ClinGen file'''
    # Split input line (tab-delimited) into a list containing five strings
    split_line = line.split("\t")
    # Print error and exit if incoming BED line does not split into 5 fields
    if len(split_line) != 5:
        print("Error: BED entry incorrect length")
        print(line)
        exit()
    # Store each item of the split line in a separate variable
    chrom, start, stop, ID, score = split_line
    # Clean the newline character from the final entry (score) and convert to an integer
    score_cln = int(score.strip())
    # Return the new BED-detail string by merging a list of field strings using tab characters ("\t")
    return "\t".join([chrom, start, stop, ID, "0", ".", start, stop, RGB[score_cln]["color"], ID,
                     html(ID, score_cln)])


def main():
    '''Convert input ClinGen Haploinsufficiency bed file to BED-detail format.'''
    # Set today's date in format YYMMDD for output file name
    datestring = datetime.datetime.today().strftime('%y%m%d')
    # Store output file name in a string. Output file is written to the current directory.
    outfile = "guys_ClinGenHaplo_" + datestring + ".bed"
    # Open output file stream and write hard-coded output file header for UCSC track.
    stdout = open(outfile, "w")
    stdout.write('browser position chr9:140513444-140678845\n' +
                 'track name="ClinGen Haploinsufficiency" type=bedDetail db=hg19 visibility=3 itemRgb=On\n')
    # Open ClinGen haploinsufficiency input file for reading. This is provided as the first command
    # line argument to the script.
    with open(sys.argv[1], 'r') as f:
        # Store all lines of the input file, except the header (line 1), in a list.
        line_list = f.readlines()[1:]
        # Loop over input file lines and pass lines to parse(). Write the returned
        # BED-detail line to the output file.
        for line in line_list:
            stdout.write(parse(line) + "\n")
    # Close the output file
    stdout.close()


if __name__ == '__main__':
    main()
