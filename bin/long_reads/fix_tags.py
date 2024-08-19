#!/usr/bin/python

import sys
import pysam


##


def main():

    path_input = sys.argv[1]
    path_output = 'CB_UB.bam'

    input_bam = pysam.AlignmentFile(path_input, 'rb')
    output_bam = pysam.AlignmentFile(path_output, 'wb', template=input_bam)

    for read in input_bam:
        read_name = read.qname
        CB, UB = read_name[:29].split('_')
        read.qname = read_name[30:]
        read = read.set_tag('CB', CB)
        read = read.set_tag('UB', UB)
        output_bam.write(read)

    input_bam.close()
    output_bam.close()


##


# Run
if __name__ == '__main__':
    main()