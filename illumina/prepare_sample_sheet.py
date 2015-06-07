#-*- coding: utf-8 -*-

'''
Prepare sample sheet file for Illumina MiSEq sequencers.

'''

__version__ = "0.1"


import os
import click
import sys
from base import get_casava_sample_sheet


def log(category, message, *args, **kwargs):
    click.echo('%s: %s' % (click.style(category.ljust(10), fg='cyan'),
        message.replace('{}', click.style('{}', fg='yellow')).format(*args, **kwargs)))


@click.command(context_settings=dict(
               help_option_names=['-h', '--help'],
               ignore_unknown_options=True,))
@click.option('--view', '-v',
                is_flag=True,
                help='view contents of the SampleSheet')
@click.option('--fix-spaces', '-f',
                is_flag=True,
                help='replacce spaces in Sample ID and SampleProject fields with underscores')
@click.option('--fix-empty-projects', '-fep',
                is_flag=True,
                help='create SampleProject names where these are blank in the original')
@click.option('--ignore-warnings', '-w',
                is_flag=True,
               help='ignore warnings about spaces and duplicate SampleID/SampleProject')
@click.option('--truncate-barcodes', '-t',
                type=int,
                default = 0,
               help='trim barcode sequences in sample sheet to number of bases specified.'
                        'Default is to leave barcode sequence unaltered.')
@click.version_option(__version__)
@click.argument('samplesheet', type=click.Path(exists=True))
@click.argument('output',
                type=click.File(mode='w'))
@click.argument('prepare_ss_args', nargs=-1, type=click.UNPROCESSED)
def prepare_sample_sheet(samplesheet, output, view, fix_spaces, fix_empty_projects,
                           ignore_warnings, truncate_barcodes, prepare_ss_args):

    if not os.path.isfile(samplesheet):
        log("Error", "sample sheet not found:" , samplesheet)
        sys.exit(1)

    #Read the data as CSV
    data = get_casava_sample_sheet(samplesheet)
    print data

    '''
    #Truncate barcodes
    if truncate_barcodes:
        barcode_length = truncate_barcodes


    # Truncate barcodes
    if options.barcode_len is not None:
        barcode_len = options.barcode_len
        for line in data:
            barcode = truncate_barcode(line['Index'],options.barcode_len)
            print "Lane %d '%s/%s': barcode '%s' -> '%s'" \
                % (line['Lane'],
                   line['SampleProject'],
                   line['SampleID'],
                   line['Index'],
                   barcode)
            line['Index'] = barcode
    # Fix spaces
    if options.fix_spaces:
        data.fix_illegal_names()
    # Fix empty projects
    if options.fix_empty_projects:
        for line in data:
            if not line['SampleProject']:
                line['SampleProject'] = line['SampleID']
    # Fix duplicates
    if options.fix_duplicates:
        data.fix_duplicated_names()
    # Print transposed data in tab-delimited format
    if options.view:
        data.transpose().write(fp=sys.stdout,delimiter='\t')
    # Check for non-unique id/project combinations, spaces and empty names
    check_status = 0
    # Duplicated names
    duplicates = data.duplicated_names
    if len(duplicates) > 0:
        check_status = 1
        for duplicate_set in duplicates:
            for lane in duplicate_set:
                logging.warning("Duplicated SampleID/SampleProject in lane %s (%s/%s)" %
                                (lane['Lane'],lane['SampleID'],lane['SampleProject']))
    # Illegal characters/spaces in names
    illegal_names = data.illegal_names
    if len(illegal_names) > 0:
        check_status = 1
        for lane in illegal_names:
            logging.warning("Spaces in SampleID/SampleProject in lane %s (%s/%s)" %
                            (lane['Lane'],lane['SampleID'],lane['SampleProject']))
    # Empty names
    empty_names = data.empty_names
    if len(empty_names) > 0:
        check_status = 1
        for lane in empty_names:
            logging.warning("Empty SampleID and/or SampleProject name in lane %s (%s/%s)" %
                            (lane['Lane'],lane['SampleID'],lane['SampleProject']))
    # Predict outputs
    if check_status == 0 or options.ignore_warnings:
        projects = data.predict_output()
        print "Predicted output:"
        for project in projects:
            print "%s (%d samples)" % (project,len(projects[project]))
            for sample in projects[project]:
                print "\t%s" % sample
                for sub_sample in projects[project][sample]:
                    print "\t\t%s" % sub_sample

    # Write out new sample sheet
    if options.samplesheet_out:
        if check_status and not options.ignore_warnings:
            logging.error("please fix above errors in sample sheet data")
        else:
            data.write(options.samplesheet_out)
    # Finish
    sys.exit(check_status)
    '''


if __name__ == '__main__':
    prepare_sample_sheet()

