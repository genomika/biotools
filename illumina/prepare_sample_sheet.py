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


def truncate_barcode(sequence, length):
    """Return barcode sequence truncated to requested length.
    """
    try:
        index5, index7 = sequence.split('-')
        #dual index
        if len(index5) >= length:
            newindex5 = index5[:length]
        if len(index7) >= length:
            newindex7 = index7[:length]

        return '-'.join([newindex5, newindex7])
    except ValueError:
        # No hyphen: single index barcode
        return sequence[:length]

def verify_spreadsheet(ss):
    '''
    Check if the spreadsheet has invalid format
    '''
    log('info', 'Checking for any invalid parameters in samplesheet')
    #Check for duplicated names
    status = 0
    duplicated_names = ss.duplicated_names
    if len(duplicated_names) > 0:
        status = 1
        for duplicate in duplicates:
            log('warning', 'Duplicated SampleID/SampleProject in lane %s (%s/%s)'
                            % (lane['Lane'], lane['SampleID'], lane['SampleProject']))

    # Check for illegal names
    illegal_names = ss.illegal_names
    if len(illegal_names) > 0:
        status = 1
        for line in illegal_names:
            log('warning', 'Spaces in SampleID/SampleProject in lane %s (%s/%s)'
                                % (lane['Lane'], lane['SampleID'], lane['SampleProject']))

    #Check for Empty names
    empty_names = ss.empty_names
    if len(empty_names) > 0:
        status = 1
        for line in empty_names:
            log('warning', "Empty SampleID and/or SampleProject name in lane %s (%s/%s)" %
                            (line['Lane'],line['SampleID'],line['SampleProject']))

    return status


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

    #Truncate barcodes
    if truncate_barcodes:
        barcode_length = truncate_barcodes
        log('Info', 'Truncating the barcodes to length: %s' % barcode_length)
        for line in data:
            barcode = truncate_barcode(line['Index'], barcode_length)
            log('Info', '\tLane %d "%s/%s": barcode "%s" -> "%s"' %
                            (line['Lane'],
                             line['SampleProject'],
                             line['SampleID'],
                             line['Index'],
                             barcode))
            line['Index'] = barcode

    #Fix Spaces
    if fix_spaces:
        log('Info', 'Fixing illegal characters in SampleID and SampleProject pairs')
        data.fix_illegal_names()

    if fix_empty_projects:
        log('Info', 'Fixing missing SampleProject in samples')
        for line in data:
            if not line['SampleProject']:
                line['SampleProject'] = line['SampleID']
    if view:
        log('info', 'viewing data')
        print data

    invalid = verify_spreadsheet(data)

    #Check how it will be the outputs
    if not invalid or ignore_warnings:
        projects = data.expected_output()
        print 'Output:'
        for project in projects:
            print '%s (%d samples)' % (project, len(projects[project]))
            for sample in projects[project]:
                print '\t%s' % sample
                for subsample in projects[project][sample]:
                    print "\t\t%s" % subsample

    if invalid and not ignore_warnings:
        log('error', 'Please fix above errors in sample sheet')
    else:
        data.write(output)

    sys.exit(invalid)

if __name__ == '__main__':
    prepare_sample_sheet()

