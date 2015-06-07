"""
Provides classes for extracting data about runs of Illumina-based sequencers,
in our case MiSEQ from directory structure, data files and naming
conventions.

"""

__version__ =   '0.1'

from collections import OrderedDict

import csv
import re
import logging
import string
from utils import extract_initials

class IlluminaError(Exception):
    """Base class for errors with Illumina code"""

class IEMSampleSheet(object):
    '''
    Class tha encapsulates the IEM Experimental Manager format Sample sheet.
    '''

    def __init__(self, sample_sheet=None):
        self._headers = OrderedDict()
        self._reads = list()
        self._settings = OrderedDict()
        self._content = []
        self._content_headers = []
        self._data = []

        if sample_sheet is not None:
            ss = open(sample_sheet, 'rU')
            self._load_data(ss)

    def set_header(self, key, value):
        self._headers[key] = value

    def add_reads(self, key):
        self._reads.append(key)

    def set_settings(self, key, value):
        self._settings[key] = value

    def _load_data(self, ss):
        '''
        populate the data from external file
        '''
        section = None
        reader = csv.reader(ss)
        for idx, line in enumerate(reader):
            if line:  #skip blank lines
                if line[0].startswith('['):
                    # new section
                    try:
                        end = line[0].index(']')
                        section = line[0][1:end]
                        continue
                    except ValueError:
                        logging.error('Bad line %d: %s' % (idx+1, line[0]))
                if section == 'Header':
                    #Header lines
                    self.set_header(*line)
                elif section == 'Reads':
                    #Read lines
                    self.add_reads(*line)
                elif section == 'Data':
                    #Store the data
                    if not self._content_headers:
                        #Store the header
                        self._content_headers = line
                    else:
                        self._content.append(line)
                elif section == 'Settings':
                    # Settings lines are comma-separated PARAM,VALUE lines
                    self.set_settings(*line)
                elif section is None:
                    raise IlluminaError("Not a valid IEM sample sheet")
                else:
                    raise IlluminaError(
                        "Unrecognised section '%s': not a valid IEM sample sheet" % section)
        #Clean up data items: remove surrounding whitespace
        if self._content:
           self._data = [ OrderedDict(zip(self._content_headers, \
                            map(string.strip, sample))) for sample in self._content]

    @property
    def header_items(self):
        """Return list of items listed in the  [Header] section
        """
        return self._headers.keys()

    @property
    def header(self):
        """Return ordered dictionary for the [Header] section
        """
        return self._headers

    @property
    def reads(self):
        """Return list of values from the [Reads] section
        """
        return self._reads

    @property
    def settings_items(self):
        """Return list of items listed in the [Settings] section
        """
        return self._settings.keys()

    @property
    def settings(self):
        """Return ordered dictionary for the [Settings] section
        """
        return self._settings

    @property
    def samples(self):
        """Return Samples for the [Data] section
        """
        return self._data

    def show(self):
        """Reconstructed version of original sample sheet
        """
        s = []
        s.append('[Header]')
        for param in self._headers:
            s.append('%s,%s' % (param,self._headers[param]))
        s.append('')

        s.append('[Reads]')
        for value in self._reads:
            s.append(value)

        s.append('')
        s.append('[Settings]')

        for param in self._settings:
            s.append('%s,%s' % (param,self._settings[param]))

        s.append('')
        s.append('[Data]')

        s.append(','.join(self._content_headers))
        for line in self._content:
            s.append(','.join((line)))
        return '\n'.join(s)

    def casava_samplesheet(self, FCID='FC1', fix_empty_projects=False):
        """Return data as a CASAVA formatted sample sheet
        """
        samplesheet = CasavaSampleSheet()
        for line in self._data:
            #Set the lane
            try:
                lane = line['Lane']
            except KeyError:
                "No lane column"
                lane = 1

            #set the index tag
            try:
                index_tag = '%s-%s' % (line['index'].strip(),
                                       line['index2'].strip())
            except KeyError:
                #Assume no index 2
                try:
                    index_tag = line['index'].strip()
                except KeyError:
                    #No index
                    index_tag = ''

            sample = OrderedDict()
            sample['FCID'] = FCID
            sample['Lane'] = lane
            sample['Index'] = index_tag
            sample['SampleID'] = line['Sample_ID']
            sample['Description'] = line['Description']

            samplesheet.append(sample)

            # Fix project name
            if line['Sample_Project'] == '' and fix_empty_projects:
                # No project name - try to use initials from sample name.
                sample['SampleProject'] =  extract_initials(line['Sample_ID'])
            else:
                sample['SampleProject'] = line['Sample_Project']

        return samplesheet

def get_casava_sample_sheet(samplesheet,FCID_default='FC1'):
    '''
     It reads the CSV samplesheet from Illuminat platform and populates and
     returns a CasavaSampleSheet object with can be used to make a Casava
     Samplesheet format suitable for bcl2fastq conversion.
    '''
    #Open the file for reading.
    iem_ss = IEMSampleSheet(samplesheet)

    return iem_ss.casava_samplesheet()


class CasavaSampleSheet(object):
    '''
    Class tha encapsulates the CASAVA format Sample sheet.
    '''
    def __init__(self, samplesheet=None):
        '''
        Create a new CasavaSampleSheet instance
        '''
        self._content = []
        self._content_headers = ('FCID','Lane','SampleID','SampleRef',
                                 'Index','Description','Control',
                                     'Recipe','Operator','SampleProject')
        self._data = []

        self.illegal_characters = "?()[]/\=+<>:;\"',*^|&. \t"

        if samplesheet is not None:
            ss = open(samplesheet, 'rU')
            self._load_data(ss)

    def append(self, line):
        self._data.append(line)

    def __repr__(self):
        return '\n'.join([str(x) for x in self._data])

    def _load_data(self, ss):
        '''
        populate the data from external file
        '''
        headers = None
        reader = csv.reader(ss)
        for idx, line in enumerate(reader):
            if line:  #skip blank lines
                #Store the data
                if headers is None:
                    #Store the header
                    headers = line
                else:
                    self._content.append([el.strip('"') for el in line
                                            if not el.strip('"').startswith('#')])

        if headers is None:
            raise IlluminaError("Not a valid IEM sample sheet")
        elif tuple(headers) == self._content_headers:
            raise IlluminaError(
                "Unrecognised header '%s': not a valid IEM sample sheet" % headers)

        #Clean up data items: remove surrounding whitespace
        if self._content:
           self._data = [OrderedDict(zip(self._content_headers, sample)) for sample in self._content]

    @property
    def duplicated_names(self):
        """List lines where the SampleID/SampleProject pairs are identical
        """
        samples = {}
        for line in self._data:
            name = ((line['SampleID'],line['SampleProject'],line['Index'],line['Lane']))
            if name not in samples:
                samples[name] = [line]
            else:
                samples[name].append(line)
        duplicates = []
        for name in samples:
            if len(samples[name]) > 1: duplicates.append(samples[name])
        return duplicates

    @property
    def empty_names(self):
        """List lines with blank SampleID or SampleProject names
        """
        empty_names = []
        for line in self._data:
            if str(line['SampleID']).strip() == ''  or str(line['SampleProject']).strip() == '':
                empty_names.append(line)
        return empty_names


    @property
    def illegal_names(self):
        """List lines with illegal characters in SampleID or SampleProject names
        """
        illegal_names = []
        for line in self._data:
            for c in self.illegal_characters:
                illegal = (str(line['SampleID']).count(c) > 0) \
                          or (str(line['SampleProject']).count(c) > 0)
                if illegal:
                    illegal_names.append(line)
                    break
        return illegal_names

