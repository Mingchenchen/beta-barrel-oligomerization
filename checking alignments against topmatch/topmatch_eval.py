from __future__ import division
import urllib
import json
from HTMLParser import HTMLParser

topmatch_address = 'http://topmatch.services.came.sbg.ac.at/json-topmatch/'
def topmatch_comparison(query, target):
    '''Given two PDB ID's, compare the structure on topmatch.
    Returns a'''
    
    # TopMatch only understands lower-case PDB ID's, so decapitalize them:
    query = query.lower()
    target = target.lower()

    # Retrieve the JSON-formatted results
    params = urllib.urlencode({'query': query, 'target': target})
    result_file = urllib.urlopen('{}?{}'.format(topmatch_address, params))
    
    # Parse the results
    result_dict = json.loads(result_file.read())

    return result_dict

class SeqReporter(HTMLParser):
    def __init__(self, instrumentation = True):
        self.instrumentation = instrumentation 

        self.inseq = False
        self.inseq_span_depth = 0
        self.accumulated_data = ''
        HTMLParser.__init__(self)
    
    def handle_starttag(self, tag, attrs):
        # If in a sequence, keep track of how many new span tags have
        # opened and closed
        if self.inseq and tag == 'span':
            self.inseq_span_depth += 1
            if self.instrumentation:
                print('down a level...')

        # If not previously within a sequence, determine whether we're
        # in one now
        if not self.inseq and tag == 'span':
            attr_dict = dict(attrs)
            # There'll be an exception here if there's a span tag with
            # no "class" attribute. I don't know if those exist though.
            if 'alig_seq' in attr_dict['class']:
                self.inseq = True
                if self.instrumentation:
                    print('in a sequence!')

    def handle_endtag(self, tag):
        # If in a sequence within additional span tags, keep track of
        # how many span tags have opened and closed
        if self.inseq and tag == 'span' and self.inseq_span_depth > 0:
            self.inseq_span_depth -= 1
            if self.instrumentation:
                print('up a level.')

        # If in a sequence in which no span tags are left open,
        # and this is a span tag, then mark that we've left the sequence
        elif self.inseq and tag == 'span' and self.inseq_span_depth == 0:
            self.inseq = False
            if self.instrumentation:
                print('not in a sequence.')
    
    def handle_data(self, data):
        # The only data I'm interested in is sequence data. If we're in
        # a sequence, save the data.
        if self.inseq:
            self.accumulated_data += data
            if self.instrumentation:
                print('found data: ' + data)

# testing
box = SeqReporter()
box.feed(y)
