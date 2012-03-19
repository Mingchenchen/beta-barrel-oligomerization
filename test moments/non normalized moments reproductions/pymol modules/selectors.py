"""percent_sasa_and_ss(minsasa,targetss)
percent_sasa_and_ss(.2,'E') and percent_sasa_and_ss(.2,'[^E]?$') and
percent_sasa_and_ss(.2,'.*') are what I usually use"""

import functools
import re

def sel_by_percent_sasa_and_ss_generator(spreadsheet, minsasa, targetss):
    for resi, sasa,ss in spreadsheet.get_columns(['resi','%sasa','ss']):
        if sasa != '' and float(sasa) >= float(minsasa) and re.match(targetss,ss) is not None:
            yield 'i. ' + resi

def percent_sasa_and_ss(spreadsheet, minsasa, targetss):
    """percent_sasa_and_ss(minsasa,targetss)
    percent_sasa_and_ss(.2,'E') and percent_sasa_and_ss(.2,'[^E]?$') are
    what I usually use
    Returns a selection string"""
    return '(' + '|'.join(sel_by_percent_sasa_and_ss_generator(spreadsheet, minsasa, targetss)) + ')'