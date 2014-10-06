'''
xmfa | xmfa.py
~~~~~~~~~~~~~~

Lightweight XMFA parser. No dependencies. Files are parsed into simple 
namedtuples as documented below.

'''
from __future__ import print_function

import re
from collections import namedtuple


_sub_aligment_group_re = re.compile(r'(?P<alignments>^>[^=]+)=(?P<notes>.*)$',
                                    flags=re.M)

_sub_alignment_re = re.compile(r'^>\s*(?P<seq_num>\d+):(?P<start_idx>\d+)-'
                               r'(?P<end_idx>\d+)\s*(?P<strand>[+|-])\s*'
                               r'(?P<notes>[ \S]+)\s*?$(?P<seq>(?:[^=>]+)+)', 
                               flags=re.M)

_ws_split = re.compile('\s+')

SubAlignmentGroup = namedtuple('SubAlignmentGroup', ['alignments', 'notes'])

SubAlignment = namedtuple('SubAlignment',['seq_num', 'start_idx', 
                                          'end_idx', 'strand', 'notes', 'seq'])


def removeWhitespace(raw_string):
    ''' Remove all whitespace (including newlines/carriage returns) from 
    a string.
    '''
    return ''.join(re.split(_ws_split, raw_string))


def parseXMFA(xmfa_fp):

    with open(xmfa_fp) as xmfa_fd:
        xmfa_data = xmfa_fd.read()

    sub_alignment_groups = re.finditer(_sub_aligment_group_re, xmfa_data)
    sub_alignment_group_list = []

    for sub_alignment_group in sub_alignment_groups:

        sub_alignments = re.finditer(_sub_alignment_re, 
                                     sub_alignment_group.group('alignments'))
        sub_alignment_list = []

        for sub_alignment in sub_alignments:

            sub_alignment_list.append(
                SubAlignment(
                    seq_num=        int(sub_alignment.group('seq_num')),
                    start_idx=      int(sub_alignment.group('start_idx')),
                    end_idx=        int(sub_alignment.group('end_idx')),
                    strand=         sub_alignment.group('strand'),
                    notes=          sub_alignment.group('notes'),
                    seq=            removeWhitespace(sub_alignment.group('seq'))
            ))

        sub_alignment_group_list.append(
            SubAlignmentGroup(
                alignments=         sub_alignment_list,
                notes=              sub_alignment_group.group('notes')
        ))

    return sub_alignment_group_list


if __name__ == '__main__':
    import os

    DIR = os.path.dirname(os.path.realpath(__file__))
    print()
    # output = parseXMFA(os.path.join(DIR, 'example.xmfa'))
    output = parseXMFA(os.path.join(DIR, '1swap100_output', '1swap100.xmfa'))
    print()
    print(output)

