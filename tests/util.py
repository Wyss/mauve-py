
import re

_fasta_re = re.compile(r'>(.*)[\n\r]([^>]+)')

_ALPHABETS = {
    'DNA': 'ACGTWSMKRYBDHVN',
    'RNA': 'ACGUWSMKRYBDHVN',
    'AMINO_ACID': 'ABCDEFGHIJKLMNOPQRSTUVWYZX*-'
}


def sanitizeRecSeq(rec_id, rec_seq, alphabet, not_allowed):
    if alphabet == None and not_allowed:
        for idx, c in enumerate(rec_seq):
            if c in not_allowed:
                raise ValueError('Record: {}\n Contains illegal char {} at '
                                 'position {}'.format(rec_id, c, idx))
    if alphabet:
        ab = _ALPHABETS.get(alphabet.upper(), alphabet)
        if not_allowed:
            for c in not_allowed:
                ab = ab.replace(c, '')
        for idx, c in enumerate(rec_seq):
            if not c in ab:
                raise ValueError('Record: {}\n Contains illegal char {} at '
                                 'position {}'.format(rec_id, c, idx))
# end def

def parseFasta(fasta_fn, alphabet=None, not_allowed=None):
    ''' Parse a standard fasta file and return a list of tuples containing
    the record id, sequence. Optionally check each sequence against an alphabet
    (DNA, RNA, AMINO_ACID, or a custom alphabet) and/or against a list of
    characters that are not allowed:

        alphabet = 'DNA' and not_allowed = 'N' - check against the DNA alphabet
                                                  but do not allow degeneracy
        alphabet = 'ATGC' - custom alphabet (insure that every character in
                            the sequence is A, T, G, or C)

    If a sequence fails the alphabet / not_allowed check a ValueError is
    raised.
    '''
    with open(fasta_fn, 'r') as fd:
        d = [ (match.group(1).strip().split(' ')[0], \
            ''.join(match.group(2).strip().split())) \
                for match in re.finditer(_fasta_re, fd.read()) ]
        if alphabet or not_allowed:
            for rec_id, rec_seq in d:
                sanitizeRecSeq(rec_id, rec_seq, alphabet, not_allowed)
    return d
# end def

def parseFastaGen(fasta_fn, alphabet=None, not_allowed=None):
    ''' Generator that returns parsed records (ID, sequence) from a FASTA
    file.
    '''
    rec_id = ''
    rec_seq = ''
    with open(fasta_fn, 'r') as fd:
        for line in fd:
            if '>' in line:
                if rec_id != '':
                    if alphabet or not_allowed:
                        sanitizeRecSeq(rec_id, rec_seq, alphabet, not_allowed)
                    yield rec_id, rec_seq
                rec_id = line.strip().split(' ')[0][1:]
                rec_seq = ''
            else:
                rec_seq += line.strip()
    if alphabet or not_allowed:
        sanitizeRecSeq(rec_id, rec_seq, alphabet, not_allowed)
    yield rec_id, rec_seq
# end def
