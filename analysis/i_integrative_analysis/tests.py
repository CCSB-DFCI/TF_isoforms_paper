import re
from Bio.Alphabet import IUPAC

from data_loading import load_valid_isoform_clones


def test_valid_aa_seqs(df):
    aa_letters = IUPAC.IUPACProtein.letters
    aa_re = re.compile('^[' + aa_letters + ']*$')
    assert df['aa_seq'].str.match(aa_re).all()


valid_clones = load_valid_isoform_clones()
test_valid_aa_seqs(valid_clones)