from . import myutils
# Import pandas library
import pandas as pd

class TestFileInput:
    GZ_FILE = 'test-files/founders.vcf.gz'
    NON_GZ_FILE = 'test-files/founders.vcf'
    NAMES = ['#CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 
                     'INFO', 'FORMAT', 'F1', 'F2', 'F3', 'F4']
    VCF_ALL = pd.read_csv(NON_GZ_FILE, comment='#', 
                                  delim_whitespace = True, 
                                  header = None, names = NAMES)
    
    def test_is_gz_true(self):
        assert myutils._is_gz_file(self.GZ_FILE)

    def test_is_gz_false(self):
        assert not myutils._is_gz_file(self.NON_GZ_FILE)
    
    def test_names(self):
        assert myutils._get_vcf_names(self.NON_GZ_FILE) == self.NAMES
    
    def test_names_gz(self):
        assert myutils._get_vcf_names(self.GZ_FILE) == self.NAMES

    def test_load_all(self):
        # should actually system exit
        assert myutils._load_vcf(self.NON_GZ_FILE, None) == self.VCF_ALL
    
    def test_load_Y(self):
        assert myutils._load_vcf(self.NON_GZ_FILE, 'Y').equals(self.VCF_ALL[self.VCF_ALL['#CHROM'] == 'Y'])
    
    def test_load_gz_all(self):
        # should actually system exit
        assert myutils._load_vcf(self.GZ_FILE, None).equals(self.VCF_ALL)
    
    def test_load_gz_Y(self):
        assert myutils._load_vcf(self.GZ_FILE, 'Y').equals(self.VCF_ALL[self.VCF_ALL['#CHROM'] == 'Y'])