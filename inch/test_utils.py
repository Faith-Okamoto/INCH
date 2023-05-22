from . import myutils
import pytest
import pandas as pd

class TestFileInput:
    NON_GZ_FOUNDER_FILE = 'test-files/founders.vcf'
    GZ_FOUNDER_FILE = 'test-files/founders.vcf.gz'
    DESC_FILE = 'test-files/descendents.vcf'
    FOUNDER_FILES = [NON_GZ_FOUNDER_FILE, GZ_FOUNDER_FILE]
    
    NONEXIST_FILE = 'test-files/none.vcf'
    BAD_NAMES_FILE = 'test-files/bad_names.vcf'

    FOUNDER_NAMES = ['#CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 
                      'INFO', 'FORMAT', 'F1', 'F2', 'F3', 'F4']
    FOUNDER_VCF = pd.read_csv(NON_GZ_FOUNDER_FILE, comment='#', 
                               delim_whitespace = True, header = None, 
                               names = FOUNDER_NAMES)

    DESC_NAMES = ['#CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 
                  'INFO', 'FORMAT', 'D1', 'D2', 'D3', 'D4']
    DESC_VCF = pd.read_csv(DESC_FILE, comment='#', delim_whitespace = True, 
                           header = None, names = DESC_NAMES)
    
    # TODO: add smaller files to test parts of geno loading
    
    def test_is_gz_true(self):
        assert myutils._is_gz_file(self.GZ_FOUNDER_FILE)

    def test_is_gz_false(self):
        assert not myutils._is_gz_file(self.NON_GZ_FOUNDER_FILE)
    
    def test_is_gz_nonexist_file(self):
        with pytest.raises(FileNotFoundError) as e_info:
            assert myutils._is_gz_file(self.NONEXIST_FILE)
        assert e_info.type == FileNotFoundError
    
    def test_names(self):
        for file in self.FOUNDER_FILES:
            assert myutils._get_vcf_names(file) == self.FOUNDER_NAMES
    
    def test_names_nonexist_file(self):
        with pytest.raises(FileNotFoundError) as e_info:
            assert myutils._get_vcf_names(self.NONEXIST_FILE)
        assert e_info.type == FileNotFoundError
        
    def test_load_single_chr(self):
        for chr in [None, 'Y']:
            assert myutils._load_vcf(self.DESC_FILE, chr).equals(self.DESC_VCF)
    
    def test_load_wrong_single_chr(self):
        with pytest.raises(SystemExit) as e_info:
            myutils._load_vcf(self.DESC_FILE, 'X')
        assert e_info.type == SystemExit
        assert e_info.value.code != 0

    def test_load_multi_chr(self):
        for file in self.FOUNDER_FILES:
            with pytest.raises(SystemExit) as e_info:
                myutils._load_vcf(file, None)
            assert e_info.type == SystemExit
            assert e_info.value.code != 0
    
    def test_load_existing_multi_chr(self):
        for file in self.FOUNDER_FILES:
            for chr in ['Y', 'MT']:
                assert myutils._load_vcf(file, chr).equals(
                    self.FOUNDER_VCF[self.FOUNDER_VCF['#CHROM'] == chr])
    
    def test_load_nonexist_multi_chr(self):
        for file in self.FOUNDER_FILES:
            with pytest.raises(SystemExit) as e_info:
                myutils._load_vcf(file, 'X')
            assert e_info.type == SystemExit
            assert e_info.value.code != 0
    
    def test_load_nonexist_file(self):
        with pytest.raises(FileNotFoundError) as e_info:
            assert myutils._load_vcf(self.NONEXIST_FILE, None)
        assert e_info.type == FileNotFoundError
    
    def test_load_names_malformed(self):
        with pytest.raises(SystemExit) as e_info:
            myutils._load_vcf(self.BAD_NAMES_FILE, None)
        assert e_info.type == SystemExit
        assert e_info.value.code != 0

    #def test_get_geno_single_chr(self):
    #    for chr in [None, 'Y']:
    #        assert myutils._load_vcf(self.DESC_FILE, chr).equals(self.DESC_VCF)
    
    def test_get_geno_wrong_single_chr(self):
        with pytest.raises(SystemExit) as e_info:
            myutils._get_geno(self.DESC_FILE, 'X')
        assert e_info.type == SystemExit
        assert e_info.value.code != 0

    def test_get_geno_multi_chr(self):
        for file in self.FOUNDER_FILES:
            with pytest.raises(SystemExit) as e_info:
                myutils._get_geno(file, None)
            assert e_info.type == SystemExit
            assert e_info.value.code != 0
    
    #def test_get_geno_existing_multi_chr(self):
    #    for file in self.FOUNDER_FILES:
    #        for chr in ['Y', 'MT']:
    #            assert myutils._load_vcf(file, chr).equals(
    #                self.FOUNDER_VCF[self.FOUNDER_VCF['#CHROM'] == chr])
    
    def test_get_geno_multi_chr(self):
        for file in self.FOUNDER_FILES:
            with pytest.raises(SystemExit) as e_info:
                myutils._get_geno(file, 'X')
            assert e_info.type == SystemExit
            assert e_info.value.code != 0
    
    def test_get_geno_nonexist_file(self):
        with pytest.raises(FileNotFoundError) as e_info:
            assert myutils._get_geno(self.NONEXIST_FILE, None)
        assert e_info.type == FileNotFoundError
    
    def test_get_geno_names_malformed(self):
        with pytest.raises(SystemExit) as e_info:
            myutils._get_geno(self.BAD_NAMES_FILE, None)
            assert e_info.type == SystemExit
        assert e_info.value.code != 0