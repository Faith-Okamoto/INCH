from . import myutils
import pytest
import pandas as pd

class TestFileInput:
    NON_GZ_FILE = 'test-files/single_chr.vcf'
    GZ_FILE = 'test-files/single_chr.vcf.gz'
    MULTI_FILE = 'test-files/multi_chr.vcf'
    SINGLE_CHR_FILES = [NON_GZ_FILE, GZ_FILE]
    
    NONEXIST_FILE = 'test-files/none.vcf'
    BAD_NAMES_FILE = 'test-files/bad_names.vcf'
    BAD_FORMAT_FILE = 'test-files/bad_format.vcf'
    BAD_GT_CODE_FILE = 'test-files/bad_gt_code.vcf'

    VCF_NAMES = ['#CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 
             'FILTER', 'INFO', 'FORMAT', 'S1', 'S2']
    SAMPLE_NAMES = ['S1', 'S2']
    
    SINGLE_CHR_VCF = [
        ['Y', 1, '.', 'A', 'T', 20, 'PASS', 'AF=0.5', 'GT', 0, 1],
        ['Y', 2, '.', 'C', 'G', 20, 'PASS', 'AF=0.5', 'GT', 0, 1]
    ]
    SINGLE_CHR_VCF = pd.DataFrame(SINGLE_CHR_VCF, columns = VCF_NAMES)
    SINGLE_CHR_GENO = pd.DataFrame([[0, 3], [1, 2]], columns = SAMPLE_NAMES)
    SINGLE_CHR_POS = pd.Series([1, 2])

    MULTI_CHR_VCF = [
        ['Y', 1, '.', 'A', '*,T', 20, 'PASS', 'AF=0', 'GT', './1', '1'],
        ['MT', 2, '.', 'CGC', 'G,CC,CAC', 20, 
         'PASS', 'AF=0.5', 'GT:GQ', '3|3:.', '2:.']
    ]
    MULTI_CHR_VCF = pd.DataFrame(MULTI_CHR_VCF, columns = VCF_NAMES)
    MULTI_CHR_GENO = {'Y' : pd.DataFrame([[-10, -1]], columns = SAMPLE_NAMES),
                      'MT' : pd.DataFrame([[17, 5]], columns = SAMPLE_NAMES)}
    MULTI_CHR_POS = {'Y' : pd.Series([1]), 'MT' : pd.Series([2])}
    
    def test_is_gz_true(self):
        assert myutils._is_gz_file(self.GZ_FILE)

    def test_is_gz_false(self):
        assert not myutils._is_gz_file(self.NON_GZ_FILE)
    
    def test_is_gz_nonexist_file(self):
        with pytest.raises(FileNotFoundError) as e_info:
            assert myutils._is_gz_file(self.NONEXIST_FILE)
        assert e_info.type == FileNotFoundError
    
    def test_names(self):
        for file in self.SINGLE_CHR_FILES:
            assert myutils._get_vcf_names(file) == self.VCF_NAMES
    
    def test_names_nonexist_file(self):
        with pytest.raises(FileNotFoundError) as e_info:
            assert myutils._get_vcf_names(self.NONEXIST_FILE)
        assert e_info.type == FileNotFoundError
    
    def test_names_malformed(self):
        with pytest.raises(SystemExit) as e_info:
            myutils._get_vcf_names(self.BAD_NAMES_FILE)
        assert e_info.type == SystemExit
        
    def test_load_single_chr(self):
        for file in self.SINGLE_CHR_FILES:
            for chr in [None, 'Y']:
                assert myutils._load_vcf(file, chr).equals(self.SINGLE_CHR_VCF)
    
    def test_load_wrong_single_chr(self):
        for file in self.SINGLE_CHR_FILES:
            with pytest.raises(SystemExit) as e_info:
                myutils._load_vcf(file, 'X')
            assert e_info.type == SystemExit

    def test_load_all_multi_chr(self):
        with pytest.raises(SystemExit) as e_info:
            myutils._load_vcf(self.MULTI_FILE, None)
        assert e_info.type == SystemExit   
    
    def test_load_existing_multi_chr(self):
        for chr in ['Y', 'MT']:
            correct_geno = self.MULTI_CHR_VCF[
                self.MULTI_CHR_VCF['#CHROM'] == chr
            ].reset_index(drop = True)
            assert myutils._load_vcf(self.MULTI_FILE, chr).equals(correct_geno)
    
    def test_load_nonexist_multi_chr(self):
        with pytest.raises(SystemExit) as e_info:
            myutils._load_vcf(self.MULTI_FILE, 'X')
        assert e_info.type == SystemExit
    
    def test_load_nonexist_file(self):
        with pytest.raises(FileNotFoundError) as e_info:
            assert myutils._load_vcf(self.NONEXIST_FILE, None)
        assert e_info.type == FileNotFoundError
    
    def test_load_names_malformed(self):
        with pytest.raises(SystemExit) as e_info:
            myutils._load_vcf(self.BAD_NAMES_FILE, None)
        assert e_info.type == SystemExit

    def test_get_geno_single_chr(self):
        for file in self.SINGLE_CHR_FILES:
            for chr in [None, 'Y']:
                geno, pos, chr_found = myutils._get_geno(file, chr)
                assert geno.equals(self.SINGLE_CHR_GENO)
                assert pos.equals(self.SINGLE_CHR_POS)
                assert chr_found == 'Y'
    
    def test_get_geno_wrong_single_chr(self):
        for file in self.SINGLE_CHR_FILES:
            with pytest.raises(SystemExit) as e_info:
                myutils._get_geno(file, 'X')
            assert e_info.type == SystemExit
    
    def test_get_geno_all_multi_chr(self):
        with pytest.raises(SystemExit) as e_info:
            myutils._get_geno(self.MULTI_FILE, None)
        assert e_info.type == SystemExit   
    
    def test_get_geno_existing_multi_chr(self):
        for chr in ['Y', 'MT']:
            geno_found, pos, chr_found = myutils._get_geno(self.MULTI_FILE, chr)
            correct_geno = self.MULTI_CHR_GENO[chr]
            assert geno_found.values.tolist() == correct_geno.values.tolist()
            assert geno_found.columns.equals(correct_geno.columns)
            assert pos.equals(self.MULTI_CHR_POS[chr])
            assert chr_found == chr
    
    def test_get_geno_nonexist_multi_chr(self):
        with pytest.raises(SystemExit) as e_info:
            myutils._get_geno(self.MULTI_FILE, 'X')
        assert e_info.type == SystemExit
    
    def test_get_geno_names_malformed(self):
        with pytest.raises(SystemExit) as e_info:
            myutils._get_geno(self.BAD_NAMES_FILE, None)
        assert e_info.type == SystemExit
    
    def test_get_geno_bad_format(self):
        with pytest.raises(SystemExit) as e_info:
            myutils._get_geno(self.BAD_FORMAT_FILE, None)
        assert e_info.type == SystemExit
    
    def test_get_geno_bad_gt_code(self):
        with pytest.raises(SystemExit) as e_info:
            myutils._get_geno(self.BAD_GT_CODE_FILE, None)
        assert e_info.type == SystemExit

class TestUtilities:
    # test _to_code, _make_groups, _merge_matrix_groups, _geno_dists, ERROR
    pass

class TestAnalysis:
    # test dist_matrix, pca, identify_founders 
    # use test-files/founders.vcf and test-files/descendents.vcf
    pass