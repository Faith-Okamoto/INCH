from . import myutils
import pytest
import pandas as pd
from itertools import combinations, product

class TestUtilities:
    MULTI_FILE = 'test-files/multi_chr.vcf'
    SINGLE_CHR_FILES = ['test-files/single_chr.vcf', 
                        'test-files/single_chr.vcf.gz']
    SAMPLE_NAMES = ['S1', 'S2']
    SINGLE_CHR_GENO = pd.DataFrame([[1, 4], [2, 3]], columns = SAMPLE_NAMES,
                                   index = pd.Series([1, 2]))

    MULTI_CHR_GENO = {'Y' : pd.DataFrame([[-2, -1]], columns = SAMPLE_NAMES,
                                         index = pd.Series([1])),
                      'MT' : pd.DataFrame([[57, 12]], columns = SAMPLE_NAMES,
                                          index = pd.Series([2]))}
    
    IDS = ['F1', 'F2', 'F3']
    GENOS = pd.DataFrame([[1, 1, 1], [1, 1, 2], [1, 1, 2], [1, 2, 2]], 
                         columns = IDS)
    MATRIX = pd.DataFrame([[0, 0.25, 0.75], [0.25, 0, 0.5], [0.75, 0.5, 0]],
                          index = IDS, columns = IDS)
    IDS = set(IDS)

    def test_get_geno_single_chr(self):
        for file in self.SINGLE_CHR_FILES:
            for chr in [None, 'Y']:
                geno, chr_found = myutils._get_geno(file, chr)
                assert (geno.values == self.SINGLE_CHR_GENO.values).all()
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
            geno_found, chr_found = myutils._get_geno(self.MULTI_FILE, chr)
            correct_geno = self.MULTI_CHR_GENO[chr]
            assert (geno_found.values == correct_geno.values).all()
            assert geno_found.columns.equals(correct_geno.columns)
            assert chr_found == chr
    
    def test_get_geno_nonexist_multi_chr(self):
        with pytest.raises(SystemExit) as e_info:
            myutils._get_geno(self.MULTI_FILE, 'X')
        assert e_info.type == SystemExit

    def test_error(self):
        for msg in ['ERROR', 'test message', '']:
            with pytest.raises(SystemExit) as e_info:
                myutils.ERROR(msg)
            assert e_info.type == SystemExit
    
    def test_to_code_missing(self):
        assert myutils._to_code("*") == -1
    
    def test_to_code_single(self):
        for base, code in myutils.BASES.items():
            assert myutils._to_code(base) == code
    
    def test_to_code_multi(self):
        assert myutils._to_code('AA') == 6
        assert myutils._to_code('AC') == 7
        assert myutils._to_code('GCT') == 89
    
    def test_to_code_illegal(self):
        for allele in ['.', 'A*', '2']:
            with pytest.raises(SystemExit) as e_info:
                myutils._to_code(allele)
            assert e_info.type == SystemExit
    
    def test_make_groups_nomissing(self):
        assert myutils._make_groups(self.IDS, ['F1', 'F2', 'F3']).sort() == \
            [['F1'], ['F2'], ['F3']].sort() 
        assert myutils._make_groups(self.IDS, ['F1,F2', 'F3']).sort() == \
            [['F1', 'F2'], ['F3']].sort()
        assert myutils._make_groups(self.IDS, ['F3,F1,F2']).sort() == \
            [['F3', 'F1', 'F2']].sort()
    
    def test_make_groups_missing(self):
        assert myutils._make_groups(self.IDS, ['F1', 'F2']).sort() == \
            [['F1'], ['F2'], ['F3']].sort()
        assert myutils._make_groups(self.IDS, ['F1,F2']).sort() == \
            [['F1', 'F2'], ['F3']].sort()
        assert myutils._make_groups(self.IDS, []).sort() == \
            [['F1'], ['F2'], ['F3']].sort()
    
    def test_make_groups_duplicate(self):
        dup_groups = [
            ['F1', 'F1', 'F2', 'F3'], ['F2', 'F2'], ['F2,F2'],
            ['F1,F3', 'F3'], ['F2,F3', 'F1,F3']
        ]
        for groups in dup_groups:
            with pytest.raises(SystemExit) as e_info:
                myutils._make_groups(self.IDS, groups)
            assert e_info.type == SystemExit

    def test_make_groups_nonexist(self):
        nonexist_groups = [['X'], ['F1,F2', 'F3', 'X'], ['F2,', 'F1'], ['F3,X']]
        for groups in nonexist_groups:
            with pytest.raises(SystemExit) as e_info:
                myutils._make_groups(self.IDS, groups)
            assert e_info.type == SystemExit
    
    def test_merge_singles(self):
        merged = myutils._merge_matrix_groups(
            self.MATRIX, [['F1'], ['F2'], ['F3']]
        )
        assert merged.shape == (3, 3)
        for i in self.IDS:
            for j in self.IDS:
                assert merged.loc[i, j] == self.MATRIX.loc[i, j]
    
    def test_merge_pair(self):
        merged = myutils._merge_matrix_groups(
            self.MATRIX, [['F1', 'F2'], ['F3']]
        )
        assert merged.shape == (2, 2)
        assert merged.loc['F1,F2', 'F1,F2'] == 0
        assert merged.loc['F3', 'F3'] == 0
        assert merged.loc['F1,F2', 'F3'] == 0.625
        assert merged.loc['F3', 'F1,F2'] == 0.625
    
    def test_merge_all(self):
        merged = myutils._merge_matrix_groups(self.MATRIX, [['F1', 'F2', 'F3']])
        assert merged.shape == (1, 1)
        assert merged.loc['F1,F2,F3', 'F1,F2,F3'] == 0
    
    def test_geno_dist_identical(self):
        dists = myutils._geno_dists(self.GENOS[['F1']], self.GENOS[['F1']])
        assert dists.shape == (1, 1)
        assert dists.loc['F1', 'F1'] == 0

        F1_copy = self.GENOS[['F1']]
        F1_copy['COPY'] = F1_copy['F1'].copy()
        dists = myutils._geno_dists(self.GENOS[['F1']], F1_copy)
        assert dists.shape == (1, 2)
        assert dists.loc['F1', 'F1'] == 0
        assert dists.loc['F1', 'COPY'] == 0
    
    def test_geno_dist_nonidentical(self):
        dists = myutils._geno_dists(self.GENOS, self.GENOS)
        for i in self.IDS:
            for j in self.IDS:
                assert dists.loc[i, j] == self.MATRIX.loc[i, j]

class TestAnalysis:
    FOUNDER_FILES = ['test-files/founders.vcf', 'test-files/founders.vcf.gz']
    DESC_FILE = 'test-files/descendents.vcf'
    FOUNDER_MAX_PCS = {'Y' : 4, 'MT': 3}
    CORRECT_ID = {'D1': 'F1', 'D2': 'F2', 'D3': 'F3', 'D4': 'F3'}
    FOUNDER_GROUPS = [['F1,F2,F3,F4'], ['F1', 'F2', 'F4'], ['F1,F2', 'F3']]
    FOUNDER_BAD_GROUPS = [['F1,F2,F3,F3'], ['F1', 'X', 'F4'], ['F1,F2', 'F1']]

    # test dist_matrix
    def test_identify_all_multi_chr(self):
        for f, d in product(self.FOUNDER_FILES + [self.DESC_FILE], repeat = 2): 
            if f != self.DESC_FILE or d != self.DESC_FILE:
                with pytest.raises(SystemExit) as e_info:
                    myutils.identify_founders(f, d, None, None)
                assert e_info.type == SystemExit
    
    def test_identify_identical_single_chr(self):
        for chr in ['Y', None]:
            ids = myutils.identify_founders(self.DESC_FILE, self.DESC_FILE, 
                                            chr, None)
            for i in range(len(ids)):
                assert ids[i] == ids.index[i]
    
    def test_identify_identical_single_chr_groups(self):
        for groups in [['D1,D2,D3,D4'], ['D1', 'D2', 'D4'], ['D1,D2', 'D3']]:
            for chr in ['Y', None]:
                ids = myutils.identify_founders(self.DESC_FILE, self.DESC_FILE, 
                                                chr, groups)
                for i in range(len(ids)):
                    assert ids.index[i] in ids[i]
    
    def test_identify_identical_single_chr_bad_groups(self):
        for groups in [['D1,D2,D3,D3'], ['D1', 'X', 'D4'], ['D1,D2', 'D1']]:
            for chr in ['Y', None]:
                with pytest.raises(SystemExit) as e_info:
                    myutils.identify_founders(self.DESC_FILE, self.DESC_FILE, 
                                              chr, groups)
                assert e_info.type == SystemExit
    
    def test_identify_identical_multi_chr(self):
        for f, d in product(self.FOUNDER_FILES, repeat = 2):
            for chr in ['Y', 'MT']:
                ids = myutils.identify_founders(f, d, chr, None)
                for i in range(len(ids)):
                    assert ids[i] == ids.index[i]
    
    def test_identify_identical_multi_chr_groups(self):
        for f, d in product(self.FOUNDER_FILES, repeat = 2):
            for groups in self.FOUNDER_GROUPS:
                for chr in ['Y', 'MT']:
                    ids = myutils.identify_founders(f, d, chr, groups)
                    for i in range(len(ids)):
                        assert ids.index[i] in ids[i]
    
    def test_identify_identical_multi_chr_bad_groups(self):
        for f, d in product(self.FOUNDER_FILES, repeat = 2):
            for groups in self.FOUNDER_BAD_GROUPS:
                for chr in ['Y', 'MT', 'X', None]:
                    with pytest.raises(SystemExit) as e_info:
                        myutils.identify_founders(f, d, chr, groups)
                    assert e_info.type == SystemExit

    def test_identify_diff(self):
        for f_file in self.FOUNDER_FILES:
            ids = myutils.identify_founders(f_file, self.DESC_FILE, 'Y', None)
            for d, f in self.CORRECT_ID.items():
                assert ids[d] == f
        
    def test_identify_diff_groups(self):
        for f_file in self.FOUNDER_FILES:
            for groups in self.FOUNDER_GROUPS:
                ids = myutils.identify_founders(f_file, self.DESC_FILE, 
                                                'Y', groups)
                for d, f in self.CORRECT_ID.items():
                    assert f in ids[d]
    
    def test_identify_diff_bad_groups(self):
        for f_file in self.FOUNDER_FILES:
            for groups in self.FOUNDER_BAD_GROUPS:
                for chr in ['Y', 'MT', 'X', None]:
                    with pytest.raises(SystemExit) as e_info:
                        myutils.identify_founders(f_file, self.DESC_FILE, 
                                                  chr, groups)
                    assert e_info.type == SystemExit
    
    def test_identify_diff_wrong_chr(self):
        for groups in self.FOUNDER_GROUPS:
            for chr in ['MT', 'X']:
                for f_file in self.FOUNDER_FILES:
                    with pytest.raises(SystemExit) as e_info:
                        myutils.identify_founders(f_file, self.DESC_FILE, 
                                                  chr, groups)
                    assert e_info.type == SystemExit

    def test_pca_single_chr(self):
        for chr in [None, 'Y']:
            for n in range(1, 5):
                evectors, evalues = myutils.pca(self.DESC_FILE, chr, n)
                assert evectors.shape == (4, n)
                assert len(evalues) == n
                for i in range(1, n):
                    assert evalues[i - 1] > evalues[i]
                d3_d4 = abs(evectors.loc['D3', 'PC1'] - 
                            evectors.loc['D4', 'PC1'])
                for i, j in combinations(['D1', 'D2', 'D3', 'D4'], 2):
                    assert d3_d4 <= abs(evectors.loc[i, 'PC1'] - 
                                        evectors.loc[j, 'PC1'])
    
    def test_pca_wrong_single_chr(self):
        for n in range(1, 5):
            with pytest.raises(SystemExit) as e_info:
                myutils.pca(self.DESC_FILE, 'X', n)
            assert e_info.type == SystemExit
    
    def test_pca_bad_n_single_chr(self):
        for n in [-2, -1, 0, 5, 6, 7]:
            with pytest.raises(SystemExit) as e_info:
                myutils.pca(self.DESC_FILE, 'Y', n)
            assert e_info.type == SystemExit
    
    def test_pca_all_multi_chr(self):
        for file in self.FOUNDER_FILES:
            for n in range(1, 4):
                with pytest.raises(SystemExit) as e_info:
                    myutils.pca(file, None, n)
                assert e_info.type == SystemExit   
    
    def test_pca_existing_multi_chr(self):
        for file in self.FOUNDER_FILES:
            for chr in ['Y', 'MT']:
                for n in range(1, self.FOUNDER_MAX_PCS[chr]):
                    evectors, evalues = myutils.pca(file, chr, n)
                    assert evectors.shape == (4, n)
                    assert len(evalues) == n
                    for i in range(1, n):
                        assert evalues[i - 1] > evalues[i]
                    f3_f4 = abs(evectors.loc['F3', 'PC1'] - 
                                evectors.loc['F4', 'PC1'])
                    for i, j in combinations(['F1', 'F2', 'F3', 'F4'], 2):
                        assert f3_f4 <= abs(evectors.loc[i, 'PC1'] - 
                                            evectors.loc[j, 'PC1'])
    
    def test_pca_nonexist_multi_chr(self):
        for file in self.FOUNDER_FILES:
            for n in range(1, 4):
                with pytest.raises(SystemExit) as e_info:
                    myutils.pca(file, 'X', n)
                assert e_info.type == SystemExit
    
    def test_pca_bad_n_multi_chr(self):
        for file in self.FOUNDER_FILES:
            for chr in ['Y', 'MT']:
                max_pcs = min(self.FOUNDER_MAX_PCS[chr], 4)
                for n in [-2, -1, 0] + list(range(max_pcs + 1, max_pcs + 4)):
                    with pytest.raises(SystemExit) as e_info:
                        myutils.pca(file, chr, n)
                    assert e_info.type == SystemExit