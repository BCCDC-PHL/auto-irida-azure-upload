#!/usr/bin/env python3

import unittest

import auto_irida_azure_upload.samplesheet as samplesheet

class SampleSheetTest(unittest.TestCase):

    def test_parse_samplesheet_correct_num_samples(self):
        samplesheet_path = 'tests/data/run_parent_dirs/M00123/22/220101_M00123_0001_000000000-AAA11/SampleSheet.csv'
        parsed_samplesheet = samplesheet.parse_samplesheet(samplesheet_path, 'miseq')

        sample_data = parsed_samplesheet.get('data', [])
        actual_num_samples = len(sample_data)
        expected_num_samples = 6
        self.assertEqual(actual_num_samples, expected_num_samples)

if __name__ == '__main__':
    unittest.main()
