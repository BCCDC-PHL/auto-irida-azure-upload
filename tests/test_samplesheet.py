import unittest
import os
import json

import auto_irida_azure_upload.samplesheet as samplesheet

TEST_DATA_DIR = os.path.abspath(os.path.join(os.path.dirname(__file__), 'data'))

class TestFindSampleSheets(unittest.TestCase):
    """
    Test the samplesheet.find_sample_sheets function
    """

    def test_find_samplesheets(self):
        """
        Test that the function returns a list of sample sheets
        """
        config = {
            'run_parent_dirs': [
                os.path.join(TEST_DATA_DIR, 'run_parent_dirs', 'M00123', '22'),
            ],
            'upload_staging_dir': os.path.join(TEST_DATA_DIR, 'upload_staging'),
        }

        sequencing_run = {
            'path': os.path.join(TEST_DATA_DIR, 'run_parent_dirs/M00123/22/220101_M00123_0001_000000000-AAA11'),
            'sequencing_run_id': '220101_M00123_0001_000000000-AAA11',
            'instrument_type': 'miseq',
        }
        expected_samplesheet_paths = [
            os.path.join(TEST_DATA_DIR, 'run_parent_dirs/M00123/22/220101_M00123_0001_000000000-AAA11/SampleSheet.csv'),
        ]
        actual_samplesheet_paths = samplesheet.find_samplesheets(sequencing_run['path'], sequencing_run['instrument_type'])

        self.assertCountEqual(expected_samplesheet_paths, actual_samplesheet_paths)

if __name__ == '__main__':
    unittest.main()
