#!/usr/bin/env python3

import json
import os
import unittest

import auto_irida_azure_upload.core as core

TEST_DATA_DIR = os.path.abspath(os.path.join(os.path.dirname(__file__), 'data'))

class TestFindRunDirs(unittest.TestCase):
    """
    Test the core.find_run_dirs function
    """

    def test_find_run_dirs(self):
        """
        Test that the function returns a list of run directories
        """
        self.maxDiff = None
        config = {
            'run_parent_dirs': [
                os.path.join(TEST_DATA_DIR, 'run_parent_dirs', 'M00123', '22'),
            ],
            'upload_staging_dir': os.path.join(TEST_DATA_DIR, 'upload_staging'),
        }

        expected_run_dirs = [
            {
                'path': os.path.join(TEST_DATA_DIR, 'run_parent_dirs/M00123/22/220101_M00123_0001_000000000-AAA11'),
                'sequencing_run_id': '220101_M00123_0001_000000000-AAA11',
                'instrument_type': 'miseq',
            },
            {
                'path': os.path.join(TEST_DATA_DIR, 'run_parent_dirs/M00123/22/220112_M00123_0002_000000000-AAA12'),
                'sequencing_run_id': '220112_M00123_0002_000000000-AAA12',
                'instrument_type': 'miseq',
            },
        ]
        
        actual_run_dirs = []
        for run_dir in core.find_run_dirs(config):
            actual_run_dirs.append(run_dir)

        self.assertCountEqual(expected_run_dirs, actual_run_dirs)


class TestScan(unittest.TestCase):
    """
    Test the core.scan function.
    """
    def test_scan(self):
        """
        Test that the function returns a list of runs
        """
        # self.maxDiff = None
        config = {
            'run_parent_dirs': [
                os.path.join(TEST_DATA_DIR, 'run_parent_dirs', 'M00123', '22'),
            ],
            'upload_staging_dir': os.path.join(TEST_DATA_DIR, 'upload_staging'),
        }

        expected_runs = [
            {
                'path': os.path.join(TEST_DATA_DIR, 'run_parent_dirs/M00123/22/220101_M00123_0001_000000000-AAA11'),
                'sequencing_run_id': '220101_M00123_0001_000000000-AAA11',
                'instrument_type': 'miseq',
            },
            {
                'path': os.path.join(TEST_DATA_DIR, 'run_parent_dirs/M00123/22/220112_M00123_0002_000000000-AAA12'),
                'sequencing_run_id': '220112_M00123_0002_000000000-AAA12',
                'instrument_type': 'miseq',
            },
        ]

        actual_runs = []
        for run in core.scan(config):
            actual_runs.append(run)

        self.assertCountEqual(expected_runs, actual_runs)

        
class TestFindFastq(unittest.TestCase):
    """
    """
    def test_find_fastq(self):
        """
        Test that the function returns the correct fastq path
        """
        run = {
            'path': os.path.join(TEST_DATA_DIR, 'run_parent_dirs/M00123/22/220101_M00123_0001_000000000-AAA11'),
            'instrument_type': 'miseq',
            'sequencing_run_id': '220101_M00123_0001_000000000-AAA11',
        }
        library_id = 'sample-01-A01'
        read_type = 'R1'
        expected_fastq_path = os.path.join(
            TEST_DATA_DIR,
            'run_parent_dirs/M00123/22/220101_M00123_0001_000000000-AAA11',
            'Data/Intensities/BaseCalls',
            library_id + '_S1_L001_R1_001.fastq.gz'
        )
        actual_fastq_path = core.find_fastq(run, library_id, read_type)

        self.assertEqual(expected_fastq_path, actual_fastq_path)


class TestPrepareDownsamplingInputs(unittest.TestCase):
    """
    """
    def test_prepare_downsampling_samplesheet(self):
        """
        Test that the function returns the correct downsampling inputs
        """
        self.maxDiff = None
        config = {
            'run_parent_dirs': [
                os.path.join(TEST_DATA_DIR, 'run_parent_dirs', 'M00123', '22'),
            ],
            'upload_staging_dir': os.path.join(TEST_DATA_DIR, 'upload_staging'),
            'projects': [
                {
                    'local_project_id': 'data_sharing',
                    'local_project_name': 'Data Sharing',
                    'remote_project_id': '42',
                    'remote_project_name': 'BC Data',
                    'downsample_reads': True,
                    'genome_size_mb': 5.0,
                    'max_depth': 50,
                },
                {
                    'local_project_id': 'new_project',
                    'local_project_name': 'New Project',
                    'remote_project_id': '48',
                    'remote_project_name': 'BC New Test Project',
                    'downsample_reads': False,
                    'genome_size_mb': 5.0,
                    'max_depth': 50,
                },
            ]
        }
        run = {
            'sequencing_run_id': '220101_M00123_0001_000000000-AAA11',
            'instrument_type': 'miseq',
            'path': os.path.join(TEST_DATA_DIR, 'run_parent_dirs/M00123/22/220101_M00123_0001_000000000-AAA11'),
        }
        expected_downsampling_samplesheet = [
            {
                'ID': 'sample-01-A01',
                'R1': os.path.join(
                    TEST_DATA_DIR,
                    'run_parent_dirs/M00123/22/220101_M00123_0001_000000000-AAA11',
                    'Data/Intensities/BaseCalls',
                    'sample-01-A01_S1_L001_R1_001.fastq.gz'
                ),
                'R2': os.path.join(
                    TEST_DATA_DIR,
                    'run_parent_dirs/M00123/22/220101_M00123_0001_000000000-AAA11',
                    'Data/Intensities/BaseCalls',
                    'sample-01-A01_S1_L001_R2_001.fastq.gz'
                ),
                'GENOME_SIZE': '5.0m',
                'COVERAGE': '50',
            },
            {
                'ID': 'sample-02-B01',
                'R1': os.path.join(
                    TEST_DATA_DIR,
                    'run_parent_dirs/M00123/22/220101_M00123_0001_000000000-AAA11',
                    'Data/Intensities/BaseCalls',
                    'sample-02-B01_S2_L001_R1_001.fastq.gz'
                ),
                'R2': os.path.join(
                    TEST_DATA_DIR,
                    'run_parent_dirs/M00123/22/220101_M00123_0001_000000000-AAA11',
                    'Data/Intensities/BaseCalls',
                    'sample-02-B01_S2_L001_R2_001.fastq.gz'
                ),
                'GENOME_SIZE': '5.0m',
                'COVERAGE': '50',
            },
            {
                'ID': 'sample-05-E01',
                'R1': os.path.join(
                    TEST_DATA_DIR,
                    'run_parent_dirs/M00123/22/220101_M00123_0001_000000000-AAA11',
                    'Data/Intensities/BaseCalls',
                    'sample-05-E01_S5_L001_R1_001.fastq.gz'
                ),
                'R2': os.path.join(
                    TEST_DATA_DIR,
                    'run_parent_dirs/M00123/22/220101_M00123_0001_000000000-AAA11',
                    'Data/Intensities/BaseCalls',
                    'sample-05-E01_S5_L001_R2_001.fastq.gz'
                ),
                'GENOME_SIZE': '5.0m',
                'COVERAGE': '50',
            },
            {
                'ID': 'sample-06-F01',
                'R1': os.path.join(
                    TEST_DATA_DIR,
                    'run_parent_dirs/M00123/22/220101_M00123_0001_000000000-AAA11',
                    'Data/Intensities/BaseCalls',
                    'sample-06-F01_S6_L001_R1_001.fastq.gz'
                ),
                'R2': os.path.join(
                    TEST_DATA_DIR,
                    'run_parent_dirs/M00123/22/220101_M00123_0001_000000000-AAA11',
                    'Data/Intensities/BaseCalls',
                    'sample-06-F01_S6_L001_R2_001.fastq.gz'
                ),
                'GENOME_SIZE': '5.0m',
                'COVERAGE': '50',
            },
        ]

        actual_downsampling_samplesheet = core.prepare_downsampling_samplesheet(config, run)

        self.assertEqual(expected_downsampling_samplesheet, actual_downsampling_samplesheet)


if __name__ == '__main__':
    unittest.main()
