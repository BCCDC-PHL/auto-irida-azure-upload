import datetime
import json
import logging
import os
import re
import shutil
import subprocess
import time
import uuid

from typing import Iterator, Optional

import auto_irida_upload.samplesheet as samplesheet


def find_run_dirs(config, check_symlinks_complete=True):
    """
    """
    miseq_run_id_regex = "\d{6}_M\d{5}_\d+_\d{9}-[A-Z0-9]{5}"
    nextseq_run_id_regex = "\d{6}_VH\d{5}_\d+_[A-Z0-9]{9}"
    run_parent_dirs = config['run_parent_dirs']

    for run_parent_dir in run_parent_dirs:
        subdirs = os.scandir(run_parent_dir)

        for subdir in subdirs:
            run_id = subdir.name
            matches_miseq_regex = re.match(miseq_run_id_regex, run_id)
            matches_nextseq_regex = re.match(nextseq_run_id_regex, run_id)
            instrument_type = 'unknown'
            if matches_miseq_regex:
                instrument_type = 'miseq'
            elif matches_nextseq_regex:
                instrument_type = 'nextseq'
            ready_to_upload = os.path.exists(os.path.join(subdir, 'upload_complete.json'))
            upload_not_already_initiated = not os.path.exists(os.path.join(subdir, 'irida_upload_started.json'))
            conditions_checked = {
                "is_directory": subdir.is_dir(),
                "matches_illumina_run_id_format": ((matches_miseq_regex is not None) or
                                                   (matches_nextseq_regex is not None)),
                "ready_to_upload": ready_to_upload,
                "upload_not_already_initiated": upload_not_already_initiated,
            }
            conditions_met = list(conditions_checked.values())
            run_dir = {}
            if all(conditions_met):
                logging.info(json.dumps({"event_type": "run_directory_found", "sequencing_run_id": run_id, "run_directory_path": os.path.abspath(subdir.path)}))
                run_dir['path'] = os.path.abspath(subdir.path)
                run_dir['sequencing_run_id'] = run_id
                run_dir['instrument_type'] = instrument_type
            yield run_dir
        else:
            logging.debug(json.dumps({"event_type": "directory_skipped", "run_directory_path": os.path.abspath(subdir.path), "conditions_checked": conditions_checked}))
            yield None
    

def scan(config: dict[str, object]) -> Iterator[Optional[dict[str, object]]]:
    """
    Scanning involves looking for all existing runs and storing them to the database,
    then looking for all existing symlinks and storing them to the database.
    At the end of a scan, we should be able to determine which (if any) symlinks need to be created.

    :param config: Application config.
    :type config: dict[str, object]
    :return: A run directory to analyze, or None
    :rtype: Iterator[Optional[dict[str, object]]]
    """
    logging.info(json.dumps({"event_type": "scan_start"}))
    for run_dir in find_run_dirs(config):    
        yield run_dir


def prepare_samplelist_file(config, run):
    """
    
    """
    samplelist_header = [
        '[Data]',
        'Sample_Name,Project_ID,File_Forward,File_Reverse',
    ]
    samples_to_upload = []
    run_id = run['sequencing_run_id']
    run_dir = run['path']
    instrument_type = run['instrument_type']

    run_samplesheets = samplesheet.find_samplesheets(run_dir, instrument_type)
    samplesheet_to_parse = samplesheet.choose_samplesheet_to_parse(run_samplesheets, instrument_type, run_id)
    if os.path.exists(samplesheet_to_parse):
        parsed_samplesheet = samplesheet.parse_samplesheet(samplesheet_to_parse, instrument_type)
    if instrument_type == 'nextseq':
        if 'cloud_data' in parsed_samplesheet:
            for library in parsed_samplesheet['cloud_data']:
                if 'project_name' in library:
                    for project in config['projects']:
                        if library['project_name'] == project['local_project_id']:
                            samplelist_library = {}
                            samplelist_library['Sample_Name'] = library['sample_id']
                            samplelist_library['Project_ID'] = project['remote_project_id']
                            samples_to_upload.append(samplelist_library)
    elif instrument_type == 'miseq':
        if 'data' in parsed_samplesheet:
            for library in parsed_samplesheet['data']:
                if 'sample_project' in library:
                    for project in config['projects']:
                        if library['sample_project'] == project['local_project_id']:
                            samplelist_library = {}
                            samplelist_library['Sample_Name'] = library['sample_id']
                            samplelist_library['Project_ID'] = project['remote_project_id']
                            samples_to_upload.append(samplelist_library)

    if len(samples_to_upload) > 0:
        upload_tmpdir = config['upload_tmpdir']
        samplelist_path = os.path.join(upload_tmpdir, run_id, 'SampleList.csv')
        if os.path.exists(upload_tmpdir):
            run_upload_tmpdir = os.path.join(upload_tmpdir, run_id)
            if os.path.exists(run_upload_tmpdir):
                shutil.rmtree(run_upload_tmpdir)
            os.makedirs(run_upload_tmpdir)
            with open(samplelist_path, 'w') as f:
                for line in samplelist_header:
                    f.write(line + '\n')
            with open(samplelist_path, 'a') as f:
                for sample in samples_to_upload:
                    f.write(','.join([sample['Sample_Name'], sample['Project_ID']]) + '\n')
        

def upload_run(config, run):
    """
    Initiate an analysis on one directory of fastq files.
    """
    run_id = run['sequencing_run_id']
    
    azcopy_command = [
        'azcopy',
        'cp',
        '--recursive',
        '--folow-symlinks',
    ]
    
    logging.info(json.dumps({"event_type": "upload_started", "sequencing_run_id": run_id, "azcopy_command": " ".join(azcopy_command)}))
    try:
        # subprocess.run(azcopy_command, capture_output=True, check=True)
        time.sleep(10)
        logging.info(json.dumps({"event_type": "upload_completed", "sequencing_run_id": run_id, "azcopy_command": " ".join(azcopy_command)}))
    except subprocess.CalledProcessError as e:
        logging.error(json.dumps({"event_type": "upload_failed", "sequencing_run_id": run_id, "azcopy_command": " ".join(azcopy_command)}))
    except OSError as e:
        logging.error(json.dumps({"event_type": "delete_analysis_work_dir_failed", "sequencing_run_id": analysis_run_id, "analysis_work_dir_path": analysis_work_dir}))
