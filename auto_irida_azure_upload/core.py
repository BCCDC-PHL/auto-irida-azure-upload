import csv
import datetime
import glob
import hashlib
import json
import logging
import os
import re
import shutil
import subprocess
import time
import uuid

from typing import Iterator, Optional

import auto_irida_azure_upload.samplesheet as samplesheet


def find_run_dirs(config, check_upload_complete=True):
    """
    Find sequencing run directories under the 'run_parent_dirs' listed in the config.

    :param config: Application config.
    :type config: dict[str, object]
    :param check_upload_complete: Check for presence of 'upload_complete.json' file.
    :type check_upload_complete: bool
    :return: Run directory. Keys: ['sequencing_run_id', 'path', 'instrument_type']
    :rtype: Iterator[Optional[dict[str, str]]]
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
            upload_complete = os.path.exists(os.path.join(subdir, 'upload_complete.json'))
            qc_check_complete = os.path.exists(os.path.join(subdir, 'qc_check_complete.json'))
            qc_check_passed = False
            qc_check = {}
            with open(os.path.join(subdir, 'qc_check_complete.json'), 'r') as f:
                qc_check = json.load(f)
            if re.match("PASS", qc_check.get("overall_pass_fail", "")):
                qc_check_passed = True
            upload_not_already_initiated = not os.path.exists(os.path.join(config['upload_staging_dir'], run_id))
            not_excluded = False
            if 'excluded_runs' in config:
                not_excluded = not run_id in config['excluded_runs']
            else:
                not_excluded = True

            conditions_checked = {
                "is_directory": subdir.is_dir(),
                "matches_illumina_run_id_format": ((matches_miseq_regex is not None) or
                                                   (matches_nextseq_regex is not None)),
                "upload_complete": upload_complete,
                "qc_check_complete": qc_check_complete,
                "qc_check_passed": qc_check_passed,
                "upload_not_already_initiated": upload_not_already_initiated,
                "not_excluded": not_excluded,
            }

            conditions_met = list(conditions_checked.values())
            run = {}
            if all(conditions_met):
                logging.info(json.dumps({"event_type": "run_directory_found", "sequencing_run_id": run_id, "run_directory_path": os.path.abspath(subdir.path)}))
                run['path'] = os.path.abspath(subdir.path)
                run['sequencing_run_id'] = run_id
                run['instrument_type'] = instrument_type
                yield run
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


def find_fastq(run, library_id, read_type):
    """
    Find the fastq file for a specific library on a specific run.

    :param run: Sequencing run. Keys: ['sequencing_run_id', 'path', 'instrument_type']
    :type run: dict[str, str]
    :param library_id: Library ID
    :type library_id: str
    :param read_type: Read type ('R1' or 'R2')
    :type read_type: str
    :return: Path to fastq file
    :rtype: Optional[str]
    """
    fastq_path = None
    run_path = run.get('path', None)
    run_id = run.get('sequencing_run_id', None)
    instrument_type = run.get('instrument_type', None)
    if instrument_type == 'miseq':
        # If the library ID contains underscores or periods, replace them with dashes.
        # Occasionally a library ID is entered into the SampleSheet with underscores or periods,
        # but the fastq files are named with dashes instead.
        library_id_cleaned = library_id.replace('_', '-').replace('.', '-')
        # New-style MiSeq runs have an 'Alignment' directory for each demultiplexing.
        if os.path.exists(os.path.join(run_path, 'Alignment_1')):
            alignment_dirs_glob = os.path.join(run_path, 'Alignment_*')
            alignment_dirs = list(sorted(glob.glob(alignment_dirs_glob), reverse=True))
            if len(alignment_dirs) > 0:
                last_alignment_dir = alignment_dirs[0]
                last_alignment_subdirs_glob = os.path.join(last_alignment_dir, '*')
                last_alignment_subdirs = list(sorted(glob.glob(last_alignment_subdirs_glob), reverse=True))
                if len(last_alignment_subdirs) > 0:
                    last_alignment_subdir = last_alignment_subdirs[0]
                    fastq_dir_path = os.path.join(last_alignment_subdir, 'Fastq')
                    fastq_glob = os.path.join(fastq_dir_path, library_id_cleaned +
                                              '_*' + read_type + '_*.fastq*')
                    fastq_paths = glob.glob(fastq_glob)
                    if len(fastq_paths) > 0:
                        fastq_path = fastq_paths[0]
                        logging.debug(json.dumps({"event_type": "fastq_file_found",
                                                  "sequencing_run_id": run_id,
                                                  "library_id": library_id,
                                                  "read_type": read_type,
                                                  "fastq_path": fastq_path}))
        # Old-style MiSeq runs store fastqs in 'Data/Intensities/BaseCalls' directory.
        elif os.path.exists(os.path.join(run_path, 'Data', 'Intensities', 'BaseCalls')):
            fastq_dir_path = os.path.join(run_path, 'Data', 'Intensities', 'BaseCalls')
            fastq_glob = os.path.join(fastq_dir_path, library_id_cleaned +
                                      '_*' + read_type + '_*.fastq*')
            fastq_paths = glob.glob(fastq_glob)
            if len(fastq_paths) > 0:
                fastq_path = fastq_paths[0]
                logging.debug(json.dumps({"event_type": "fastq_file_found",
                                          "sequencing_run_id": run_id,
                                          "library_id": library_id,
                                          "read_type": read_type,
                                          "fastq_path": fastq_path}))
        else:
            logging.error(json.dumps({"event_type": "fastq_file_not_found",
                                      "sequencing_run_id": run_id,
                                      "library_id": library_id,
                                      "read_type": read_type,
                                      "fastq_path": fastq_path}))
    elif instrument_type == 'nextseq':
        library_id_cleaned = library_id.replace('.', '-')
        if os.path.exists(os.path.join(run_path, 'Analysis')):
            analysis_dirs_glob = os.path.join(run_path, 'Analysis', '*')
            analysis_dirs = list(sorted(glob.glob(analysis_dirs_glob), reverse=True))
            if len(analysis_dirs) > 0:
                last_analysis_dir = analysis_dirs[0]
                fastq_dir_path = os.path.join(last_analysis_dir, 'Data', 'fastq')
                fastq_glob = os.path.join(fastq_dir_path, library_id_cleaned +
                                          '_*' + read_type + '_*.fastq*')
                fastq_paths = glob.glob(fastq_glob)
                if len(fastq_paths) > 0:
                    fastq_path = fastq_paths[0]
                    logging.debug(json.dumps({"event_type": "fastq_file_found",
                                              "sequencing_run_id": run_id,
                                              "library_id": library_id,
                                              "read_type": read_type,
                                              "fastq_path": fastq_path}))
        else:
            logging.error(json.dumps({"event_type": "fastq_file_not_found",
                                      "sequencing_run_id": run_id,
                                      "library_id": library_id,
                                      "read_type": read_type,
                                      "fastq_path": fastq_path}))

    return fastq_path


def prepare_downsampling_samplesheet(config, run):
    """
    Prepare a SampleSheet to use for downsampling.

    :param config: Application config. 
    :type config: dict[str, object]
    :param run: Sequencing run to prepare SampleSheet.csv file for. Keys: ['sequencing_run_id', 'path', 'instrument_type']
    :type run: dict[str, str]
    :return: Downsampling samplesheet [ { ID: str, R1: str, R2: str, GENOME_SIZE: str, COVERAGE: str } ]
    :rtype: list[dict[str, str]]
    """
    downsampling_samplesheet = []
    run_id = run['sequencing_run_id']
    run_dir = run['path']
    instrument_type = run['instrument_type']
    run_samplesheets = samplesheet.find_samplesheets(run_dir, instrument_type)
    samplesheet_to_parse = samplesheet.choose_samplesheet_to_parse(run_samplesheets, instrument_type, run_id)

    parsed_samplesheet = None
    if samplesheet_to_parse is not None and os.path.exists(samplesheet_to_parse):
        parsed_samplesheet = samplesheet.parse_samplesheet(samplesheet_to_parse, instrument_type)

    if parsed_samplesheet is not None and instrument_type == 'nextseq':
        if 'cloud_data' in parsed_samplesheet:
            for library in parsed_samplesheet['cloud_data']:
                if 'project_name' in library:
                    for project in config['projects']:
                        if library['project_name'] == project['local_project_id'] and project.get('downsample_reads', False):
                            local_project_id = project['local_project_id']
                            samplesheet_library = {}
                            library_id = library['sample_id']
                            library_r1_fastq = find_fastq(run, library_id, 'R1')
                            library_r2_fastq = find_fastq(run, library_id, 'R2')
                            samplesheet_library['ID'] = library_id
                            samplesheet_library['R1'] = library_r1_fastq
                            samplesheet_library['R2'] = library_r2_fastq
                            try:
                                samplesheet_library['GENOME_SIZE'] = str(project['genome_size_mb']) + 'm'
                            except KeyError as e:
                                samplesheet_library['GENOME_SIZE'] = None
                            samplesheet_library['COVERAGE'] = str(project.get('max_depth', None))
                            downsampling_samplesheet.append(samplesheet_library)

    elif parsed_samplesheet is not None and instrument_type == 'miseq':
        if 'data' in parsed_samplesheet:
            for library in parsed_samplesheet['data']:
                if 'sample_project' in library:
                    for project in config['projects']:
                        if library['sample_project'] == project['local_project_id']:
                            local_project_id = project['local_project_id']
                            samplesheet_library = {}
                            # There are two fields that can be used for a sample ID in the SampleSheet:
                            # Sample_ID and Sample_Name. Sometimes Sample_ID is populated with
                            # 'S1', 'S2', 'S3', etc. instead of the actual sample ID.
                            # If that's the case, take the 'Sample_Name' from the SampleSheet.
                            if re.match('S\d+$', library['sample_id']) and 'sample_name' in library:
                                library_id = library['sample_name']
                            else:
                                library_id = library['sample_id']
                            library_r1_fastq = find_fastq(run, library_id, 'R1')
                            library_r2_fastq = find_fastq(run, library_id, 'R2')
                            samplesheet_library['ID'] = library_id
                            samplesheet_library['R1'] = library_r1_fastq
                            samplesheet_library['R2'] = library_r2_fastq
                            try:
                                samplesheet_library['GENOME_SIZE'] = str(project['genome_size_mb']) + 'm'
                            except KeyError as e:
                                samplesheet_library['GENOME_SIZE'] = None
                            samplesheet_library['COVERAGE'] = str(project.get('max_depth', None))
                            downsampling_samplesheet.append(samplesheet_library)

    return downsampling_samplesheet


def downsample_reads(config, run_id, samplesheet):
    """
    Downsample reads.

    :param config: Application config.
    :type config: dict[str, object]
    :param run_id: Sequencing run ID.
    :type run_id: str
    :param samplesheet: samplesheet: [ { ID: str, R1: str, R2: str , GENOME_SIZE: str, COVERAGE: str} ]
    """
    downsampled_reads = {}

    
    output_dir = os.path.join(config['downsampling']['output_dir'], run_id)
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    downsampling_samplesheet_path = os.path.join(output_dir, '_'.join([run_id, 'samplesheet.csv']))
    with open(downsampling_samplesheet_path, 'w') as f:
        writer = csv.DictWriter(f, fieldnames=['ID', 'R1', 'R2', 'GENOME_SIZE', 'COVERAGE'])
        writer.writeheader()
        for library in samplesheet:
            writer.writerow(library)
    logging.info(json.dumps({"event_type": "downsampling_samplesheet_created", "sequencing_run_id": run_id, "samplesheet_path": downsampling_samplesheet_path}))

    if len(samplesheet) < 1:
        return downsampled_reads

    pipeline_name = config['downsampling']['pipeline_name']
    pipeline_version = config['downsampling']['pipeline_version']
    analysis_log_path = os.path.abspath(os.path.join(output_dir, '_'.join([run_id, 'analysis.log'])))
    analysis_report_path = os.path.abspath(os.path.join(output_dir, '_'.join([run_id, 'nextflow_report.html'])))
    analysis_trace_path = os.path.abspath(os.path.join(output_dir, '_'.join([run_id, 'nextflow_trace.tsv'])))
    analysis_timeline_path = os.path.abspath(os.path.join(output_dir, '_'.join([run_id, 'nextflow_timeline.html'])))
    work_dir = config['downsampling']['work_dir']
    downsampling_command = [
        'nextflow',
        '-log', analysis_log_path,
        'run',
        pipeline_name,
        '-r',
        pipeline_version,
        '-profile', 'conda',
        '--cache', os.path.join(os.path.expanduser('~'), '.conda/envs'),
        '--samplesheet_input', downsampling_samplesheet_path,
        '--outdir', output_dir,
        '-work-dir', work_dir,
        '-with-report', analysis_report_path,
        '-with-trace', analysis_trace_path,
        '-with-timeline', analysis_timeline_path,
    ]
    logging.info(json.dumps({"event_type": "downsampling_started", "sequencing_run_id": run_id, output_dir: output_dir}))
    analysis_complete = {"timestamp_analysis_start": datetime.datetime.now().isoformat()}
    try:
        analysis_result = subprocess.run(downsampling_command, capture_output=True, check=True)
        analysis_complete['timestamp_analysis_complete'] = datetime.datetime.now().isoformat()
        analysis_complete['analysis_success'] = True
        with open(os.path.join(output_dir, 'analysis_complete.json'), 'w') as f:
            json.dump(analysis_complete, f, indent=2)
        for library in samplesheet:
            library_id = library['ID']
            max_depth = library['COVERAGE']
            fastq_path_r1 = os.path.join(output_dir, library_id, library_id + '-downsample-' + str(max_depth) + 'x_R1.fastq.gz')
            fastq_path_r2 = os.path.join(output_dir, library_id, library_id + '-downsample-' + str(max_depth) + 'x_R2.fastq.gz')
            if os.path.exists(fastq_path_r1) and os.path.exists(fastq_path_r2):
                downsampled_reads[library_id] = {}
                downsampled_reads[library_id]['library_id'] = library_id
                downsampled_reads[library_id]['fastq_path_r1'] = fastq_path_r1
                downsampled_reads[library_id]['fastq_path_r2'] = fastq_path_r2
        logging.info(json.dumps({"event_type": "downsampling_complete", "sequencing_run_id": run_id}))
    except subprocess.CalledProcessError as e:
        logging.error(json.dumps({"event_type": "downsampling_failed", "sequencing_run_id": run_id, "command": " ".join(downsampling_command)}))
        
    return downsampled_reads


def prepare_samplelist(config, run, downsampled_reads={}):
    """
    Prepare a SampleList for a specific run.

    :param config: Application config.
    :type config: dict[str, object]
    :param run: Sequencing run to prepare SampleList.csv file for. Keys: ['sequencing_run_id', 'path', 'instrument_type']
    :type run: dict[str, str]
    :return: List of samples to upload. Keys: ['Sample_Name', 'Project_ID', 'File_Forward', 'File_Forward_Absolute_Path', 'File_Reverse', ''File_Reverse_Absolute_Path']
    :rtype: list[dict[str, str]]
    """
    samples_to_upload = []

    run_id = run['sequencing_run_id']
    run_dir = run['path']
    instrument_type = run['instrument_type']

    run_samplesheets = samplesheet.find_samplesheets(run_dir, instrument_type)
    samplesheet_to_parse = samplesheet.choose_samplesheet_to_parse(run_samplesheets, instrument_type, run_id)

    parsed_samplesheet = None
    if samplesheet_to_parse is not None and os.path.exists(samplesheet_to_parse):
        parsed_samplesheet = samplesheet.parse_samplesheet(samplesheet_to_parse, instrument_type)

    if parsed_samplesheet is not None and instrument_type == 'nextseq':
        if 'cloud_data' in parsed_samplesheet:
            for library in parsed_samplesheet['cloud_data']:
                if 'project_name' in library:
                    for project in config['projects']:
                        if library['project_name'] == project['local_project_id']:
                            samplelist_library = {}
                            samplelist_library['Sample_Name'] = library['sample_id']
                            samplelist_library['Project_ID'] = project['remote_project_id']
                            samplelist_library['Project_Name'] = project['remote_project_name']
                            samplelist_library['Project_ID_Local'] = project['local_project_id']
                            samplelist_library['Project_Name_Local'] = project['local_project_name']
                            samples_to_upload.append(samplelist_library)
    elif parsed_samplesheet is not None and instrument_type == 'miseq':
        if 'data' in parsed_samplesheet:
            for library in parsed_samplesheet['data']:
                if 'sample_project' in library:
                    for project in config['projects']:
                        if library['sample_project'] == project['local_project_id']:
                            samplelist_library = {}
                            # There are two fields that can be used for a sample ID in the SampleSheet:
                            # Sample_ID and Sample_Name. Sometimes Sample_ID is populated with
                            # 'S1', 'S2', 'S3', etc. instead of the actual sample ID.
                            # If that's the case, take the 'Sample_Name' from the SampleSheet.
                            if re.match('S\d+$', library['sample_id']) and 'sample_name' in library:
                                samplelist_library['Sample_Name'] = library['sample_name']
                            else:
                                samplelist_library['Sample_Name'] = library['sample_id']
                            samplelist_library['Project_ID'] = project['remote_project_id']
                            samplelist_library['Project_Name'] = project['remote_project_name']
                            samplelist_library['Project_ID_Local'] = project['local_project_id']
                            samplelist_library['Project_Name_Local'] = project['local_project_name']
                            samples_to_upload.append(samplelist_library)

    for sample in samples_to_upload:
        library_id = sample['Sample_Name']
        fastq_absolute_path_r1 = None
        if library_id in downsampled_reads:
            fastq_absolute_path_r1 = downsampled_reads[library_id]['fastq_path_r1']
        else:
            fastq_absolute_path_r1 = find_fastq(run, sample['Sample_Name'], 'R1')
        if fastq_absolute_path_r1 is not None:
            sample['File_Forward_Absolute_Path'] = fastq_absolute_path_r1
            fastq_filename_r1 = os.path.basename(fastq_absolute_path_r1)
            sample['File_Forward'] = fastq_filename_r1
        else:
            logging.error(json.dumps({"event_type": "find_fastq_failed",
                                      "sequencing_run_id": run_id,
                                      "sequencing_run_directory": run_dir,
                                      "library_id": sample['Sample_Name'],
                                      "local_project_id": sample['Project_ID_Local'],
                                      "remote_project_id": sample['Project_ID']}))
        fastq_absolute_path_r2 = None
        if library_id in downsampled_reads:
            fastq_absolute_path_r2 = downsampled_reads[library_id]['fastq_path_r2']
        else:
            fastq_absolute_path_r2 = find_fastq(run, sample['Sample_Name'], 'R2')
        if fastq_absolute_path_r2 is not None:
            sample['File_Reverse_Absolute_Path'] = fastq_absolute_path_r2
            fastq_filename_r2 = os.path.basename(fastq_absolute_path_r2)
            sample['File_Reverse'] = fastq_filename_r2

    return samples_to_upload


def prepare_upload_dir(config, run, sample_list):
    """
    Prepare upload directory for run.

    :param config:
    :type config: dict[str, object]
    :param run: Sequencing run to prepare upload directory for. Keys: ['sequencing_run_id', 'path', 'instrument_type']
    :type run: dict[str, str]
    :param sample_list: List of samples to upload. Keys: ['Sample_Name', 'Project_ID', 'File_Forward', 'File_Forward_Absolute_Path', 'File_Reverse', ''File_Reverse_Absolute_Path']
    :type sample_list: list[dict[str, str]]
    :return: Upload dir path
    :rtype: path
    """
    run_id = run['sequencing_run_id']
    upload_staging_dir = config['upload_staging_dir']
    run_upload_staging_dir = os.path.join(upload_staging_dir, run_id)

    sample_list_headers = [
        'Sample_Name',
        'Project_ID',
        'File_Forward',
        'File_Reverse',
    ]

    if not os.path.exists(run_upload_staging_dir):
        os.makedirs(run_upload_staging_dir)
        sample_list_path = os.path.join(run_upload_staging_dir, 'SampleList.csv')
        
        with open(sample_list_path, 'w') as f:
            f.write('[Data]\n')
            f.write(','.join(sample_list_headers) + '\n')
            for sample in sample_list:
                all_headers_present = all([header in sample.keys() for header in sample_list_headers])
                if not all_headers_present:
                    logging.error(json.dumps({"event_type": "samplelist_headers_missing", "sequencing_run_id": run_id, "library_id": sample.get("Sample_Name", "unknown"), "expected_headers": sample_list_headers, "actual_keys": list(sample.keys()), "missing_headers": [header for header in sample_list_headers if header not in sample.keys()]}))
                    continue
                f.write(','.join([sample[k] for k in sample_list_headers]) + '\n')
                symlink_src_fwd = sample['File_Forward_Absolute_Path']
                symlink_dest_fwd = os.path.join(run_upload_staging_dir, sample['File_Forward'])
                symlink_src_rev = sample['File_Reverse_Absolute_Path']
                symlink_dest_rev = os.path.join(run_upload_staging_dir, sample['File_Reverse'])
                os.symlink(symlink_src_fwd, symlink_dest_fwd)
                os.symlink(symlink_src_rev, symlink_dest_rev)

        upload_prepared_path = os.path.join(run_upload_staging_dir, 'upload_prepared.json')
        libraries = []
        for sample in sample_list:
            all_headers_present = all([header in sample.keys() for header in sample_list_headers])
            if not all_headers_present:
                logging.error(json.dumps({"event_type": "samplelist_headers_missing", "sequencing_run_id": run_id, "library_id": sample.get("Sample_Name", "unknown"), "expected_headers": sample_list_headers, "actual_headers": list(sample.keys()), "missing_headers": [header for header in sample_list_headers if header not in sample.keys()]}))
                continue
            fastq_forward_path = os.path.join(run_upload_staging_dir, sample['File_Forward'])
            fastq_forward_realpath = os.path.realpath(fastq_forward_path)
            fastq_forward_md5 = None
            with open(fastq_forward_realpath, 'rb') as f:
                fastq_forward_hash = hashlib.md5()
                while chunk := f.read(8192):
                    fastq_forward_hash.update(chunk)
                fastq_forward_md5 = fastq_forward_hash.hexdigest()
            fastq_reverse_path = os.path.join(run_upload_staging_dir, sample['File_Reverse'])
            fastq_reverse_realpath = os.path.realpath(fastq_reverse_path)
            fastq_reverse_md5 = None
            with open(fastq_reverse_realpath, 'rb') as f:
                fastq_reverse_hash = hashlib.md5()
                while chunk := f.read(8192):
                    fastq_reverse_hash.update(chunk)
                fastq_reverse_md5 = fastq_reverse_hash.hexdigest()
            
            library = {}
            library['library_id'] = sample['Sample_Name']
            library['local_project_id'] = sample['Project_ID_Local']
            library['local_project_name'] = sample['Project_Name_Local']
            library['remote_project_id'] = sample['Project_ID']
            library['remote_project_name'] = sample['Project_Name']
            library['fastq_forward_md5'] = fastq_forward_md5
            library['fastq_reverse_md5'] = fastq_reverse_md5
            libraries.append(library)

        if len(libraries) < 1:
            return None

        upload_prepared = {
            'sequencing_run_id': run_id,
            'libraries': libraries,
        }
        with open(upload_prepared_path, 'w') as f:
            json.dump(upload_prepared, f, indent=2)
            f.write("\n")

        return run_upload_staging_dir
    else:
        return None


def upload_run(config, run, upload_dir):
    """
    Initiate an analysis on one directory of fastq files.
    """
    run_id = run['sequencing_run_id']
    upload_successful = False
    upload_id = str(uuid.uuid4())

    upload_url = '/'.join([
        config['container_url'],
        upload_id,
        config['sas_token'],
    ])

    azcopy_command = [
        'azcopy',
        'cp',
        '--put-md5',
        '--recursive',
        '--follow-symlinks',
        '--output-type', 'json',
        '--from-to=LocalBlob',
        '--metadata=upload_id=' + upload_id,
        '--exclude-pattern=*NML_Upload_Finished*',
        upload_dir,
        upload_url,        
    ]

    logging.info(json.dumps({"event_type": "upload_started", "sequencing_run_id": run_id, "azcopy_command": " ".join(azcopy_command)}))
    timestamp_upload_started = datetime.datetime.now().isoformat()
    upload_started_path = os.path.join(upload_dir, 'upload_started.json')
    upload_started_file_contents = {
        'timestamp_upload_started': timestamp_upload_started,
    }
    with open(upload_started_path, 'w', encoding='utf-8') as f:
        json.dump(upload_started_file_contents, f, indent=2)
        f.write("\n")
    try:
        upload_successful = False
        data_upload_result = subprocess.run(azcopy_command, capture_output=True, check=True, text=True)
        for line in data_upload_result.stdout.splitlines():
            line_json = json.loads(line)
            if "MessageContent" in line_json and line_json["MessageContent"].startswith("{"):
                parsed_message_content = json.loads(line_json["MessageContent"])
                line_json["MessageContent"] = parsed_message_content
            logging.info(json.dumps(line_json))
        if data_upload_result.returncode == 0:
            upload_successful = True
            upload_completed_timestamp = datetime.datetime.now().isoformat()
            upload_completed_path = os.path.join(upload_dir, 'upload_completed.json')
            upload_completed_file_contents = {
                'timestamp_upload_completed': upload_completed_timestamp,
            }
            with open(upload_completed_path, 'w', encoding='utf-8') as f:
                json.dump(upload_completed_file_contents, f, indent=2)
                f.write("\n")
        logging.info(json.dumps({"event_type": "upload_completed", "sequencing_run_id": run_id, "azcopy_command": " ".join(azcopy_command)}))
        time.sleep(5)
        
    except subprocess.CalledProcessError as e:
        logging.error(json.dumps({"event_type": "upload_failed", "sequencing_run_id": run_id, "azcopy_command": " ".join(azcopy_command)}))

    if upload_successful:
        nml_upload_finished_file_contents = {"action": "UPLOAD", "result": upload_successful, "job_id": upload_id}
        nml_upload_finished_filename = upload_id + "-NML_Upload_Finished.json"
        nml_upload_finished_path = os.path.join(upload_dir, nml_upload_finished_filename)
        with open(nml_upload_finished_path, "w", encoding="utf-8") as f:
            json.dump(nml_upload_finished_file_contents, f, indent=2)
            f.write("\n")

        upload_url = config['container_url'] + config['sas_token']

        azcopy_command = [
            'azcopy',
            'cp',
            '--output-type', 'json',
            '--from-to=LocalBlob',
            nml_upload_finished_path,
            upload_url,
        ]

        try:
            nml_upload_finished_file_result = subprocess.run(azcopy_command, capture_output=True, check=True, text=True)
            for line in nml_upload_finished_file_result.stdout.splitlines():
                line_json = json.loads(line)
                if "MessageContent" in line_json and line_json["MessageContent"].startswith("{"):
                    parsed_message_content = json.loads(line_json["MessageContent"])
                    line_json["MessageContent"] = parsed_message_content
                logging.info(json.dumps(line_json))
            if nml_upload_finished_file_result.returncode == 0:
                logging.info(json.dumps({"event_type": "upload_confirmation_completed", "sequencing_run_id": run_id, "azcopy_command": " ".join(azcopy_command)}))
        except subprocess.CalledProcessError as e:
            logging.error(json.dumps({"event_type": "upload_confirmation_failed", "sequencing_run_id": run_id, "azcopy_command": " ".join(azcopy_command)}))
