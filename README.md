[![Unit Tests](https://github.com/BCCDC-PHL/auto-irida-azure-upload/actions/workflows/unit_tests.yml/badge.svg)](https://github.com/BCCDC-PHL/auto-irida-azure-upload/actions/workflows/unit_tests.yml)

# auto-irida-azure-upload
Automated upload of sequence data to the IRIDA platform.

# Installation
This tool assumes that the system has `azcopy` on the `PATH`. Follow the directions [here](https://learn.microsoft.com/en-us/azure/storage/common/storage-use-azcopy-v10) to
download and install azcopy.

It's recommended to install the tool into an isolated environment using a tool such as conda, mamba or virtualenv:
```
conda create -n auto-irida-azure-upload python=3 pip
conda activate auto-irida-azure-upload
```

The tool is pip-installable. Using the `-e` flag will allow updates/hotfixes to be applied directly to the codebase and automatically installed.

```
git clone https://github.com/BCCDC-PHL/auto-irida-azure-upload.git
cd auto-irida-azure-upload
pip install -e .
```

# Usage
Start the tool as follows:

```bash
auto-irida-azure-upload --config config.json
```

See the Configuration section of this document for details on preparing a configuration file.

More detailed logs can be produced by controlling the log level using the `--log-level` flag:

```bash
auto-irida-azure-upload --config config.json --log-level debug
```

# Configuration
This tool takes a single config file, in JSON format, with the following structure:

```json
{
  "container_url": "",
  "sas_token": "",
  "projects_definition_file": "/path/to/projects.csv",
  "excluded_runs_list": "/path/to/excluded_runs.csv",
  "scan_interval_seconds": 3600,
  "upload_staging_dir": "/path/to/auto-irida-upload-tmp",
  "run_parent_dirs": [
    "/path/to/instrument_1",
    "/path/to/instrument_2"
  ]
}
```

The `projects_definition_file` should be .csv format and should include the following fields:

```
local_project_id
local_project_name
remote_project_id
remote_project_name
```

## Downsampling

If downsampling is needed prior to upload, add the following to the `config.json` file:

```json
{
  "downsampling": {
    "enabled": true,
    "output_dir": "/path/to/downsampled-reads",
    "work_dir": "/path/to/downsampling-work",
    "pipeline_name": "BCCDC-PHL/downsample-reads",
    "pipeline_version": "v0.2.0"
  }
}
```

...and the following fields to the `projects.csv` file:

```
downsample_reads
genome_size_mb
max_depth
```

# Application Flowchart

```mermaid
main('__main__.main')
main --> scan('core.scan')
```

# Logging
This tool outputs [structured logs](https://www.honeycomb.io/blog/structured-logging-and-your-team/) in [JSON Lines](https://jsonlines.org/) format:

Every log line should include the fields:

- `timestamp`
- `level`
- `module`
- `function_name`
- `line_num`
- `message`

...and the contents of the `message` key will be a JSON object that includes at `event_type`. The remaining keys inside the `message` will vary by event type.

```json
{"timestamp": "2022-09-22T11:32:52.287", "level": "INFO", "module", "core", "function_name": "scan", "line_num", 56, "message": {"event_type": "scan_start"}}
```
