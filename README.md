# auto-irida-upload
Automated upload of sequence data to the IRIDA platform.

# Installation

# Usage
Start the tool as follows:

```bash
auto-irida-upload --config config.json
```

See the Configuration section of this document for details on preparing a configuration file.

More detailed logs can be produced by controlling the log level using the `--log-level` flag:

```bash
auto-irida-upload --config config.json --log-level debug
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
  "upload_tmpdir": "/path/to/auto-irida-upload-tmp",
  "run_parent_dirs": [
    "/path/to/instrument_1",
    "/path/to/instrument_2"
  ]
}
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
