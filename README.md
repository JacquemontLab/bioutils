[![Jacquemont's Lab Header](labheader.png)](https://www.jacquemont-lab.org/)

# bioutils

A collection of commonly used bioinformatics Bash scripts developed and maintained by our lab to streamline data processing and project management tasks.

## Scripts Overview

* **summarize**
  Parses command-line arguments and uses DuckDB to compute summary statistics (mean, median, quantiles) for each column of a Parquet file.

* **show**
  Displays the schema (column names and data types) of a Parquet file using DuckDB.

* **pcat**
  Prints a specified number of rows from a Parquet file in TSV format to STDOUT with optional verbose output showing the executed SQL query using DuckDB.

* **parquet2csv**
  Converts a Parquet file to a tab-delimited CSV file with header using DuckDB.

* **open\_access**
  Recursively sets default ACL permissions (`rwx`) for a specified user on a directory, granting read, write, and execute access to all current and future files.

* **merge\_multithread**
  Merges multiple files from a directory into a single output file by processing them in parallel chunks, preserving only the header from the first file.

* **markdown\_to\_pdf**
  Converts a Markdown file to PDF by rendering it to HTML with grip and converting the HTML to PDF using wkhtmltopdf.

* **launch\_salloc**
  Submits an interactive cluster job request with specified duration, CPUs, and memory per CPU using `salloc`.

* **generate\_project\_directory**
  Interactively creates a standardized project directory structure with subfolders for data, scripts, results, logs, and software, initializes a README.md with placeholders and current date, and sets user/group permissions.

* **duckdb**
  DuckDB executable.
