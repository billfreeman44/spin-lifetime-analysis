# Spin lifetime analysis

**Copyright 2014 Evan Sosenko**

## Requirements

- [Python 3](http://www.python.org/)
  with [pip](http://www.pip-installer.org/).

## Setup

Install the required Python packages with

````bash
$ pip install -r requirements.txt
````

Depending on how your system is configured:

* You may need to run this commands with `sudo`.
* You may want to install the core SciPy packages via your package manager.
* You may want to use [virtualenv](http://www.virtualenv.org/en/latest/).

## Usage

Execute the analysis with

````bash
$ make
````

Output will be saved in the `build` directory.

Python output will be redirected to `stdout.log`
and errors to `stderr.log`.

Remove the build and logs with

````bash
$ make clean
````

Serve the build directory from a local http server on port `8000` with

````bash
$ make serve
````
