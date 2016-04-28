# Spin lifetime analysis

[![MIT License](https://img.shields.io/github/license/evansosenko/spin-lifetime-analysis.svg)](./LICENSE.txt)

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

## Viewing fits with Fitalyzer.

After generating fits with `make`, run `make serve` and load
[io.evansosenko.com/fitalyzer/?firebase=spin-lifetime&port=8000&path=/fitalyzer](http://io.evansosenko.com/fitalyzer/?firebase=spin-lifetime&port=8000&path=/fitalyzer).

Note: you must visit [https://localhost:8000](https://localhost:8000)
in your browser and accept the SSL certificate for this to work.

## License

This code is licensed under the MIT license
with the exception of any files under the path `data/PhysRevLett.105.167202`.

All files under the path `data/PhysRevLett.105.167202`
were compiled from data presented in:

> [Tunneling Spin Injection into Single Layer Graphene](http://link.aps.org/doi/10.1103/PhysRevLett.105.167202).
> Wei Han, K. Pi, K. M. McCreary, Yan Li, Jared J. I. Wong, Adrian G. Swartz, and Roland K. Kawakami, Phys. Rev. Lett. 105, 167202 (2010)

This data was used with permission and is available for download:
[Joint Laboratory for Spintronics Research, Department of Physics and Astronomy, University of California, Riverside](http://physics.ucr.edu/~kawakami/jlsrPublications.html).

## Warranty

This software is provided "as is" and without any express or
implied warranties, including, without limitation, the implied
warranties of merchantibility and fitness for a particular
purpose.
