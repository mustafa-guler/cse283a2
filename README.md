# cse283a2
assignment 2 cse283

## Question 2
All scripts were compiled into one file, `question2.py`.
It takes an argument to determine which part to run.

## 2(b) Theoretical Mass
```shell
$ python question2.py 2b -h

usage: question2.py 2b [-h] peptide

positional arguments:
  peptide     target peptide

optional arguments:
  -h, --help  show this help message and exit
```

This takes just the peptide sequence and prints out the 
theoretical mass.

```shell
$ python question2.py 2b PEPTIDE
799.3694399999999
```

## 2(c) Candidate Peptide Generation
Tolerance of 50ppm instead of 10ppm was used.

```shell
$ python question2.py 2c -h
usage: question2.py 2c [-h] fasta precursor_mass

positional arguments:
  fasta           path to fasta file to search from
  precursor_mass  molecular mass of precursor

optional arguments:
  -h, --help      show this help message and exit
```

The precursor mass from 2a was used, this prints a table.

```shell
$ python question2.py 2c CAH_BOVIN_fasta.txt 5433.78076171

  start index    end index     mass    ppm error  peptide
-------------  -----------  -------  -----------  ----------------------------------------------------
          134          186  5433.89      20.4955  QQPDGLAVVGVFLKVGDANPALQKVLDALDSIKTKGKSTDFPNFDPGSLLPN
          210          257  5433.9       22.2236  LKEPISVSSQQMLKFRTLNFNAEGEPELLMLANWRPAQPLKNRQVRG
```

## 2(d) Theoretical Spectra
This outputs a csv with just two fields: the m/z value and the ion type.
It writes it to the file `${PEPTIDE}_spectrum.csv`

```shell
$ python question2.py 2d -h
usage: question2.py 2d [-h] candidates [candidates ...]

positional arguments:
  candidates  each peptide sequence to generate to theoretical spectrum for

optional arguments:
  -h, --help  show this help message and exit
```

Running `python question2.py 2d QQPDGLAVVGVFLKVGDANPALQKVLDALDSIKTKGKSTDFPNFDPGSLLPN LKEPISVSSQQMLKFRTLNFNAEGEPELLMLANWRPAQPLKNRQVRG`
generates the spectra for both candidates.

## 2(e)
This requires the csvs from 2d, so run that first.

```shell
$ python question2.py 2e -h
usage: question2.py 2e [-h] [-E E]
                       experimental theoreticals [theoreticals ...]

positional arguments:
  experimental  path to mgf file containing experimental spectrum
  theoreticals  paths to theoretical spectra generated from 2d

optional arguments:
  -h, --help    show this help message and exit
  -E E          error tolerance in ppm, DEFAULT=0.5
```

This does not work for any generic MGF file, it's not too robust. It works
for the given MGF file though. Use below command to run with error tolerance
of .1. Also, .mz files are created for use with 2f.

```shell
$ python question2.py 2e CAH_test_01_scan1133.mgf QQPDGLAVVGVFLKVGDANPALQKVLDALDSIKTKGKSTDFPNFDPGSLLPN_spectrum.csv LKEPISVSSQQMLKFRTLNFNAEGEPELLMLANWRPAQPLKNRQVRG_spectrum.csv -E .1
QQPDGLAVVGVFLKVGDANPALQKVLDALDSIKTKGKSTDFPNFDPGSLLPN 0.192990919597
LKEPISVSSQQMLKFRTLNFNAEGEPELLMLANWRPAQPLKNRQVRG 0.231863744949
```

## 2(f)
This plots the experimental spectra for each candidate. The explained
peaks are highlighted in red. The output from 2e is required.

```shell
$ python question2.py 2f -h
usage: question2.py 2f [-h] experimental theoreticals [theoreticals ...]

positional arguments:
  experimental  path to mgf file containing experimental spectrum
  theoreticals  paths to .mz files from 2e

optional arguments:
  -h, --help    show this help message and exit
```

```shell
$ python question2.py 2f CAH_test_01_scan1133.mgf QQPDGLAVVGVFLKVGDANPALQKVLDALDSIKTKGKSTDFPNFDPGSLLPN_spectrum.csv.mz LKEPISVSSQQMLKFRTLNFNAEGEPELLMLANWRPAQPLKNRQVRG_spectrum.csv.mz 
```

This creates png files of the plots.


