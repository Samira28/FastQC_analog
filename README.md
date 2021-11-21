# FastQC Analog

## Functional 
Here is the code for executing tasks and commands by analogy with the FastQC program, namely, you can calculate the content of guanine and cytosine, construct a distribution by the number of reads containing different compositions of guanine and cytosine. You can also look at the distribution of quality and quantity (coverage) depending on the position of the nucleotide in the sequence, display the number of reads depending on the insertion distance, look at the frequencies of various substitutions, i.e. how many times did the replacement A -> G, A -> C, A-> T, etc. for each letter and find out the average number of errors in reads.

* GC-composition
The program has functionality for building the distribution of the GC-composition for all reads. See which reads have the most percentage of guanine and cytosine.

* Quality distribution
Using FASTQ file, construct the error probability distribution depending on the position of the nucleotide.

* Coverage distribution
Using this function, you can build the number of reads by positions to adhere to, i.e. coverage of reads.
 

* Distribution of insertion distance
It is possible to show the number of reads depending on the insertion distance.

* Frequencies of various substitutions
It is also possible to calculate which substitutions and in what quantity were allowed.

* Average number of errors in reads
The average percentage of errors in reading is also provided by calculating the percentage of errors in each read, dividing it by the number of reads.


## Detailed instructions for starting and using the program

First of all, install all needed dependencies, provided below and in the requirements.txt. You can download the file to your own computer and open it with Jupyter Notebook. To obtain FastQC like results you should provide program with your own path to input file and directory where all results will be saved.

## Dependencies, OS and python versions on which the program was tested.


macOS Catalina 10.15.7

pip version == 21.0.1.
Python == Python 3.8.3

```
import pysam
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from Bio import SeqIO
import seaborn as sns
import math
```
