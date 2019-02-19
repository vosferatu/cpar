# CPAR
Working with the Performance API (PAPI) -> http://icl.cs.utk.edu/papi/

Python lib: https://github.com/flozz/pypapi

## CPAR Project 1
#### Performance evaluation of a single core

In this project we will study the effect on the processor performance of the memory hierarchy when accessing large amounts of data. The product of two matrices will be used for this study.

Use the Performance API (PAPI) to collect relevant performance indicators of the program execution.

### 1. 

Download the example file from moodle that contains the basic algorithm in C/C++ that multiplies two  matrices,  i.e.  multiplies one  line  of  the  first  matrix  by  each  column  of  the second  matrix. 

Implement  the  same  algorithm  in another  programming language,  such  as JAVA, C#, Fortran, etc.

Register the  processing  time  for  the several versions  of  the  algorithm,  for  input  matrices from 600x600 to 3000x3000 elements with increments in both dimensions of 400.

### 2.

Implement   a   version   that   multiplies   an   element   from   the   first   matrix   by   the correspondent line of the second matrix. 

Register the processing time for the two versions of the algorithm, for input matrices from 600x600 to 3000x3000 elements with increments in both dimensions of 400. Register the processing time from 4000x4000to 10000x10000 with intervals of 2000.


### 3. 

Implement a block oriented algorithm that divides the matrices in blocks and uses the same sequence of computation as in 2.

Register  the  processing  time  from  4000x4000  to  10000x10000  with  intervals  of  2000, for different block sizes (e.g. 128, 256, 512).

### OUTCOMES

Write a report of up to 6 pages explaining the algorithm versions and analyzing the results obtained.

Justify  the  performance  parameters  selected  and  use  them  to  evaluate and compare the versions implemented. 

### To be delivered on: 22/03/2019 (Friday)

#### Parameters for Report Evaluation:

- Problem description and algorithms explanation;
- Performance metrics and evaluation methodology;
- Results and analysis;
- Conclusions;
- Writing quality.
