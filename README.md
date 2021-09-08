# CS412-Gene-Mutation-Detection

# Project Description
A ‘motif’ is a pattern in a sequence. For example, in DNA sequences (which are sequences over the alphabet {A,C,G,T}), an example of a motif is the pattern ‘TCACGTG’. A slightly more complex motif is the pattern TC[A/C]CGTG, which represents ‘either TCACGTG or TCCCGTG’. An occurrence of a motif in a given DNA sequence is called a ‘site’. When gene mutation or human planted gene sequences exist, there might impact the local motif table, so that it is possible to detect the planted gene locations. The task of this project is to detect the locations and sequences of planted gene sequences.

In this project, effects of three variables were evaluated:
- Information Content per Column (ICPC) = 1, 1.5, and 2
- Motif Length (ML) = 6, 7, and 8
- Sequence Count (SC) = 5, 10, and 20

# Work

Step 1: Build gene sequences with equal probability for each alphabet in {A,C,T,G} at any locations as benchmark. Randomly generate a position weight matrix (motif) . Planted gene will be generated according to this motif table and then placed in the benchmark sequences randomly.

Step 2: Develop an algorithm to find the location and sequences of those planted DNA sequences.

Step 3: Evaluate the performance of developed algorithm.

# Algorithm
Gibbs sampling algorithm was implemented to tackle the challenge. To avoid local optimum and make results approach the global optimum, phase shifting was used. 

# Code Description
"run_this.py": main function.
"Step1_2.py": functions for creating the datasets.
"Step2_Function.py": functions for implementing Gibbs sampling algorithm and phase shifting.
"Step3_Function.py": functions for evaluation.

# Performance
![alt text](https://github.com/Zhongyihhh/CS412-Gene-Mutation-Detection/blob/main/image/Screen%20Shot%202021-09-08%20at%2018.32.21.png)
