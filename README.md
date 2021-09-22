# CS412-Gene-Mutation-Detection

# Project Description
A ‘motif’ is a pattern in a sequence. For example, in DNA sequences (which are sequences over the alphabet {A,C,G,T}), an example of a motif is the pattern ‘TCACGTG’. A slightly more complex motif is the pattern TC[A/C]CGTG, which represents ‘either TCACGTG or TCCCGTG’. An occurrence of a motif in a given DNA sequence is called a ‘site’. When gene mutation or human planted gene sequences exist, there might impact the local motif table, so that it is possible to detect the planted gene locations. The task of this project is to detect the locations and sequences of planted gene sequences.

In this project, effects of three parameters on the developed algorithm need to be evaluated:
- Information Content per Column (ICPC) = 1, 1.5, and 2
- Motif Length (ML) = 6, 7, and 8
- Sequence Count (SC) = 5, 10, and 20

# Work

Step 1: Build gene sequences with equal probability for each alphabet in {A,C,T,G} at any locations as benchmark. Randomly generate a position weight matrix (motif) . Planted gene will be generated according to this motif table and then placed in the benchmark sequences randomly.

Step 2: Develop an algorithm to find the location and sequences of those planted DNA sequences.

Step 3: Evaluate the performance of developed algorithm.

# Algorithm
Gibbs sampling algorithm was implemented to tackle the challenge. To avoid local optimum and make results approach the global optimum, phase shifting was used. 

- Step 1: Choose random starting points of motifs for each sequence and record the starting points in a list named as “pst”
![image](https://user-images.githubusercontent.com/47155713/134284192-99c12e2f-a521-438b-9eab-22832e87a13b.png)
- Step 2: Choose random sequence from a set (for example sequence 1)  to search the most possible motif
![image](https://user-images.githubusercontent.com/47155713/134284217-54f12572-3500-4da6-9bfc-a864f58b920a.png)
- Step 3: Formulate position weight matrix (PWM) of width of ML by excluding chosen sequence in step 2
- Step 4: Compute the weights for each location in the sequence chosen in Step 2 using Ax = Qx/Px and normalize weights as the probability for corresponding locations
![image](https://user-images.githubusercontent.com/47155713/134284329-7f089e1b-f86a-4838-8e9a-cb1d4842010c.png)
- Step 5: Randomly sample a new location of starting point for the selected sequence in Step 2 based on the normalized weight and renew the list “pst”.
![image](https://user-images.githubusercontent.com/47155713/134284371-f625a99d-3f76-4d80-85cc-32be8669a351.png)
- Step 6: Phase Shifts: after 40 iterations, the location of starting points for each sequence will be moved a distance of m to left and right, in which m ranges from 1 to 5. Used to jump out of local optimum
![image](https://user-images.githubusercontent.com/47155713/134284398-4d57c21f-6e40-4ccc-8de3-e3cbf93953a4.png)
- Step 7 : Repeat Step1~6 to until convergence:
    Convergence Criteria: 
      1. The SC latest results are the same.
      2. Number of iteration < 20000.
      3. Number of iteration ≥ 4000 + added iteration
        Added iteration depends on the IC of result: 1000*(2*ML-Current IC)^2 
        low IC tends to have more iterations
- Step 8: Do one more phase shift comparison after iterations. Avoid the situation, which is induced by background DNA noise, that the predicted location overlapping with the actual location at the same error
![image](https://user-images.githubusercontent.com/47155713/134284600-51586104-cfe2-4da1-9e63-5849aeda3412.png)





# Code Description
- "run_this.py": main function.
- "Step1_2.py": functions for creating the datasets.
- "Step2_Function.py": functions for implementing Gibbs sampling algorithm and phase shifting.
- "Step3_Function.py": functions for evaluation.

# Performance
Detailed information can be found in the project report, "Project Report.pdf", in the repo.
![alt text](https://github.com/Zhongyihhh/CS412-Gene-Mutation-Detection/blob/main/image/Screen%20Shot%202021-09-08%20at%2018.32.21.png)
