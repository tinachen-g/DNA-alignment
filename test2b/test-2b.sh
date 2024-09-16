#!/usr/bin/env bash

# To run the test, make sure that your script (`2b.py`) 
# and the files needed (`test-2b-1.fasta`, `test-2b-2.fasta`, `test-2b-3.fasta` and `score_matrix.json`) are in the working directory.
# Or modify the commands below to accomodate for their path

# In the terminal, execute `bash test-2b.sh` to run all the test cases.
# Correct output of these test cases are written in the folder `outputs`.

### TIPS: Since multiple optimal alignment exists for some of the test cases, we recommend comparing
### each provided output file manually with your output, rather than using commands like `sdiff`.
### especially for test 1, 4 and 5.


# TEST 1
# PLEASE NOTE THAT THIS TEST CASE HAS MULTIPLE OPTIMAL ALIGNMENTS, PLEASE CHECK THE OUTPUT FILE FOR ALL RESULTS
python 2b.py -f test2b/test-2b-1.fasta -s test2b/score_matrix.json -d 50
echo -e "----------------------------------------------\n"
# The result of this test case is in `test-1.out`

# Test 2
python 2b.py -f test2b/test-2b-1.fasta -s test2b/score_matrix.json -d 100
echo -e "----------------------------------------------\n"
# The result of this test case is in `test-2.out`

# TEST 3
python 2b.py -f test2b/test-2b-1.fasta -s test2b/score_matrix.json -d 150
echo -e "----------------------------------------------\n"
# The result of this test case is in `test-3.out`

# TEST 4
# PLEASE NOTE THAT THIS TEST CASE HAS MULTIPLE OPTIMAL ALIGNMENTS, PLEASE CHECK THE OUTPUT FILE FOR ALL RESULTS
python 2b.py -f test2b/test-2b-2.fasta -s test2b/score_matrix.json -d 100
echo -e "----------------------------------------------\n"
# The result of this test case is in `test-4.out`

# TEST 5
# PLEASE NOTE THAT THIS TEST CASE HAS MULTIPLE OPTIMAL ALIGNMENTS, PLEASE CHECK THE OUTPUT FILE FOR ALL RESULTS
# FYI: This fasta file (`test-2b-3.fasta`) contains the actual gene APOE in mouse and human,
# This gene is involved in making a protein that helps carry cholesterol and other types of fat in the bloodstream,
# whose genetic variant has been shown to be related to Alzheimer's disease!
# Read more at: https://www.sciencedirect.com/topics/medicine-and-dentistry/apolipoprotein-e
# This test case will not be included in autograder.
python 2b.py -f test2b/test-2b-3.fasta -s test2b/score_matrix.json -d 100
echo -e "----------------------------------------------\n"
# The result of this test case is in `test-5.out`