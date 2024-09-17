#!/usr/bin/env bash

# To run the test, make sure that your script (`2c.py`) 
# and the files needed (`test-2c-1.fasta`, `test-2c-2.fasta`, `test-2c-3.fasta` and `score_matrix.json`) are in the working directory.
# Or modify the commands below to accomodate for their path

# In the terminal, execute `bash test-2c.sh` to run all the test cases or execute each line in your terminal separately.
# Correct output of these test cases are written in the folder `test-2c-output`.

### TIPS: Since multiple optimal alignment exists for some of the test cases, we recommend comparing
### each provided output file manually with your output, rather than using commands like `sdiff`.
### especially for test 4 and 5.

# Test 1
python 2c.py -f test2c/test-2c-1.fasta -s test2c/score_matrix.json -d 430 -e 30
# The result of this test case is in `test-1.out`
echo -e "----------------------------------------------\n"

# Test 2
python 2c.py -f test2c/test-2c-1.fasta -s test2c/score_matrix.json -d 430 -e 60
# The result of this test case is in `test-2.out`
echo -e "----------------------------------------------\n"

# Test 3
python 2c.py -f test2c/test-2c-1.fasta -s test2c/score_matrix.json -d 430 -e 430
# The result of this test case is in `test-3.out`
echo -e "----------------------------------------------\n"

# Test 4
# PLEASE NOTE THAT THIS TEST CASE HAS MULTIPLE OPTIMAL ALIGNMENTS, PLEASE CHECK THE OUTPUT FILE FOR ALL RESULTS
python 2c.py -f test2c/test-2c-2.fasta -s test2c/score_matrix.json -d 200 -e 50
# The result of this test case is in `test-4.out`
echo -e "----------------------------------------------\n"

# Test 5
# PLEASE NOTE THAT THIS TEST CASE HAS MULTIPLE OPTIMAL ALIGNMENTS, PLEASE CHECK THE OUTPUT FILE FOR ALL RESULTS
python 2c.py -f test2c/test-2c-2.fasta -s test2c/score_matrix.json -d 200 -e 200
# The result of this test case is in `test-5.out`
echo -e "----------------------------------------------\n"

# Test 6
python 2c.py -f test2c/test-2c-3.fasta -s test2c/score_matrix.json -d 430 -e 30
# The result of this test case is in `test-6.out`
echo -e "----------------------------------------------\n"

# Test 7
python 2c.py -f test2c/test-2c-3.fasta -s test2c/score_matrix.json -d 100 -e 100
# The result of this test case is in `test-7.out`
echo -e "----------------------------------------------\n"

# Secret Test
python 2c.py -f sequences.fasta -s score_matrix.json -d 430 -e 30