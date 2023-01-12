# biomedical-datamining

usage:
`python motif.py [file_path] > [output.txt]`

example:
`python motif.py project1_beta/Q1.txt > project1_beta/A1.txt`

目前僅實做 1 mutation base 的 brute force motif finding algorithm

final testing script
```
python motif.py final_testing_data/Q1.txt > final_testing_data/A1.txt
python motif.py final_testing_data/Q2.txt > final_testing_data/A2.txt
python motif.py final_testing_data/Q3.txt > final_testing_data/A3.txt
python motif.py final_testing_data/Q4.txt > final_testing_data/A4.txt
python motif.py final_testing_data/Q5.txt > final_testing_data/A5.txt
python motif.py final_testing_data/Q6.txt > final_testing_data/A6.txt
python motif.py final_testing_data/Q7.txt > final_testing_data/A7.txt
python motif.py final_testing_data/Q8.txt > final_testing_data/A8.txt
```

實驗結果於資料夾中 A[n].txt