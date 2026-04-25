#!/usr/bin/env python3

# Test parsing with the correct field mapping
line = 'SMEL4.1_03g019210.1.01 -            203 SSXT                 PF05030.17    62   2.2e-22   79.3   4.9   1   2   6.9e-26   2.2e-22   79.3   4.9     3    60    25    82    24    84 0.95 gene=SMEL4.1_03g019210.1 Name:"Similar to GIF2 GRF1-interacting factor 2 (Arabidopsis thaliana OX=3702)"'

parts = line.strip().split()

print('HMMER domtblout format:')
print('0: target_name')
print('1: accession')  
print('2: tlen')
print('3: query_name')
print('4: query_accession') 
print('5: qlen')
print('6: E-value')
print('7: score')
print('8: bias')
print('9: # (domain number)')
print('10: of (total domains)')
print('11: c-Evalue')
print('12: i-Evalue')
print('13: domain score')
print('14: domain bias')
print('15: from (hmm coord)')
print('16: to (hmm coord)')
print('17: from (ali coord)')
print('18: to (ali coord)')
print('19: from (env coord)')
print('20: to (env coord)')
print('21: acc (accuracy)')

print(f'\nActual fields for this line:')
for i, field in enumerate(parts[:22]):
    print(f'{i}: {field}')

print(f'\nAccuracy should be at index 21: {parts[21]}')