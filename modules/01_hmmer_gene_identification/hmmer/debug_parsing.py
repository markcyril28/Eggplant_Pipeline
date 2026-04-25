#!/usr/bin/env python3

# Test parsing of a sample line
line = 'SMEL4.1_03g019210.1.01 -            203 SSXT                 PF05030.17    62   2.2e-22   79.3   4.9   1   2   6.9e-26   2.2e-22   79.3   4.9     3    60    25    82    24    84 0.95 gene=SMEL4.1_03g019210.1 Name:"Similar to GIF2 GRF1-interacting factor 2 (Arabidopsis thaliana OX=3702)"'

parts = line.strip().split()
main_fields = parts[:21]

print(f'Total parts: {len(parts)}')
print(f'Main fields length: {len(main_fields)}')
print(f'Field 20 (accuracy): {main_fields[20] if len(main_fields) > 20 else "MISSING"}')

print('\nFirst 22 fields:')
for i, field in enumerate(parts[:22]):
    print(f'{i}: {field}')

print(f'\nAccuracy value should be: {parts[20]}')