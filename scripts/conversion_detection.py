#!/usr/bin/env python3

import sys

def map_files(file1_path, file2_path, output_path):
    file1_dict = {}
    with open(file1_path, 'r', encoding='utf-8') as f1:
        for idx, line in enumerate(f1, 1):
            parts = line.strip().split()
            if len(parts) == 2:
                letter = parts[0]
                number = int(parts[1])
                file1_dict[letter] = (idx, number, letter)

    with open(file2_path, 'r', encoding='utf-8') as f2:
        letters = [line.strip() for line in f2 if line.strip()]

    valid_entries = []
    for letter in letters:
        if letter in file1_dict:
            idx, num, orig_letter = file1_dict[letter]
            valid_entries.append((idx, num, orig_letter))

    valid_entries.sort(key=lambda x: x[0])

    num_groups = []
    letter_groups = []
    if valid_entries:
        current_num_group = [valid_entries[0][1]]
        current_letter_group = [valid_entries[0][2]]
        for i in range(1, len(valid_entries)):
            prev_line = valid_entries[i-1][0]
            current_line = valid_entries[i][0]
            current_num = valid_entries[i][1]
            current_letter = valid_entries[i][2]
            
            if current_line == prev_line + 1:
                current_num_group.append(current_num)
                current_letter_group.append(current_letter)
            else:
                num_groups.append(current_num_group)
                letter_groups.append(current_letter_group)
                current_num_group = [current_num]
                current_letter_group = [current_letter]
        num_groups.append(current_num_group)
        letter_groups.append(current_letter_group)

    output = []
    for nums, letters in zip(num_groups, letter_groups):
        num_str = ','.join(map(str, nums))
        letter_str = ','.join(letters)
        output.append(f"{num_str} {letter_str}")

    with open(output_path, 'w', encoding='utf-8') as out_file:
        out_file.write('\n'.join(output))

if __name__ == '__main__':
    if len(sys.argv) != 4:
        print("usage: python script.py file1 file2 output_file")
        sys.exit(1)
    file1_path = sys.argv[1]
    file2_path = sys.argv[2]
    output_path = sys.argv[3]
    map_files(file1_path, file2_path, output_path)
    