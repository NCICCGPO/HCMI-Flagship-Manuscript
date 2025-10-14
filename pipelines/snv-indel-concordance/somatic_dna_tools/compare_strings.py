#!/usr/bin/env python
# USAGE: compare_strings.py FILE STRING
# DESCRIPTION: Create an empty file if a string does not match the content of an empty file

import sys

def get_string(file):
    with open(file) as input:
        for line in input:
            string = line.rstrip()
    return string
    
def compare(file_string, string):
    if file_string == string:
        return True
    return False

#  =====================
#  Main
#  =====================
if __name__ == "__main__":
    file = sys.argv[1]
    string = sys.argv[2]
    file_string = get_string(file)
    result = compare(file_string, string)
    if result:
        sys.stdout.write('Matched: ' + string + '\n')
    else:
        sys.stdout.write('')

