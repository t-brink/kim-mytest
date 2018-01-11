#!/bin/bash
set -e

while read line; do
    echo
    echo
    echo "Running tests:"
    echo
    for item in $line; do
        echo "     $item"
    done
    echo
    echo "This could take a while. Only errors are shown. If the output"
    echo "is empty, all tests passed."
    echo
    time python3 run_tests.py --only-fast $line | { grep -B1 '!!!' || true; }
done < run_test_input.txt
