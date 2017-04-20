#!/bin/sh

echo "Running tests..."
echo
echo "This could take a while. Only errors are shown. If the output"
echo "is empty, all tests passed."
echo
python3 run_tests.py | grep -B1 '!!!'
