#!/bin/bash

# Fail on any error
set -e

# Compile the C++ program
echo "Compiling figure4b-reflective.cpp..."
g++ -std=c++17 -O3 -I "$(brew --prefix eigen)/include/eigen3" -o figure4b-reflective figure4b-reflective.cpp
echo "Compilation successful."

# Run the C++ program
echo "Running simulation..."
./figure4b-reflective
echo "Simulation completed. Results saved to results_walls.txt."

# Run the Python script using Anaconda's Python directly
echo "Running Python script for visualization using Anaconda Python..."
~/opt/anaconda3/bin/python figure4b-reflective.py
echo "Python script completed. Check the animation and plots."

echo "Script completed successfully."
