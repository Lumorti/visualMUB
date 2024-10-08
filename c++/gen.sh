#!/bin/bash

# The list of all the possible combinations of the parameters
toTry="
        2,2 
        2,3
        2,1,1 
        2,2,1 
        2,2,2
        2,1,1,1
        3,2
        3,3
        3,4
        3,1,1
        3,2,2
        3,3,3
        3,1,1,1
        3,2,1,1
        3,2,2,1
        3,2,2,2
        3,3,2,2
        3,3,3,2
        3,3,3,3
        3,1,1,1,1
        3,2,1,1,1
        3,3,1,1,1
        3,3,3,1,1
        3,2,2,1,1
        3,2,2,2,2
        3,3,2,2,2
        3,3,3,2,2
        3,3,3,3,2
        3,3,3,3,3
        4,2
        4,3
        4,4
        4,5
        4,1,1
        4,2,2
        4,3,3
        4,4,4
        4,1,1,1
        4,2,1,1
        4,2,2,1
        4,2,2,2
        4,3,2,2
        4,3,3,2
        4,3,3,3
        4,4,3,3
        4,4,4,3
        4,4,4,4
        4,1,1,1,1
        4,4,4,4,4
        4,1,1,1,1,1
        4,2,1,1,1,1
        4,2,2,1,1,1
        4,2,2,2,1,1
        4,2,2,2,2,1
        4,4,1,1,1,1
        4,4,4,1,1,1
        4,4,4,4,1,1
        4,3,3,2,1,1
        4,4,4,4,4,1
        4,4,4,3,3,3
        4,4,4,4,3,3
        4,4,4,4,4,3
        4,4,4,4,4,4
        5,2
        5,3
        5,4
        5,5
        5,6
        5,1,1
        5,2,2
        5,3,3
        5,4,4
        5,5,5
        5,1,1,1
        5,3,3,2
        5,3,3,3
        5,4,3,3
        5,4,4,3
        5,1,1,1,1
        5,1,1,1,1,1
        5,2,2,2,2,2
        5,3,3,2,2,2
        5,3,3,3,3,3
        5,4,4,3,3,3
        5,4,4,4,4,4
        5,5,5,4,4,4
        5,5,5,5,5,5
        5,1,1,1,1,1,1
        5,2,1,1,1,1,1
        5,2,2,2,2,2,2
        6,2
        6,3
        6,4
        6,5
        6,6
        6,7
        6,1,1
        6,2,2
        6,3,3
        6,4,4
        6,5,5
        6,6,6
        6,1,1,1
        6,2,2,2
        6,3,3,3
        6,4,4,4
        6,3,1,1
        6,6,1,1
        6,6,2,1
        6,6,3,1
        6,6,3,2
        6,6,4,1
        6,6,5,1
        6,3,3,2
        6,3,3,3
        6,4,3,3
        6,5,4,1
        6,5,4,4
        6,5,5,5
        6,6,6,6
        6,1,1,1,1
        6,2,2,2,2
        6,3,3,3,3
        6,4,4,4,4
        6,5,5,5,5
        6,6,6,6,6
        6,1,1,1,1,1
        6,2,2,2,2,2
        6,3,3,3,3,3
        6,4,4,4,4,4
        6,1,1,1,1,1,1
        6,2,2,2,2,2,2
        7,2
        7,3
        7,4
        7,5
        7,6
        7,7
        7,8
        7,1,1
        7,2,2
        7,7,7
        8,2
        8,3
        8,4
        8,5
        8,6
        8,7
        8,8
        8,9
        8,8,8
        8,1,1
        8,2,2
        9,9,9
        10,2
        10,3
        10,4
        10,5
        10,6
        10,7
        10,8
        10,9
        10,10
        10,11
        10,3,3,3
        10,4,4,4
        10,5,5,5
        10,10,10
        10,1,1,1,1,1,1,1,1,1,1
      "

# Set up data directory
mkdir -p data
rm -f data/temp*

# Define function to run
run () {

    # Check if we already have this one
    withDashes=$(echo $1 | sed 's/,/-/g')
    if [ -f data/$withDashes.log ]; then
        echo "Already have $1"
        return
    fi

    # Run the input
    echo "Running $1"
    ./run -N $1 -v --optimshotgun -p 5000 | tee data/temp$withDashes.log
    mv data/temp$withDashes.log data/$withDashes.log

}

# Run in parallel
export -f run
parallel -j 16 run ::: $toTry

# Clean up
rm -f data/temp*

# Push to github
git add .
git commit -m "Automatically pushed data files"
git push
