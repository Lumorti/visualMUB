#!/bin/bash

mkdir -p data
toTry="
        2,2 
        2,3
        2,1,1 
        2,2,1 
        2,2,2
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
        6,3,1,1
        6,6,1,1
        6,6,2,1
        6,6,3,1
        6,6,3,2
        6,3,3,2
        6,3,3,3
        6,6,6,6
        6,1,1,1,1,1,1
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
        16,6,6,6
        18,6,6,6
        20,6,6,6
        20,3,3,3
        20,10,2,2
        20,1,1,1,1,1,1,1,1,1,1
        50,6,6,6
        50,1,1,1,1,1,1,1,1,1,1
        100,6,6,6
        100,1,1,1,1,1,1,1,1,1,1
      "
rm -f data/temp.log
for i in $toTry; do
    withDashes=$(echo $i | sed 's/,/-/g')
    if [ -f data/$withDashes.log ]; then
        echo "Already have $i"
        continue
    fi
    echo "Running $i"
    ./run -N $i -v --shotgun | tee data/temp.log
    mv data/temp.log data/$withDashes.log
done
