#!/bin/bash

for file in $(ls *.py)
do
  sed -i '' -e 's/\r$//' $file
done
