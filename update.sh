#!/bin/bash

search='n0'
replace='nZero'
list=$(find . -name \*.php)
for i in $list; do
 sed -i "s/${search}/${replace}/g" $i
done
