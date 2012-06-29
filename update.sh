#!/bin/bash

search='rh0'
replace='rhZero'
list=$(find . -name \*.php)
for i in $list; do
 sed -i "s/${search}/${replace}/g" $i
done
