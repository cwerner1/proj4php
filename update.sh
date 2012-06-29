#!/bin/bash

search='eZero'
replace='e0'
list=$(find . -name \*.php)
for i in $list; do
 sed -i "s/${search}/${replace}/g" $i
done
