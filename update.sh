#!/bin/bash

search='qs1'
replace='qsOne'
list=$(find . -name \*.php)
for i in $list; do
 sed -i "s/${search}/${replace}/g" $i
done
