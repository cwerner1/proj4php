#!/bin/bash

search='$pjdThreeParam'
replace='pjdThreeParam'
list=$(find . -name \*.php)
for i in $list; do
 sed -i "s/${search}/${replace}/g" $i
done
