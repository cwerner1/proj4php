#!/bin/bash

search='cosP20'
replace='cosPTwenty'
list=$(find . -name \*.php)
for i in $list; do
 sed -i "s/${search}/${replace}/g" $i
done
