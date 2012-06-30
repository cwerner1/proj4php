#!/bin/bash

search='to_meter'
replace='toMeter'
list=$(find . -name \*.php)
for i in $list; do
 sed -i "s/${search}/${replace}/g" $i
done
