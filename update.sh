#!/bin/bash

search='atanTwo'
replace='atan2'
list=$(find . -name \*.php)
for i in $list; do
 sed -i "s/${search}/${replace}/g" $i
done
