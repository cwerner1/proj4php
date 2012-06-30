#!/bin/bash

search='ep2'
replace='epTwo'
list=$(find . -name \*.php)
for i in $list; do
 sed -i "s/${search}/${replace}/g" $i
done
