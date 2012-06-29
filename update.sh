#!/bin/bash

search='a2'
replace='aTwo'
list=$(find . -name \*.php)
for i in $list; do
 sed -i "s/${search}/${replace}/g" $i
done
