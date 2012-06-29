#!/bin/bash

search='Proj4'
replace='ProjFour'
list=$(find . -name \*.php)
for i in $list; do
 sed -i "s/${search}/${replace}/g" $i
done
