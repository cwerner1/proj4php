#!/bin/bash

search='Proj4php_common::'
replace='Proj4php_Common::'
list=$(find . -name \*.php)
for i in $list; do
 sed -i "s/${search}/${replace}/g" $i
done
