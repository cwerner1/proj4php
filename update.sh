#!/bin/bash

search='Proj4php::$common->asinz'
replace='Proj4php_common::asinz'
list=$(find . -name \*.php)
for i in $list; do
 sed -i "s/${search}/${replace}/g" $i
done
