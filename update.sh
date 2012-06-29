#!/bin/bash

search='Proj4php::$common->adjustLon'
replace="Proj4php_common::adjustLon"
list=$(find . -name \*.php)
for i in $list; do
 sed -i "s/${search}/${replace}/g" $i
done
