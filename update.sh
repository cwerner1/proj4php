#!/bin/bash

search="Proj4php::$common->adjustLon"
replace="Proj4php_common::adjustLon"
list=$(find . -name \*.php)
for i in $list; do
	sed "s/$search/$replace/" $i > _$i
	mv $i ${i#_}
done
