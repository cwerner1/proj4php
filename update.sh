#!/bin/bash

search='ProjFourphp_common::$maxIter'
replace='ProjFourphp_Common::$maxIter'
list=$(find . -name \*.php)
for i in $list; do
 sed -i "s/${search}/${replace}/g" $i
done
