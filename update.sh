#!/bin/bash

search='ProjFourphp::$common->tsfnz'
replace='ProjFourphp_Common::tsfnz'
list=$(find . -name \*.php)
for i in $list; do
 sed -i "s/${search}/${replace}/g" $i
done
