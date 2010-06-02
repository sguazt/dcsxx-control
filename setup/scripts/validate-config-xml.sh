#!/bin/sh

epacsim_path=$PWD
xmldoc=$1

xmllint --noout --schema  $epacsim_path/setup/conf/epacsim_config.xsd $xmldoc
