<?php

$xmlFile="./test/data/test-conf_1.xml";

$doc = new DOMDocument();
$doc->preserveWhiteSpace=false;
$doc->formatOutput = true;
$doc->load($xmlFile);


$Valid=$doc->schemaValidate("./setup/epacsim_config.xsd");  // $doc is your xml loaded file,  and  "client_test.xsd" your xml schema definition file

echo $Valid;   // Returns  1 (if valid)

?>
