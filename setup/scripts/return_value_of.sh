#!/bin/sh

$($* 2>&1 >/dev/null)

echo $?
