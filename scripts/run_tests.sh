#!/bin/sh

base_path=.
test_src_path="$base_path/test/src/dcs/control"
test_bin_path="$base_path/build/test/dcs/control"

test_files=$(ls $test_src_path/*.cpp)

for f in $test_files; do
	f=$(basename $f .cpp)

	t="$test_bin_path/$f"

	if [ -x $t ]; then
		echo -n "--- $f ==> "
		out=$($t 2>&1 | grep -i 'failed test')

		echo $out
	fi
done
