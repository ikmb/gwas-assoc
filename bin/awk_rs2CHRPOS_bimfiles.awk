#!/usr/bin/env gawk

BEGIN {
	file1 = ARGV[1]

	# uncompressed bim file
	cmd="cat " file1
	while (cmd | getline) {
	        print $1"\t"$1":"$4"\t"$3"\t"$4"\t"$5"\t"$6
	}

	exit (0)
}
