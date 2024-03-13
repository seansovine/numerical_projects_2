#!/usr/bin/env bash

echo "Running clang-format using .clang-format config file."

changes=$(clang-format --dry-run -style=file -- src/**/*.{cpp,hpp} 2>&1)

echo
if [ ! -z "$changes" ]; then
	echo "Changes will be made to fix:"
	echo
	echo "$changes"
else
	echo "No changes made."
fi

echo
clang-format -i --verbose -style=file -- src/**/*.{cpp,hpp}
