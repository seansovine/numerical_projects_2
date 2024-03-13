#!/usr/bin/env bash

echo "Running clang-format using .clang-format config file."
# clang-format --dump-config -style=file

echo
echo "Changes will be:"
clang-format --dry-run -style=file -- src/**/*.{cpp,hpp}

echo
echo "Formatting:"
clang-format -i --verbose -style=file -- src/**/*.{cpp,hpp}
