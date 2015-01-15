#!/bin/sh
cd "$2"
if [ -f outfile ]; then
  rm outfile
fi
if [ -f outtree ]; then
  rm outtree
fi
if [ -f /etc/debian_version ]; then
  echo 'Y' | phylip "$1" > /dev/null
else
  echo 'Y' | "$1" > /dev/null
fi
