#!/bin/sh

# Find the Java binary.
JAVA=$(which java)

# Find the absolute path to this script, and change to that directory.
SCRIPT_NAME=$(echo "$0" | sed "s/.*\\///g")
if [ "$0" = "./$SCRIPT_NAME" ]; then
  ABS_PATH=$(echo `pwd`/$0 | sed "s/\/$SCRIPT_NAME//g" | sed "s/\/\.//g")
else
  ABS_PATH=$(echo "$0" | sed "s/\/$SCRIPT_NAME//g")
fi
if [ "$ABS_PATH" != "$SCRIPT_NAME" ]; then
  cd "$ABS_PATH"
fi

# Use shell expansion to translate special characters (e.g. ~ --> /home/USER).
ARGS=""
for VAR in "$@"; do
  ARGS=$(printf "$ARGS $VAR")
done

# Run Ecotype Simulation.
$JAVA -jar ecosim.jar $ARGS
