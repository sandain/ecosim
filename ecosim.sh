#!/usr/bin/env sh

JAVA=$(command -v java)
if [ -z "$JAVA" ]; then
  echo "ecosim.sh: java was not found." >&2
  exit 1
fi

# Run Ecotype Simulation.
exec "$JAVA" -jar "$(dirname "$0")/ecosim.jar" "$@"
