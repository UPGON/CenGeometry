#!/bin/bash
# Double-click this file to open CenGeometry in your web browser.
# First run takes a couple of minutes while it installs; after that it is quick.

cd "$(dirname "$0")" || exit 1
clear
echo "==============================================="
echo "  CenGeometry"
echo "==============================================="
echo

PY=""
for c in python3 python; do
    if command -v "$c" >/dev/null 2>&1; then PY="$c"; break; fi
done
if [ -z "$PY" ]; then
    echo "Python was not found on this computer."
    echo
    echo "Install it from https://www.python.org/downloads/ then"
    echo "double-click this file again."
    echo
    read -r -p "Press Return to close. "
    exit 1
fi

if [ ! -d .venv ]; then
    echo "First-time setup — this happens once and takes a few minutes."
    echo
    "$PY" -m venv .venv || { echo "Could not create the workspace."; read -r -p "Press Return. "; exit 1; }
    ./.venv/bin/pip install --quiet --upgrade pip
    ./.venv/bin/pip install --quiet -r requirements.txt || {
        echo "Could not install the required packages."; read -r -p "Press Return. "; exit 1; }
    echo "Setup complete."
    echo
fi

echo "Starting… your browser will open automatically."
echo "Leave this window open while you work."
echo "To stop, close this window or press Control-C."
echo
./.venv/bin/streamlit run app.py
