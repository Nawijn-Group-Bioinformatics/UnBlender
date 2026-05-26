#!/bin/bash
# ======= Install R through pixi =========
set -e
cd "$(dirname "$0")"

if ! which pixi &> /dev/null; then
    echo "🔍 Pixi not found. Installing now..."
    curl -fsSL https://pixi.sh/install.sh | sh
    echo "✅ Pixi installation complete."
    bash manifest.sh
else
    echo "🎉 Pixi is already installed. Skipping installation."
    bash manifest.sh
fi
