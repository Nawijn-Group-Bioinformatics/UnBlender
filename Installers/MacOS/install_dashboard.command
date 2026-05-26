#!/bin/zsh

# ======= Install R through pixi =========
set -e

cd "$(dirname "$0")"

if ! which pixi &> /dev/null; then
    echo "🔍 Pixi not found. Installing now..."
    curl -fsSL https://pixi.sh/install.sh | sh
    echo "✅ Pixi installation complete."
    $(dirname "$0")/manifest.sh
else
    echo "🎉 Pixi is already installed. Skipping installation."
    $(dirname "$0")/manifest.sh
fi
