#!/usr/bin/env bash
# Install feynsage into SageMath, with FORM for Dirac traces.
# Works on macOS (Apple silicon and Intel), Linux (x86_64 and arm64) and Windows through WSL.
#
#   ./install.sh            install (and install FORM / SageMath if they are missing)
#   ./install.sh --test     install, then run the test suite
#   ./install.sh --no-deps  install feynsage only, do not touch FORM or SageMath
set -euo pipefail

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
RUN_TESTS=0; DEPS=1
for a in "$@"; do
  case "$a" in
    --test) RUN_TESTS=1 ;;
    --no-deps) DEPS=0 ;;
    -h|--help) sed -n '2,9p' "$0"; exit 0 ;;
  esac
done

say()  { printf '\033[1;35m==>\033[0m %s\n' "$*"; }
warn() { printf '\033[1;33mwarning:\033[0m %s\n' "$*"; }
die()  { printf '\033[1;31merror:\033[0m %s\n' "$*"; exit 1; }

OS="$(uname -s)"; ARCH="$(uname -m)"
case "$OS" in
  Darwin) PLATFORM=mac ;;
  Linux)  PLATFORM=linux; grep -qi microsoft /proc/version 2>/dev/null && say "running inside WSL" ;;
  MINGW*|MSYS*|CYGWIN*) die "native Windows is not supported by SageMath; install WSL (wsl --install) and run this script there" ;;
  *) die "unknown system $OS" ;;
esac
say "system: $PLATFORM ($ARCH)"

# ------------------------------------------------------------------ SageMath
find_sage() {
  if command -v sage >/dev/null 2>&1; then command -v sage; return; fi
  for c in "$HOME"/miniforge3/envs/sage/bin/sage "$HOME"/mambaforge/envs/sage/bin/sage "$HOME"/miniconda3/envs/sage/bin/sage \
           "$HOME"/miniforge3/envs/*/bin/sage "$HOME"/mambaforge/envs/*/bin/sage "$HOME"/miniconda3/envs/*/bin/sage \
           "$HOME"/anaconda3/envs/*/bin/sage /Applications/SageMath-*.app/Contents/Frameworks/Sage.framework/Versions/*/venv/bin/sage; do
    [ -x "$c" ] && { echo "$c"; return; }
  done
  return 1
}

SAGE="$(find_sage || true)"
if [ -z "$SAGE" ]; then
  [ "$DEPS" = 1 ] || die "SageMath not found (and --no-deps was given)"
  if [ "$PLATFORM" = mac ]; then
    command -v brew >/dev/null || die "SageMath not found. Install Homebrew (https://brew.sh) or the SageMath app from https://github.com/3-manifolds/Sage_macOS/releases, then rerun"
    say "installing the SageMath app with Homebrew (large download)"
    brew install --cask sage
  else
    CONDA="$(command -v mamba || command -v conda || true)"
    [ -n "$CONDA" ] || die "SageMath not found. Install Miniforge (https://github.com/conda-forge/miniforge), then rerun; the script will create a 'sage' environment"
    say "creating the conda environment 'sage' from conda-forge (large download)"
    "$CONDA" create -y -n sage -c conda-forge sage
  fi
  SAGE="$(find_sage)" || die "SageMath was installed but the sage command is not on the PATH; open a new terminal and rerun"
fi
say "SageMath: $("$SAGE" --version 2>/dev/null | head -1) ($SAGE)"

# ------------------------------------------------------------------ FORM
if command -v form >/dev/null 2>&1; then
  say "FORM: $(command -v form)"
elif [ "$DEPS" = 1 ]; then
  if [ "$PLATFORM" = mac ] && command -v brew >/dev/null; then
    say "installing FORM with Homebrew"; brew install form
  else
    case "$PLATFORM-$ARCH" in
      linux-x86_64)  asset=x86_64-linux ;;
      linux-aarch64|linux-arm64) asset=arm64-linux ;;
      mac-arm64)     asset=arm64-osx ;;
      mac-x86_64)    asset=x86_64-osx ;;
    esac
    tag="$(curl -fsSL https://api.github.com/repos/form-dev/form/releases/latest | sed -n 's/.*"tag_name": *"\(v[^"]*\)".*/\1/p' | head -1)"
    [ -n "$tag" ] || die "could not reach GitHub to download FORM"
    ver="${tag#v}"
    tmp="$(mktemp -d)"
    say "downloading FORM $ver ($asset) into ~/.local/bin"
    curl -fsSL "https://github.com/form-dev/form/releases/download/$tag/form-$ver-$asset.tar.gz" | tar -xz -C "$tmp"
    mkdir -p "$HOME/.local/bin"
    find "$tmp" -type f \( -name form -o -name tform \) -exec cp {} "$HOME/.local/bin/" \;
    chmod +x "$HOME/.local/bin/form" "$HOME/.local/bin/tform" 2>/dev/null || true
    rm -rf "$tmp"
    "$HOME/.local/bin/form" -v >/dev/null 2>&1 || warn "the downloaded FORM does not run on this machine; build it from source (https://github.com/form-dev/form)"
    case ":$PATH:" in *":$HOME/.local/bin:"*) ;; *) warn "add ~/.local/bin to your PATH so that feynsage finds FORM" ;; esac
  fi
else
  warn "FORM not found: form.trace and form.dirac_trace will not work until it is installed"
fi

# ------------------------------------------------------------------ C compiler (optional, for the fast reducer)
if command -v cc >/dev/null 2>&1 || command -v gcc >/dev/null 2>&1 || command -v clang >/dev/null 2>&1; then
  say "C compiler found: the finite-field reducer will use its compiled kernel"
else
  warn "no C compiler: reduce_ff falls back to pure Python (same answers, about ten times slower)"
  [ "$PLATFORM" = mac ] && warn "install one with: xcode-select --install"
fi

# ------------------------------------------------------------------ feynsage itself
say "installing feynsage into Sage's Python"
PREFIX="$("$SAGE" -python -c 'import sys; print(sys.prefix)')"
case "$PREFIX" in
  *.app*|/var/tmp/sage-*|/usr/*|/opt/*) USERFLAG="--user" ;;   # never write into the app bundle or a system prefix
  *) USERFLAG="" ;;                                            # a conda env or a build you own
esac
"$SAGE" -pip install $USERFLAG --no-deps "$HERE"

say "checking the installation (this also compiles the C kernel)"
( cd "$HOME" && "$SAGE" -c "
import feynsage, feynsage.ff as ff, shutil
print('feynsage from', feynsage.__file__)
print('compiled kernel:', 'yes' if ff._clib() is not None else 'no (pure Python fallback)')
print('FORM:', feynsage.form.FORM if shutil.which('form') or __import__('os').path.exists(feynsage.form.FORM) else 'not found')
" )

if [ "$RUN_TESTS" = 1 ]; then
  say "running the tests"
  cd "$HERE"
  for t in tests/test_*.sage; do
    echo "--- $t"; "$SAGE" "$t" | tail -3
  done
fi
say "done. Try the notebook: examples/feynsage_walkthrough.ipynb (SageMath kernel)"
