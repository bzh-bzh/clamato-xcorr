#!/usr/bin/zsh

source ~/.zshrc
cd /global/u1/b/bzh/clamato-xcorr
conda activate clamato-xcorr

for i in cfg/*.yaml; do
    [ -f "$i" ] || break
    echo "$i"
    python3 BootCovar_script.py $i
done