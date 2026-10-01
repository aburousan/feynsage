#!/bin/bash
# feynsage against Kira 3.1 on the same machine, same targets, same seeds, same masters.
#     ./run_all.sh <sage> <kira> <fermat (fer64)> [case ...]
# Results in work/<case>_r<r>_s<s>/: *.json (wall time, CPU time, peak memory) and cmp_*.txt
# (feynsage against Kira, coefficient by coefficient).  summarize.py makes the tables.
SAGE=$1; KIRA=$2; export FERMATPATH=$3; shift 3
CASES=${@:-kite vertex sunset doublebox}
H=$(cd "$(dirname "$0")" && pwd); W=$H/work; mkdir -p "$W"
LIMIT=3600
for case in $CASES; do
  sizes=$(python3 -c "exec(open('$H/cases.py').read()); print(' '.join('%d:%d' % rs for rs in CASES['$case']['sizes']))")
  for rs in $sizes; do
    r=${rs%:*}; s=${rs#*:}; B=$W/${case}_r${r}_s${s}; mkdir -p "$B"
    for th in 1 8; do
      for tool in fs_ff fs_trimmed kira_fermat kira_firefly; do
        [ "$tool" = fs_trimmed ] && [ "$th" != 1 ] && continue          # the trimmed reducer has one thread
        [ -e "$W/${case}_skip_${tool}_${th}" ] && continue
        out=$B/${tool}_${th}.json
        case $tool in
          fs_ff)        cmd=("$SAGE" "$H/run_feynsage.sage" "$case" "$r" "$s" "$th" ff "$B/fs_ff_$th") ;;
          fs_trimmed)   cmd=("$SAGE" "$H/run_feynsage.sage" "$case" "$r" "$s" 1 trimmed "$B/fs_trimmed_1") ;;
          kira_fermat)  cmd=(python3 "$H/run_kira.py" "$case" "$r" "$s" "$th" fermat "$B" "$KIRA") ;;
          kira_firefly) cmd=(python3 "$H/run_kira.py" "$case" "$r" "$s" "$th" firefly "$B" "$KIRA") ;;
        esac
        echo "== $(date +%H:%M:%S) $case r=$r s=$s $tool threads=$th"
        python3 "$H/measure.py" "$out" "$LIMIT" "${cmd[@]}"
        grep -q '"rc": "timeout"' "$out" && touch "$W/${case}_skip_${tool}_${th}" && echo "   timeout: skipping larger sizes"
        if [ "${tool#kira}" != "$tool" ] && [ -e "$B/fs_ff_1/feynsage.txt" ]; then
          "$SAGE" "$H/compare.sage" "$case" "$B/fs_ff_1" "$B/kira_${tool#kira_}_$th" | tee "$B/cmp_${tool}_${th}.txt"
        fi
      done
    done
  done
done
echo "ALL DONE"
