#!/bin/bash
# Map+quant chunk submission driver: keeps at most 4 of your jobs in the system
# (site AssocMaxJobsLimit ~5). Runs in a persistent screen session.
# Safe to stop and re-run: driver_submitted.txt prevents duplicates;
# chunk scripts are idempotent (finished outputs are skipped).
LOG=~/september/logs/driver.log
SUBMITTED=~/september/logs/driver_submitted.txt
LIST=~/september/list_sample.txt
MAXJOBS=4
touch "$SUBMITTED"

log() { echo "[$(date '+%F %T')] $*" | tee -a "$LOG"; }
is_submitted() { grep -qx "$1" "$SUBMITTED" 2>/dev/null; }

queue_remaining() {
  # print "M:s-e" / "Q:s-e" entries for chunks that still have missing outputs
  # and have not been submitted in this pass
  local s e n line proj sid type strand missing entry
  s=1; while [ $s -le 407 ]; do
    e=$((s+24)); [ $e -gt 407 ] && e=407; missing=0
    for n in $(seq $s $e); do
      line=$(sed -n "${n}p" "$LIST")
      [[ -z "$line" || "$line" =~ ^# ]] && continue
      IFS=$'\t' read -r proj sid type strand <<< "$line"
      [[ -s ~/september/map/$proj/$sid/${sid}Aligned.sortedByCoord.out.bam ]] || { missing=1; break; }
    done
    entry="M:$s-$e"
    { [ $missing -eq 1 ] && ! is_submitted "$entry"; } && echo "$entry"
    s=$((e+1))
  done
  s=1; while [ $s -le 407 ]; do
    e=$((s+49)); [ $e -gt 407 ] && e=407; missing=0
    for n in $(seq $s $e); do
      line=$(sed -n "${n}p" "$LIST")
      [[ -z "$line" || "$line" =~ ^# ]] && continue
      IFS=$'\t' read -r proj sid type strand <<< "$line"
      [[ -s ~/september/quant/$proj/$sid/quant.sf ]] || { missing=1; break; }
    done
    entry="Q:$s-$e"
    { [ $missing -eq 1 ] && ! is_submitted "$entry"; } && echo "$entry"
    s=$((e+1))
  done
}

PASS=1
while [ $PASS -le 4 ]; do
  mapfile -t QUEUE < <(queue_remaining)
  log "PASS $PASS: ${#QUEUE[@]} chunks to submit"
  [ ${#QUEUE[@]} -eq 0 ] && { log "DONE: all 407 map+quant outputs present"; exit 0; }
  idx=0
  while [ $idx -lt ${#QUEUE[@]} ]; do
    running=$(squeue -h -u "$USER" 2>/dev/null | wc -l)
    if [ "$running" -lt $MAXJOBS ]; then
      entry=${QUEUE[$idx]}; t=${entry%%:*}; r=${entry#*:}; s=${r%-*}; e=${r#*-}
      if [ "$t" = "M" ]; then CMD=~/september/scripts/map_chunk.sbatch; else CMD=~/september/scripts/quant_chunk.sbatch; fi
      if sbatch "$CMD" "$s" "$e" >>"$LOG" 2>&1; then log "submit $entry OK"; echo "$entry" >> "$SUBMITTED"; idx=$((idx+1)); else log "submit $entry REJECTED, waiting"; fi
    fi
    sleep 300
  done
  log "PASS $PASS all submitted, waiting for queue to drain..."
  while [ "$(squeue -h -u "$USER" 2>/dev/null | wc -l)" -gt 0 ]; do sleep 600; done
  : > "$SUBMITTED"   # all jobs in this pass finished -> reset duplicate markers
  PASS=$((PASS+1))
done
log "4-pass limit reached"