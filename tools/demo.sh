#!/usr/bin/env bash
# End-to-end demonstration of the paired-result publish guarantees against a
# running service (default http://localhost:8000).
#
# Scenarios:
#   1. normal job succeeds and both ends can be downloaded;
#   2. controlled R2 publish failure -> failed/publish_failed, the staged R1
#      is rolled back, no download is offered, the job never stays "running";
#   3. slow publish -> status queries stay responsive, a stop is accepted and
#      the job only reports "stopped" after the publish action exits and the
#      temporary files are cleaned (it can never flip back to succeeded).
#
# Usage: tools/demo.sh [base_url]
set -euo pipefail

BASE="${1:-http://localhost:8000}"
ROOT="$(cd "$(dirname "$0")/.." && pwd)"
TMPDIR="$(mktemp -d)"
trap 'rm -rf "$TMPDIR"' EXIT

py() { python3 -c "$1"; }
field() { py "import json,sys;print(json.load(sys.stdin)$1)"; }

wait_status() { # wait_status <job_id> <status1,status2,...> [timeout_s]
  local id="$1" wanted="$2" timeout="${3:-30}"
  local deadline=$(( $(date +%s) + timeout ))
  while (( $(date +%s) < deadline )); do
    local body status
    body="$(curl -sf "$BASE/jobs/$id")"
    status="$(field '["status"]' <<<"$body")"
    if [[ ",$wanted," == *",$status,"* ]]; then
      echo "$body"
      return 0
    fi
    sleep 0.1
  done
  echo "timeout waiting for job $id to reach $wanted" >&2
  return 1
}

submit() { # submit <r1> <r2> [form fields...]
  local r1="$1" r2="$2"; shift 2
  local args=(-X POST "$BASE/jobs" -F "r1=@$r1" -F "r2=@$r2" -F "min_quality=20")
  for extra in "$@"; do args+=(-F "$extra"); done
  curl -sf "${args[@]}"
}

echo "== 0. docs / health"
echo "docs HTTP: $(curl -s -o /dev/null -w '%{http_code}' "$BASE/docs")"
curl -sf "$BASE/health"; echo

echo
echo "== 1. normal success (mixed sample: 7 total / 4 kept / 3 discarded)"
JOB="$(submit "$ROOT/samples/mixed_R1.fastq.gz" "$ROOT/samples/mixed_R2.fastq.gz")"
ID="$(field '["job_id"]' <<<"$JOB")"
FINAL="$(wait_status "$ID" "succeeded,failed,stopped")"
echo "$FINAL" | py "import json,sys;d=json.load(sys.stdin);print('status =',d['status'],'| kept =',d['kept_pairs'],'| downloads =',d['downloads'])"
curl -sf -o "$TMPDIR/r1.gz" "$BASE/jobs/$ID/download/r1"
curl -sf -o "$TMPDIR/r2.gz" "$BASE/jobs/$ID/download/r2"
py "import gzip;print('downloaded R1 records =', gzip.open('$TMPDIR/r1.gz').read().count(b'\n')//4)"

echo
echo "== 2. controlled R2 publish failure (fail_r2_publish=true)"
JOB="$(submit "$ROOT/samples/good_R1.fastq.gz" "$ROOT/samples/good_R2.fastq.gz" "fail_r2_publish=true")"
ID="$(field '["job_id"]' <<<"$JOB")"
FINAL="$(wait_status "$ID" "succeeded,failed,stopped")"
echo "$FINAL" | py "import json,sys;d=json.load(sys.stdin);print('status =',d['status'],'| error =',d['error'],'| downloads =',d['downloads'])"
for side in r1 r2; do
  echo "download $side HTTP: $(curl -s -o /dev/null -w '%{http_code}' "$BASE/jobs/$ID/download/$side") (expect 409)"
done
echo "stop on failed job HTTP: $(curl -s -o /dev/null -w '%{http_code}' -X POST "$BASE/jobs/$ID/stop") (expect 409)"

echo
echo "== 3. slow publish (slow_publish_seconds=20): responsive queries + stop"
JOB="$(submit "$ROOT/samples/good_R1.fastq.gz" "$ROOT/samples/good_R2.fastq.gz" "slow_publish_seconds=20")"
ID="$(field '["job_id"]' <<<"$JOB")"

# Wait until the worker enters the publishing phase.
for _ in $(seq 1 100); do
  PHASE="$(curl -sf "$BASE/jobs/$ID" | field '.get("phase")')"
  [[ "$PHASE" == "publishing" ]] && break
  sleep 0.1
done
echo "observed phase: $PHASE"

# Measure query latency while the publish is held.
LAT="$(for _ in 1 2 3; do
  START="$(date +%s.%N)"
  curl -sf -o /dev/null "$BASE/jobs/$ID"
  END="$(date +%s.%N)"
  py "print(f'{($END-$START)*1000:.1f}')"
done | py "import sys;v=[float(x) for x in sys.stdin];print('query latency ms: min=%.1f max=%.1f'%(min(v),max(v)))")"
echo "$LAT"

STOP="$(curl -sf -X POST "$BASE/jobs/$ID/stop")"
echo "stop accepted, status = $(field '["status"]' <<<"$STOP") (expect stopping)"
START="$(date +%s.%N)"
FINAL="$(wait_status "$ID" "stopped,succeeded,failed")"
END="$(date +%s.%N)"
echo "$FINAL" | py "import json,sys;d=json.load(sys.stdin);print('final status =',d['status'],'| downloads =',d['downloads'])"
py "print('settled in %.2fs after stop (interrupted the 20s hold)'%($END-$START))"
for side in r1 r2; do
  echo "download $side HTTP: $(curl -s -o /dev/null -w '%{http_code}' "$BASE/jobs/$ID/download/$side") (expect 409)"
done

echo
echo "== demo complete"
