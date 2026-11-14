#!/usr/bin/env bash
# End-to-end smoke test against a running service (default :8000).
# Usage: tools/smoke.sh [base_url]
set -euo pipefail

BASE="${1:-http://localhost:8000}"
ROOT="$(cd "$(dirname "$0")/.." && pwd)"

echo "# docs:  $(curl -sf -o /dev/null -w '%{http_code}' "$BASE/docs")"
echo "# health: $(curl -sf "$BASE/health")"

submit() { # submit <r1> <r2> <threshold>
  curl -sf -X POST "$BASE/jobs" \
    -F "r1=@$1" -F "r2=@$2" -F "min_quality=$3"
}

wait_job() { # wait_job <job_id>
  for _ in $(seq 1 100); do
    body="$(curl -sf "$BASE/jobs/$1")"
    status="$(python3 -c 'import json,sys;print(json.load(sys.stdin)["status"])' <<<"$body")"
    case "$status" in
      succeeded|failed|stopped) echo "$body"; return 0 ;;
    esac
    sleep 0.1
  done
  echo "timeout waiting for job $1" >&2; return 1
}

echo "# mixed sample (expect 7 total, 4 kept, 3 discarded)"
JOB="$(submit "$ROOT/samples/mixed_R1.fastq.gz" "$ROOT/samples/mixed_R2.fastq.gz" 20)"
ID="$(python3 -c 'import json,sys;print(json.load(sys.stdin)["job_id"])' <<<"$JOB")"
wait_job "$ID" | python3 -m json.tool
curl -sf -o /tmp/smoke_r1.gz "$BASE/jobs/$ID/download/r1"
curl -sf -o /tmp/smoke_r2.gz "$BASE/jobs/$ID/download/r2"
python3 -c "import gzip;print('# downloaded R1 has', gzip.open('/tmp/smoke_r1.gz').read().count(b'\n')//4, 'records')"

echo "# truncated sample (expect failed gzip_corrupt, no download)"
JOB="$(submit "$ROOT/samples/truncated_R1.fastq.gz" "$ROOT/samples/truncated_R2.fastq.gz" 20)"
ID="$(python3 -c 'import json,sys;print(json.load(sys.stdin)["job_id"])' <<<"$JOB")"
wait_job "$ID" | python3 -c "import json,sys;d=json.load(sys.stdin);print(d['status'],d['error'])"
echo "# download status: $(curl -s -o /dev/null -w '%{http_code}' "$BASE/jobs/$ID/download/r1")"
