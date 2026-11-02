#!/usr/bin/env bash
# =============================================================================
# 00_fetch_lta.sh -- monthly passenger volumes from LTA DataMall
#
# PV/ODTrain gives trips between pairs of stations by day type and hour.
# PV/Train gives tap-ins and tap-outs by station, and is used only to check
# the totals built from the OD file. DataMall keeps the three latest months,
# and each request returns a download link that expires after five minutes.
#
# The key is read from the environment so that it never lands in the repository:
#
#   LTA_KEY=your-account-key bash scripts/00_fetch_lta.sh 202607 202608 202609
#
# Run from take-home-ex/take-home-ex3/.
# =============================================================================
set -u

if [ -z "${LTA_KEY:-}" ]; then
  echo "Set LTA_KEY to your DataMall AccountKey first." >&2
  exit 1
fi
if [ "$#" -eq 0 ]; then
  echo "Give one or more months, for example: 202607 202608 202609" >&2
  exit 1
fi

base="https://datamall2.mytransport.sg/ltaodataservice"
out="data/raw/lta"
mkdir -p "$out"

for month in "$@"; do
  for endpoint in PV/ODTrain PV/Train; do
    name="$(echo "$endpoint" | tr -d '/')_${month}.zip"
    response="$(curl -s -H "AccountKey: ${LTA_KEY}" -H "accept: application/json" \
      "${base}/${endpoint}?Date=${month}")"
    link="$(printf '%s' "$response" | grep -o 'https://[^"]*' | head -1 | sed 's/\\u0026/\&/g')"
    if [ -z "$link" ]; then
      echo "FAILED ${endpoint} ${month}: ${response}" >&2
      continue
    fi
    if curl -s -f -o "${out}/${name}" "$link"; then
      echo "saved ${out}/${name} ($(wc -c < "${out}/${name}") bytes)"
    else
      echo "FAILED download of ${endpoint} ${month}" >&2
    fi
  done
done
