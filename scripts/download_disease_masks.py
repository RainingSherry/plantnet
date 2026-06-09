#!/usr/bin/env python3
"""Download All Disease Masks from Google Drive using direct export URLs - parallel version."""

import os
import sys
import time
import urllib.request
import urllib.error
import re
import threading
from concurrent.futures import ThreadPoolExecutor, as_completed

LOG_FILE = "/home/luolie/.cursor/projects/home-luolie-biopipeline-dimension-reduction-plantnet/agent-tools/91fbc849-ce6a-4447-af65-1b01a283d6f4.txt"
OUT_DIR = "/home/luolie/biopipeline/dimension-reduction/plantnet/data/AllDiseaseMasks/all_disease_masks_WACV_submission"

with open(LOG_FILE) as f:
    content = f.read()

entries = re.findall(r'Processing file ([A-Za-z0-9_-]+) (.+)', content)
seen = {}
for fid, fname in entries:
    if fid not in seen:
        seen[fid] = fname.strip()

ALL_FILES = [(fid, fname) for fid, fname in seen.items()]
print(f"Total unique files: {len(ALL_FILES)}")

os.makedirs(OUT_DIR, exist_ok=True)

# Count existing files
existing = [f for f in os.listdir(OUT_DIR) if os.path.getsize(os.path.join(OUT_DIR, f)) > 0]
print(f"Already downloaded: {len(existing)}")

to_download = [(fid, fname) for fid, fname in ALL_FILES
               if not os.path.exists(os.path.join(OUT_DIR, fname)) or os.path.getsize(os.path.join(OUT_DIR, fname)) == 0]
print(f"Need to download: {len(to_download)}")

# Thread-local storage for headers
import random

headers = {
    'User-Agent': 'Mozilla/5.0 (X11; Linux x86_64) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/120.0.0.0 Safari/537.36'
}

success = [0]
failed_list = []
lock = threading.Lock()

def download_one(args):
    idx, total, file_id, filename = args
    out_path = os.path.join(OUT_DIR, filename)
    url = f"https://drive.google.com/uc?export=download&id={file_id}"

    for attempt in range(2):
        try:
            req = urllib.request.Request(url, headers=headers)
            resp = urllib.request.urlopen(req, timeout=30)
            content = resp.read()
            resp.close()

            # Check if HTML (confirm page)
            if b'<html' in content[:200].lower():
                confirm_url = f"https://drive.usercontent.google.com/download?id={file_id}&export=download&confirm=t"
                req2 = urllib.request.Request(confirm_url, headers=headers)
                resp2 = urllib.request.urlopen(req2, timeout=30)
                content = resp2.read()
                resp2.close()

            if len(content) < 100:
                if attempt == 0:
                    time.sleep(1)
                    continue
                return ('fail', idx, total, file_id, filename, 'too small')

            with open(out_path, 'wb') as f:
                f.write(content)

            return ('ok', idx, total, file_id, filename, len(content))

        except Exception as e:
            if attempt == 0:
                time.sleep(2)
                continue
            return ('fail', idx, total, file_id, filename, str(e)[:50])

    return ('fail', idx, total, file_id, filename, 'max retries')

# Progress tracking
done_count = [0]
total_to_dl = len(to_download)
start_time = time.time()
last_print = [0]

def print_progress():
    now = time.time()
    elapsed = now - start_time[0]
    done = done_count[0]
    rate = done / elapsed if elapsed > 0 else 0
    remaining = (total_to_dl - done) / rate if rate > 0 else 0
    bar_len = 30
    filled = int(bar_len * done / total_to_dl) if total_to_dl > 0 else 0
    bar = '█' * filled + '░' * (bar_len - filled)
    print(f"\r[{bar}] {done}/{total_to_dl} ({done*100//total_to_dl}%) | {elapsed:.0f}s elapsed | ~{remaining:.0f}s remaining | {rate:.1f}/s", end='', flush=True)

print(f"Starting parallel download with 10 workers...")
print(f"Estimated time: ~{(total_to_dl * 0.5) / 60:.0f} minutes at ~2 files/sec")

with ThreadPoolExecutor(max_workers=10) as executor:
    futures = []
    for i, (fid, fname) in enumerate(to_download):
        futures.append(executor.submit(download_one, (i, total_to_dl, fid, fname)))

    for future in as_completed(futures):
        result = future.result()
        status = result[0]
        idx, total, file_id, filename, extra = result[1:]

        done_count[0] += 1

        if status == 'ok':
            success[0] += 1
        else:
            with lock:
                failed_list.append((file_id, filename, extra))
            print(f"\n  FAIL [{idx+1}/{total}]: {filename} - {extra}")

        if done_count[0] - last_print[0] >= 10 or done_count[0] == total_to_dl:
            print_progress()
            last_print[0] = done_count[0]

print(f"\n\n=== Summary ===")
print(f"Total files: {len(ALL_FILES)}")
print(f"Already existed: {len(existing)}")
print(f"Success: {success[0]}")
print(f"Failed: {len(failed_list)}")

# List failed files
if failed_list:
    print(f"\nFailed files:")
    for fid, fname, err in failed_list:
        print(f"  {fid} {fname}: {err}")

    # Save failed list to retry later
    retry_path = os.path.join(OUT_DIR, "failed_retry.txt")
    with open(retry_path, 'w') as f:
        for fid, fname, err in failed_list:
            f.write(f"{fid} {fname}\n")
    print(f"\nFailed list saved to: {retry_path}")
