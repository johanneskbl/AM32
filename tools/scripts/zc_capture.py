#!/usr/bin/env python3
"""Read and evaluate the zero-crossing snapshot of a firmware built with USE_ZC_CAPTURE (see Inc/zc_capture.h).

  zc_capture.py read [--elf obj/debug.elf] [--out zc.bin] [--wait 3]
      Arms the snapshot over SWD with openocd, waits, and dumps it. openocd only reads and writes RAM; the core keeps
      running. Do not halt it (gdb attach, reset halt) while the motor is powered: the phase outputs would stay
      switched on one commutation step.
  zc_capture.py analyze zc.bin [--pole-pairs 7] [--window-ms 42]
      Checks the record sequence and, for every continuous run, reports the per-phase timing offsets, the jitter,
      late or early crossings, the dependence on the PWM phase, and the rpm rate of a fit to the crossing times next
      to the former SG41 estimate (e_com_time) on the same crossings.

Only numpy is needed. Hold a constant speed while recording to read the rate statistics as noise.
"""
import argparse
import os
import shutil
import struct
import subprocess
import sys

import numpy as np

HEADER = struct.Struct("<BBBBHHHH")  # arm, state, max_pending, version, count, lost, length, pwm_arr
STEP_MASK, FORWARD, POLLING, STARTUP, GAP = 0x07, 0x08, 0x10, 0x20, 0x80
RING_SIZE = 64            # ZC_CAPTURE_RING_SIZE in Src/zc_capture.c
MAX_INTERVAL_US = 22500   # the stall restart in main.c fires at INTERVAL_TIMER_COUNT > 45000 (0.5 us ticks)
CPU_MHZ = 64              # PWM phase unit: 32 timer counts
SG_DT = 21 / 20000        # rate sample period of rpm_rate.c (1.05 ms)
SG_HALF = 20              # RPM_RATE_SG_HALF_WINDOW
REPO = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", ".."))


def find_tool(name):
    local = os.path.join(REPO, "tools", "linux", "xpack-arm-none-eabi-gcc-10.3.1-2.3", "bin", name)
    return local if os.path.exists(local) else shutil.which(name)


def read(args):
    nm = find_tool("arm-none-eabi-nm")
    if nm is None:
        sys.exit("arm-none-eabi-nm not found")
    symbols = subprocess.run([nm, "-S", args.elf], check=True, capture_output=True, text=True).stdout.split()
    try:
        i = symbols.index("zc_snapshot")
    except ValueError:
        sys.exit(f"zc_snapshot not in {args.elf}: build with USE_ZC_CAPTURE")
    addr, size = int(symbols[i - 3], 16), int(symbols[i - 2], 16)
    cfg = os.path.join(os.path.dirname(args.elf), "openocd.cfg")
    bundled = os.path.join(REPO, "tools", "linux", "openocd", "bin", "openocd")
    openocd = bundled if os.path.exists(bundled) else shutil.which("openocd")
    if openocd is None:
        sys.exit("openocd not found (expected tools/linux/openocd/bin/openocd, see 'make tools')")
    cmd = [openocd, "-f", cfg, "-c", f"mwb {addr:#x} 1", "-c", f"sleep {int(args.wait * 1000)}",
           "-c", f"dump_image {args.out} {addr:#x} {size:#x}", "-c", "shutdown"]
    print(" ".join(cmd))
    if args.dry_run:
        return
    subprocess.run(cmd, check=True)
    head = parse(args.out)[0]
    print(f"state {head['state']} ({'complete' if head['state'] == 2 else 'incomplete, use a longer --wait'}), "
          f"{head['count']} of {head['length']} records")


def parse(path):
    raw = open(path, "rb").read()
    arm, state, max_pending, version, count, lost, length, pwm_arr = HEADER.unpack_from(raw)
    if version != 1:
        sys.exit(f"unknown snapshot version {version} (not armed yet?)")
    count = min(count, length, (len(raw) - HEADER.size) // 4)
    rec = np.frombuffer(raw, "<u4", count, HEADER.size)
    head = dict(state=state, max_pending=max_pending, count=count, lost=lost, length=length, pwm_arr=pwm_arr)
    return head, (rec >> 16).astype(np.int64), ((rec >> 8) & 0xFF).astype(np.int64), (rec & 0xFF).astype(np.int64)


def segments(dt, flags):
    """continuous runs of settled crossings in one mode and direction, with consecutive commutation steps"""
    step = flags & STEP_MASK
    fwd = (flags & FORWARD) != 0
    dstep = np.where(fwd[1:], 1, -1)
    breaks = ((flags[1:] & (GAP | STARTUP)) != 0) | ((flags[1:] ^ flags[:-1]) & (POLLING | FORWARD | STARTUP) != 0) \
        | ((step[1:] - step[:-1] - dstep) % 6 != 0) | (dt <= 0) | (dt > MAX_INTERVAL_US)
    bounds = np.concatenate([[0], np.nonzero(breaks)[0] + 1, [len(flags)]])
    return [(a, b) for a, b in zip(bounds[:-1], bounds[1:]) if not flags[a] & STARTUP]


def sg_smooth_residual(t, m):
    """t minus a centred quadratic least-squares fit over 2m + 1 crossings (edges dropped)"""
    j = np.arange(-m, m + 1, dtype=float)
    w = np.linalg.pinv(np.vstack([np.ones_like(j), j, j * j]).T)[0]
    return t[m:-m] - np.convolve(t, w[::-1], "valid"), slice(m, len(t) - m)


def zc_fit_rate(t_us, keep, ts, window_us, pole_pairs):
    """rate [rpm/s] at the window centre of a quadratic fit t(k) over the crossings in (ts - W, ts]"""
    out = np.full(len(ts), np.nan)
    k_all = np.arange(len(t_us), dtype=float)
    i1s = np.searchsorted(t_us, ts, "right")
    i0s = np.searchsorted(t_us, ts - window_us, "left")
    for n, (i0, i1) in enumerate(zip(i0s, i1s)):
        sel = keep[i0:i1]
        if sel.sum() < 8:
            continue
        u = k_all[i0:i1][sel]
        u -= u.mean()
        t = t_us[i0:i1][sel] * 1e-6
        t -= t.mean()
        a = np.array([[len(u), 0, (u * u).sum()], [0, (u * u).sum(), (u ** 3).sum()],
                      [(u * u).sum(), (u ** 3).sum(), (u ** 4).sum()]])
        _, b, c = np.linalg.solve(a, [t.sum(), (u * t).sum(), (u * u * t).sum()])
        out[n] = -2 * c / b ** 3 * 10 / pole_pairs  # rpm = 60 / (6 * pole_pairs * b)
    return out


def sg41_rate(dt_us, step, t_us, ts, pole_pairs):
    """the former firmware path on the same crossings: commutation_interval filter, e_com_time, SG41 derivative"""
    ticks = np.round(2 * dt_us).astype(np.int64)
    c = int(ticks[0])
    last = c
    arr = [c] * 6
    ecom = np.empty(len(ticks), np.int64)
    for i, x in enumerate(ticks):
        arr[int(step[i]) - 1] = c
        c = (c + ((last + int(x)) >> 1)) >> 1
        last = int(x)
        ecom[i] = (sum(arr) + 4) >> 1
    idx = np.searchsorted(t_us, ts, "right") - 1
    z = 60e6 / ecom[np.clip(idx, 0, None)] / pole_pairs
    j = np.arange(-SG_HALF, SG_HALF + 1)
    h = j / (SG_DT * (j * j).sum())
    rate = np.full(len(ts), np.nan)
    rate[2 * SG_HALF:] = np.convolve(z, h[::-1], "valid")  # centred: describes ts - SG_HALF * SG_DT
    rate[: 2 * SG_HALF + 6] = np.nan  # filter still filling with crossings from before the run
    return rate


def analyze(args):
    head, t16, phase, flags = parse(args.file)
    pp = args.pole_pairs
    print(f"{args.file}: {head['count']} records (state {head['state']}), lost {head['lost']}, "
          f"ring high-water mark {head['max_pending']} of {RING_SIZE}, PWM period {head['pwm_arr'] + 1} counts "
          f"({CPU_MHZ * 1e3 / (head['pwm_arr'] + 1):.1f} kHz)")
    if head["max_pending"] > RING_SIZE * 3 // 4 or head["lost"]:
        print("  WARNING: the main loop came close to (or lost) records: increase ZC_CAPTURE_RING_SIZE before relying on it")
    dt = (np.diff(t16) & 0xFFFF).astype(float)
    print(f"  flags: {np.count_nonzero(flags & STARTUP)} start-up, {np.count_nonzero(flags & POLLING)} polling, "
          f"{np.count_nonzero(flags & GAP)} gaps, {np.count_nonzero(flags & FORWARD)} forward")
    period_pwm = (head["pwm_arr"] + 1) / 32.0
    for a, b in segments(dt, flags):
        n = b - a
        if n < 200:
            print(f"  run {a}-{b}: {n} crossings, too short")
            continue
        t = np.concatenate([[0.0], np.cumsum(dt[a:b - 1])])
        step = flags[a:b] & STEP_MASK
        interval = t[-1] / (n - 1)
        rpm = 60e6 / (interval * 6 * pp)
        te = 6 * interval
        mode = "polling" if flags[a] & POLLING else "interrupt"
        print(f"  run {a}-{b}: {n} crossings, {t[-1] / 1e3:.1f} ms, {mode} mode, mean {rpm:.0f} rpm "
              f"(interval {interval:.1f} us)")
        # timing residuals against a smooth quadratic over +-60 crossings
        m = min(60, (n - 1) // 4)
        r, sl = sg_smooth_residual(t, m)
        st = step[sl]
        off = np.array([r[st == s].mean() if np.any(st == s) else np.nan for s in range(1, 7)])
        rr = r - off[st - 1]
        mad = 1.4826 * np.median(np.abs(rr - np.median(rr)))
        print("    per-step offset [us]: " + " ".join(f"{o:+6.2f}" for o in off) +
              "   [deg el]: " + " ".join(f"{o * 360 / te:+5.2f}" for o in off))
        print(f"    jitter after offsets: {rr.std():.2f} us RMS, {mad:.2f} us robust (median abs dev), "
              f"max {np.abs(rr).max():.1f} us")
        # stateless gate as proposed for the fit: off the midpoint of the neighbours by a quarter interval
        e = t[1:-1] - 0.5 * (t[:-2] + t[2:])
        bad = np.abs(e) > interval / 4
        print(f"    late/early crossings (> 1/4 interval off the neighbour midpoint): {bad.sum()} "
              f"({100 * bad.mean():.2f} %)")
        # detection delay against the PWM phase at the crossing
        ph = np.clip(phase[a:b][sl] / period_pwm, 0, 0.999)
        bins = np.floor(ph * 8).astype(int)
        means = np.array([rr[bins == i].mean() if np.count_nonzero(bins == i) > 10 else np.nan for i in range(8)])
        explained = np.nanvar(means[bins]) / rr.var() if rr.var() > 0 else 0
        print("    residual vs PWM phase (8 bins) [us]: " + " ".join(f"{x:+5.1f}" for x in means) +
              f"   peak-peak {np.nanmax(means) - np.nanmin(means):.1f} us, {100 * explained:.0f} % of the variance")
        # Beat: the PWM phase at the crossing advances by frac(interval / PWM period) per crossing (q crossings for
        # the q-th subharmonic). A slow advance turns any phase-dependent delay into a slow timing error that no rate
        # estimate can tell from acceleration. The bins above only resolve the dependence when the phase sweeps.
        ratio = interval / (period_pwm * 32 / CPU_MHZ)
        f_zc = 1e6 / interval
        f_beat, q_beat = min((round(abs(q * ratio - round(q * ratio)) * f_zc / q, 1), q) for q in (1, 2, 3))
        note = "  WARNING: in the rate passband, the bins above under-read the dependence" if f_beat < 60 else ""
        print(f"    PWM: {ratio:.3f} periods per crossing, slowest beat {f_beat:.0f} Hz (every {q_beat} crossings){note}")
        # rate: fit over the crossing times vs the former SG41 on the same crossings, every 1.05 ms
        keep = np.ones(n, bool)
        keep[1:-1] = ~bad
        win = args.window_ms * 1e3
        ts = np.arange(t[0] + win, t[-1], SG_DT * 1e6)
        if len(ts) < 10:
            print(f"    rate: run shorter than the {args.window_ms} ms window plus 10 samples")
            continue
        fit = zc_fit_rate(t, keep, ts, win, pp)
        sg = sg41_rate(np.diff(t, prepend=t[0] - interval), step, t, ts, pp)
        for name, x in (("ZC fit", fit), ("SG41 (former)", sg)):
            v = x[np.isfinite(x)]
            if len(v) > 5:
                print(f"    rate {name:16s}: mean {v.mean():8.0f} rpm/s, std {v.std():6.0f} rpm/s over {len(v)} samples")


def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = p.add_subparsers(dest="cmd", required=True)
    r = sub.add_parser("read")
    r.add_argument("--elf", default=os.path.join(REPO, "obj", "debug.elf"))
    r.add_argument("--out", default="zc.bin")
    r.add_argument("--wait", type=float, default=3.0, help="seconds between arming and the dump")
    r.add_argument("--dry-run", action="store_true", help="only print the openocd command")
    a = sub.add_parser("analyze")
    a.add_argument("file")
    a.add_argument("--pole-pairs", type=int, default=7)
    a.add_argument("--window-ms", type=float, default=42.0)
    args = p.parse_args()
    read(args) if args.cmd == "read" else analyze(args)


if __name__ == "__main__":
    main()
