# AGENTS.md

## Project Overview

DTMF2NUM extracts DTMF and MF (Multi-Frequency) tones from WAV audio files and prints the detected digits. Author: Luigi Auriemma. Detects both DTMF (touch-tones: 0-9, *, #, A-D) and Bell MF (telephony signaling: 0-9, *, #, A-C, KP/ST).

## Build

### Local build
```bash
gcc -O2 -s -o dtmf2num -lm dtmf2num.c
```

`dtmf2num.c` is the **only** compilation unit. `dsp.c` and `mywav.h` are `#include`d directly — there are no separate object files or link steps beyond `-lm`.

### Cross-compilation (CI)
The CI workflow (`.github/workflows/release.yml`) cross-compiles for:
- **Linux**: x86_64, i686, aarch64, arm (gnueabihf)
- **Windows**: x86_64, i686, aarch64 (mingw-w64)
- **macOS**: arm64, x86_64

All use the same single-command pattern:
```bash
${cross}gcc -O2 -s -o dtmf2num-${suffix} -lm dtmf2num.c
```

Releases are only triggered on tags matching `v*`.

## Architecture

```
dtmf2num.c  ──#include──▶  mywav.h     (WAV reading/writing, self-contained header+impl)
             ──#include──▶  dsp.c       (Goertzel-based DTMF/MF detection, modified from Asterisk)
```

**Data flow:**
1. Parse WAV header via `mywav_data()` → get format chunk (`mywav_fmtchunk`) and data size
2. Read raw samples with `do_samples()` → allocates `int16_t[]` array
3. Convert to mono via `do_mono()` (averages channels in-place)
4. Optionally `do_dcbias()` + `do_normalize()` to clean audio (skip with `-o`)
5. Downsample to 8000 Hz via `do_downsampling()` (linear interpolation, in-place)
6. Run MF detection → `mf_detect()` using Goertzel algorithm
7. Run DTMF detection → `dtmf_detect()` using Goertzel algorithm
8. Output: detected digits printed to stdout

## Code Organization

| File | Role |
|------|------|
| `dtmf2num.c` | `main()`, audio preprocessing (mono conversion, DC bias, normalize, resampling), CLI argument parsing |
| `dsp.c` | Goertzel algorithm engine, DTMF/MF detection state machines, digit output logic. Modified from Asterisk `dsp.c` |
| `mywav.h` | Self-contained WAV file I/O — reading/writing RIFF chunks, format chunks, little-endian binary I/O |
| `Makefile` | Simple build + install targets |
| `.github/workflows/release.yml` | Multi-platform CI/CD release pipeline |

## Key Conventions & Gotchas

### Include-as-source pattern
**`dsp.c` is `#include`d, not compiled separately.** The entire tone detection engine lives in `dsp.c` and is pulled into `dtmf2num.c` via `#include "dsp.c"`. Similarly, `mywav.h` contains full function implementations, not just declarations. Any change to `dsp.c` or `mywav.h` only requires recompiling `dtmf2num.c`.

### Global `SAMPLE_RATE` variable
`SAMPLE_RATE` is declared in `dsp.c` (line ~109) and is **set from `main()` before calling any detection functions**:
```c
SAMPLE_RATE = fmt.dwSamplesPerSec;  // in dtmf2num.c:200
```
The Goertzel filter coefficients depend on this value. If you change the detection call order or add new DSP functions, ensure `SAMPLE_RATE` is set first.

### Modified Asterisk DSP
`dsp.c` is derived from Asterisk's DSP code but has been modified ("aluigi work-around" comments). Key differences from upstream:
- **Higher detection threshold**: `DTMF_THRESHOLD` and `BELL_MF_THRESHOLD` set to `800000000.0` (vs Asterisk's ~`8.0e7`)
- **Relaxed twist tests**: Reverse twist and total-energy tests are commented out
- **`DTMF_OPTIMIZED_VALUE`** set to 102 (Asterisk used 102 too, but this is explicitly tunable)
- `SAMPLE_RATE` is no longer hardcoded to 8000 — it's set dynamically

### Static parameters in dsp.c
Detection sensitivity parameters are file-scope statics in `dsp.c` (lines 79-109). If detection is too strict/loose, these are the values to adjust:
- `DTMF_THRESHOLD`, `BELL_MF_THRESHOLD` — minimum energy thresholds
- `DTMF_NORMAL_TWIST`, `BELL_MF_TWIST` — twist (level difference between tones)
- `DTMF_RELATIVE_PEAK_ROW/COL` — relative peak rejection
- `DTMF_OPTIMIZED_VALUE` — samples per detection block (affects frequency resolution)

### Platform: `stricmp` emulation
```c
#ifdef WIN32
#else
    #define stricmp strcasecmp
#endif
```
`stricmp` is used in `dtmf2num.c`. On non-Windows it's mapped to `strcasecmp`. No `<strings.h>` include is added — this works because `strcasecmp` is typically declared in `<string.h>` on glibc systems.

### Resampling is in-place
`do_downsampling()` calls `resampleData()` with the same source and destination buffer. This only works when downsampling (new rate < old rate). If `sampleRate < newSampleRate`, `resampleData()` returns early with `srcSize` to avoid corruption.

### WAV format assumptions
- Only PCM (`wFormatTag == 1`) is supported
- Supports 8, 16, 24, and 32-bit samples
- 8-bit samples are unsigned (0-255), converted to signed by `(tmp8 << 8) - 32768`
- 24-bit samples are read as little-endian 24-bit and right-shifted by 8 to fit int16

### CLI error handling
- Options are single-char flags preceded by `-` or `/` (Windows-style `/` accepted)
- `-r` consumes 3 subsequent arguments (frequency, channels, bits) — no bounds checking on `argv` index
- Exit code always 1 on error

## Testing

No test suite exists. Manual validation:
```bash
./dtmf2num test.wav                # Process WAV file
./dtmf2num -o test.wav             # Skip DC bias / normalize optimizations
./dtmf2num -r 8000 1 16 raw.pcm   # Process raw headerless PCM
./dtmf2num -w debug.wav test.wav  # Dump processed audio for debugging
```

## Usage

```
dtmf2num [options] <file.WAV>

Options:
  -r F C B   Raw headerless PCM: F=Frequency, C=Channels, B=Bits
  -o         Disable DC bias adjust and normalize
  -w FILE    Dump processed mono 8000Hz samples to FILE (debug WAV)
```

Input can be `-` for stdin. Outputs detected MF and DTMF digit sequences to stdout.
