# DTMF2NUM

**DTMF2NUM** is a command-line tool that extracts DTMF (Dual-Tone Multi-Frequency) and MF (Multi-Frequency) digits from WAV audio files. It uses a modified Goertzel algorithm derived from the Asterisk DSP engine, tuned to work with damaged or low-quality audio recordings.

**DTMF2NUM** 是一个命令行工具，用于从 WAV 音频文件中提取 DTMF（双音多频）和 MF（多频）数字信号。它基于从 Asterisk DSP 引擎修改而来的 Goertzel 算法，并针对受损或低质量音频进行了优化。

---

## Features / 功能

- Detects both **DTMF** tones (0-9, \*, #, A, B, C, D) and **Bell MF** tones (0-9, \*, #, A, B, C, KP/ST)
- Supports standard WAV files (8/16/24/32-bit PCM)
- Supports raw headerless PCM data with manual format specification
- Automatic audio preprocessing: multi-channel to mono, DC bias removal, volume normalization, resampling to 8kHz
- Dump processed audio for debugging
- Cross-platform: Linux, Windows (mingw), macOS

- 同时检测 **DTMF** 拨号音（0-9, \*, #, A, B, C, D）和 **Bell MF** 信号音（0-9, \*, #, A, B, C, KP/ST）
- 支持标准 WAV 文件（8/16/24/32 位 PCM）
- 支持原始无头 PCM 数据，可手动指定格式
- 自动音频预处理：多声道转单声道、去除直流偏置、音量归一化、重采样至 8kHz
- 可导出处理后的音频用于调试
- 跨平台：Linux、Windows（mingw）、macOS

---

## Installation / 安装

### From Source / 从源码编译

```bash
git clone https://github.com/aluigi/DTMF2NUM.git
cd DTMF2NUM
gcc -O2 -s -o dtmf2num -lm dtmf2num.c
```

Requirements: GCC or compatible C compiler, `libm` (math library).

依赖：GCC 或兼容的 C 编译器、`libm`（数学库）。

### Install / 安装到系统

```bash
sudo make install   # installs to /usr/local/bin
```

### Pre-built Binaries / 预编译二进制

Download from the [Releases](https://github.com/aluigi/DTMF2NUM/releases) page. Pre-built binaries are available for:

从 [Releases](https://github.com/aluigi/DTMF2NUM/releases) 页面下载预编译二进制，支持以下平台：

| Platform / 平台 | Architectures / 架构 |
|-----------------|---------------------|
| Linux           | amd64, 386, arm64, arm |
| Windows         | amd64, 386, arm64 |
| macOS           | amd64, arm64 |

---

## Usage / 用法

```
dtmf2num [options] <file.WAV>
```

### Options / 参数

| Option | Description / 说明 |
|--------|-------------------|
| `-r F C B` | Treat input as raw headerless PCM. F=Frequency(Hz), C=Channels, B=Bits per sample / 将输入视为原始无头 PCM 数据 |
| `-o` | Disable automatic DC bias adjustment and normalization / 禁用自动直流偏置调整和归一化 |
| `-w FILE` | Dump processed audio samples to a WAV file (debug) / 将处理后的音频导出为 WAV 文件（调试用） |

### Examples / 示例

```bash
# Basic DTMF detection from a WAV file
# 从 WAV 文件中检测 DTMF
./dtmf2num recording.wav

# Skip audio preprocessing optimizations
# 跳过音频预处理优化
./dtmf2num -o recording.wav

# Process raw PCM data (8000Hz, mono, 16-bit)
# 处理原始 PCM 数据（8000Hz，单声道，16 位）
./dtmf2num -r 8000 1 16 raw_audio.pcm

# Read from stdin
# 从标准输入读取
cat recording.wav | ./dtmf2num -

# Dump processed audio for debugging
# 导出处理后的音频用于调试
./dtmf2num -w debug.wav recording.wav
```

### Output / 输出示例

```
- open recording.wav
  wave size      160000
  format tag     1
  channels:      1
  samples/sec:   44100
  avg/bytes/sec: 88200
  block align:   2
  bits:          16
  samples:       80000
  bias adjust:   -12
  volume peaks:  -24576 28672
  normalize:     4063
  resampling to: 8000hz

- MF numbers:    1234567890*#
- DTMF numbers:  1234567890*#
```

---

## How It Works / 工作原理

### Signal Processing Pipeline / 信号处理流程

1. **WAV Parsing** — Parse the RIFF/WAV header to extract format info (sample rate, channels, bit depth) and locate audio data.

   **WAV 解析** — 解析 RIFF/WAV 文件头，提取格式信息（采样率、声道数、位深度）并定位音频数据。

2. **Sample Loading** — Read all audio samples into a signed 16-bit integer array. 8-bit unsigned samples are converted; 24/32-bit samples are downscaled.

   **采样加载** — 将所有音频样本读入有符号 16 位整数数组。8 位无符号样本会被转换；24/32 位样本会缩小。

3. **Mono Conversion** — Multi-channel audio is mixed down to mono by averaging all channels per sample.

   **单声道转换** — 多声道音频通过对每个采样点的所有声道求平均值，混合为单声道。

4. **DC Bias Removal** — Compute the midpoint between max positive and max negative amplitude, subtract it from all samples to center the waveform around zero.

   **直流偏置去除** — 计算最大正振幅和最大负振幅的中点，从所有采样中减去该值，使波形围绕零对称。

5. **Normalization** — Scale the waveform so peak amplitude reaches ±32767 (full 16-bit range), improving detection on quiet recordings.

   **归一化** — 缩放波形使峰值振幅达到 ±32767（16 位全范围），提升对低音量录音的检测效果。

6. **Downsampling** — Resample to 8000 Hz using linear interpolation. DTMF/MF standards are based on 8kHz telephony sample rate.

   **降采样** — 使用线性插值重采样至 8000 Hz。DTMF/MF 标准基于 8kHz 电话采样率。

7. **Tone Detection** — Apply the Goertzel algorithm to measure energy at each DTMF/MF frequency. A digit is detected when two specific frequencies (one row, one column) dominate the signal above thresholds.

   **音调检测** — 使用 Goertzel 算法测量每个 DTMF/MF 频率的能量。当两个特定频率（一个行频、一个列频）在信号中占主导地位且超过阈值时，判定检测到一个数字。

### DTMF Frequencies / DTMF 频率表

|         | 1209 Hz | 1336 Hz | 1477 Hz | 1633 Hz |
|---------|---------|---------|---------|---------|
| **697 Hz**  | 1       | 2       | 3       | A       |
| **770 Hz**  | 4       | 5       | 6       | B       |
| **852 Hz**  | 7       | 8       | 9       | C       |
| **941 Hz**  | \*      | 0       | #       | D       |

### Bell MF Frequencies / Bell MF 频率表

Bell MF uses 6 tones (700, 900, 1100, 1300, 1500, 1700 Hz), with each digit composed of exactly two frequencies. Bell MF 使用 6 个频率（700, 900, 1100, 1300, 1500, 1700 Hz），每个数字由恰好两个频率组成。

---

## License / 许可证

GNU General Public License v2.0. See the source files for full license text.

GNU 通用公共许可证 v2.0。完整许可文本见源文件。

### Third-party Code / 第三方代码

- **dsp.c**: Modified from Asterisk (GPL v2, originally by Digium, Mark Spencer, Steve Underwood). Detection thresholds and twist tests have been adjusted for robustness with damaged audio.

  **dsp.c**: 修改自 Asterisk（GPL v2，原作者 Digium、Mark Spencer、Steve Underwood）。检测阈值和扭曲测试已调整以增强对受损音频的鲁棒性。

- **mywav.h**: Custom lightweight WAV parsing library by Luigi Auriemma.

  **mywav.h**: Luigi Auriemma 自制的轻量级 WAV 解析库。

---

## Author / 作者

Luigi Auriemma

- Email: aluigi@autistici.org
- Website: [aluigi.org](https://aluigi.org)
- GitHub: [@aluigi](https://github.com/aluigi)

## Contributors / 贡献者

**[Moxin1044](https://github.com/Moxin1044)** — CI/CD pipeline setup, multi-platform automated builds (GitHub Actions). / CI/CD 流水线搭建、多平台自动化构建（GitHub Actions）。

---

## Contributing / 贡献

Bug reports and pull requests are welcome. When modifying `dsp.c`, note that detection parameters are tuned for damaged audio — thresholds are intentionally higher than standard Asterisk values. Test against real-world recordings before adjusting.

欢迎提交 Bug 报告和 Pull Request。修改 `dsp.c` 时请注意，检测参数已针对受损音频进行了调优——阈值故意设得比标准 Asterisk 值更高。调整前请用真实录音进行测试。
