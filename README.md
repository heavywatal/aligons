# aligons: pipeline for genome alignment and conservation scoring


## Installation

```sh
git clone https://github.com/heavywatal/aligons.git
pip3 install -v -e ./aligons
```


## Requirements

- Unix-like environment (macOS, Linux, WSL, etc.)
- Command line tools: `lastz`, `multiz`, `phastCons`, etc.
  The easiest way to install them is to use [Homebrew](https://brew.sh/):
- [UCSC kent tools](https://github.com/ucscGenomeBrowser/kent).
  Get the pre-built binaries from <https://hgdownload.soe.ucsc.edu/admin/exe/>.
  It is very difficult to compile them from source.
  Do not install them via `brewsci/bio/kent-tools` because its "bottle" is broken and old.

```sh
brew install lastz
brew install samtools bedtools
brew install zstd
brew install heavywatal/last-bin --HEAD
brew install heavywatal/multiz --HEAD
brew install heavywatal/phast --HEAD
```
