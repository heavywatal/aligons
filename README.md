# aligons: pipeline for genome alignment and conservation scoring


## Installation

```
git clone https://github.com/heavywatal/aligons.git
pip3 install -v -e ./aligons
```


## Requirements

- Unix-like environment (macOS, Linux, WSL, etc.)
- Command line tools: `lastz`, `multiz`, `phastCons`, etc.
  The easiest way to install them is to use [Homebrew](https://brew.sh/):

```
brew install lastz
brew install samtools seqkit bedtools
brew install brewsci/bio/kent-tools brewsci/bio/blat
brew install heavywatal/last-bin
brew install heavywatal/multiz
brew install heavywatal/phast
```
