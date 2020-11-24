# VRE Executor for ProGeny

Example pipelines file that is ready to run in the VRE matching the code in the HowTo documentation.

This repo structure fits well any R-based tool and might be used as the base template for future tools. It should have all of the base functionality and is set up for unit testing and with pylint to ensure code clarity.

## Requirements

- Python 3.6 or later (Recommended 3.7)
- Python3.6-pip, Python3.6-dev and Python3.6-venv or later
- R-4.0
- R library to be integrated: [PROGENy](https://github.com/saezlab/progeny)

Install python enviroment and R:

```bash
sudo apt update
sudo apt install python3.6 
sudo apt install python3.6-pip python3.6-dev python3.6-venv
sudo apt install r-base
```

Install the particular R library and its dependencies:

```bash
sudo su
su www-data -s /bin/bash
./Rprogeny.sh
su user
```

## Installation

Directly from GitHub:

```bash
cd $HOME
git clone $THIS_REPO
```

Create the Python environment:

```bash
python3 -m venv $HOME/$THIS_REPO/venv
source venv/bin/activate
pip install -r requirements.txt
```

## Run the Wrapper
```bash
./VRE_RUNNER --config tests/basic/config.json --in_metadata tests/basic/in_metadata.json --out_metadata out_metadata.json --log_file VRE_RUNNER.log
```
