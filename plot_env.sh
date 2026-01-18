#!/bin/bash

ENV_NAME="plot_env"

python3 -m venv $ENV_NAME

source $ENV_NAME/bin/activate

pip install numpy pandas matplotlib seaborn
