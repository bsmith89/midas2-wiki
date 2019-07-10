# Prerequisites

assert system has python 3.7 or above


# Systemwide (containerwide) install

pip3 install 'git+git://github.com/czbiohub/iggtools' --upgrade

iggtools --version


# Lint and run locally during development

cd /path/to/iggtools

pylint iggtools

python3 -m iggtools
