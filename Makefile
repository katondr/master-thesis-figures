SHELL := /bin/bash
all: 
	source ../venv/bin/activate ; \
	python3 module/fig-1.py ; \
	python3 module/fig-2.py
clean:
	rm -f figure/*
