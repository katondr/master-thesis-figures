SHELL := /bin/bash
all: 
	source ../venv/bin/activate ; \
	python3 example/example-1.py ; \
	python3 example/example-2.py ; \
	python3 example/example-3.py ; \
	python3 example/example-4.py ; \
	python3 example/example-5.py ; \
	python3 example/example-6.py ; \
	python3 example/example-7.py ; \
	python3 example/example-8.py 
clean:
	rm -f figure/*
