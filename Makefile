.PHONY: build verify test clean fast

build:
	python scripts/build.py

fast:
	python scripts/build.py --fast

verify:
	python scripts/verify.py

test:
	python -m pytest tests/ -v

clean:
	rm -rf results/*.json results/*.png figures/*.png .pytest_cache __pycache__
