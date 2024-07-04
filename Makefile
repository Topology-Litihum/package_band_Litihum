.PHONY: install clean

install:
	pip install -e .

clean:
	find . -name "*.pyc" -exec rm -f {} \;
	find . -name "__pycache__" -exec rm -rf {} \;

uninstall:
	rm ~/.local/bin/band_lithium
