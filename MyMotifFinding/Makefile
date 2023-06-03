install:
	pip install -r requirements.txt

build: setup.py
	python ../setup.py build bdist_wheel

ifeq ($(OS), Windows_NT)
clean:
	if exist "./build" rd /s /q build
	if exist "./dist" rd /s /q dist
	
else
clean:
	rm -rf build
	rm -rf dist
endif
