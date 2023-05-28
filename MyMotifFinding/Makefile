run:
	python MyMotifFinding.py -f "C:\\Users\\Charles Choi\\Downloads\\GRCm38.chr17.fa" \
	-transfac "C:\\Users\\Charles Choi\\Downloads\\MA0143.1.transfac" \
	-O testFolder
	"C:\\Users\\Charles Choi\\Documents\\GitHub\\MyMotifFinding\\MyMotifFinding\\Test\\peaks.txt"

install:
	pip install -r requirements.txt

build: setup.py
	python setup.py build bdist_wheel

ifeq ($(OS), Windows_NT)
clean:
	if exist "./build" rd /s /q build
	if exist "./dist" rd /s /q dist

else
clean:
	rm -rf build
	rm -rf dist
endif
