TEST_PATH=spin

clean: 
	find . -name '__pycache__' -exec rm -r {} +
	find . -name '.ipynb_checkpoints' -exec rm -r {} +
	find . -name '.pytest_cache' -exec rm -r {} +
	rm -rf build/ dist/

lint:
	pydocstyle
	pycodestyle

init:
	pip install pipenv
	pipenv install --dev

build:
	python setup.py sdist bdist_wheel

test:
	pytest --cov --codestyle --docstyle ${TEST_PATH}
