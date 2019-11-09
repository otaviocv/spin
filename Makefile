TEST_PATH=spin

clean: 
	find . -name '__pycache__' -exec rm -r {} +
	find . -name '.ipynb_checkpoints' -exec rm -r {} +
	find . -name '.pytest_cache' -exec rm -r {} +
	rm -rf build/ dist/ spin_clustering.egg-info

lint:
	pydocstyle
	pycodestyle

init:
	pip install pipenv
	pipenv install --dev

build:
	python setup.py sdist bdist_wheel

test:
	pytest --cov-report xml:coverage.xml --cov=${TEST_PATH} --codestyle --docstyle ${TEST_PATH}
