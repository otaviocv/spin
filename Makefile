TEST_PATH=.

clean: 
	find . -name '__pycache__' -exec rm -r {} +
	find . -name '.ipynb_checkpoints' -exec rm -r {} +
	find . -name '.pytest_cache' -exec rm -r {} +

lint:
	pydocstyle
	pycodestyle

test:
	pytest --cov --codestyle --docstyle ${TEST_PATH}
