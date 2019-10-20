TEST_PATH=.

clean: 
	find . -name '__pycache__' -exec rm -r {} +

lint:
	pydocstyle
	pycodestyle

test:
	pytest --cov --codestyle --docstyle {TEST_PATH}
