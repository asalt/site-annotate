install:
	pip install .

unit_tests:
	pytest tests/unit_tests/

integration_tests:
	pytest tests/unit_tests/
