all:
	find ./src -name "*.py" -exec python "{}" ";"


clean:
	rm *.pdf
