all:
	find -name "*.py" -exec python "{}" ";"


clean:
	rm *.pdf
