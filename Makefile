all: analysis

analysis:
	@python analysis >> stdout.log 2>> stderr.log

clean:
	@rm -rf build
	@rm stdout.log
	@rm stderr.log

serve:
	@python analysis/server.py

.PHONY: build analysis
