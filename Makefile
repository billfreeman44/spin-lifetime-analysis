all: analysis

analysis:
	@python analysis >> stdout.log 2>> stderr.log

clean:
	@rm -rf build
	@rm stdout.log
	@rm stderr.log

serve:
	@cd build
	@python -m http.server 8000

.PHONY: build analysis
