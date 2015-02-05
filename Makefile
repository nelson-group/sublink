all:
	@cd Descendants; \
	make; \
	cd ../SubhaloTrees; \
	make; \

clean:
	@cd Descendants; \
	make clean; \
	cd ../SubhaloTrees; \
	make clean
