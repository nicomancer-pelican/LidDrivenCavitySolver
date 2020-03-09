.PHONY: clean All

All:
	@echo "----------Building project:[ LidDrivenCavitySolver - Debug ]----------"
	@cd "LidDrivenCavitySolver" && "$(MAKE)" -f  "LidDrivenCavitySolver.mk"
clean:
	@echo "----------Cleaning project:[ LidDrivenCavitySolver - Debug ]----------"
	@cd "LidDrivenCavitySolver" && "$(MAKE)" -f  "LidDrivenCavitySolver.mk" clean
