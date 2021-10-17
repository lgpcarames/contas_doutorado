################################################################################
# Automatically-generated file. Do not edit!
################################################################################

-include ../makefile.init

RM := rm -rf

# All of the sources participating in the build are defined here
-include sources.mk
-include venv/Lib/site-packages/pandas/_libs/tslibs/src/datetime/subdir.mk
-include venv/Lib/site-packages/pandas/_libs/src/ujson/python/subdir.mk
-include venv/Lib/site-packages/pandas/_libs/src/ujson/lib/subdir.mk
-include venv/Lib/site-packages/pandas/_libs/src/parser/subdir.mk
-include venv/Lib/site-packages/numpy/f2py/tests/src/array_from_pyobj/subdir.mk
-include venv/Lib/site-packages/numpy/f2py/src/subdir.mk
-include venv/Lib/site-packages/numpy/distutils/mingw/subdir.mk
-include venv/Lib/site-packages/numpy/distutils/checks/subdir.mk
-include Rascunhos/subdir.mk
-include .metadata/.plugins/org.eclipse.cdt.make.core/subdir.mk
-include subdir.mk
-include objects.mk

ifneq ($(MAKECMDGOALS),clean)
ifneq ($(strip $(C_DEPS)),)
-include $(C_DEPS)
endif
endif

-include ../makefile.defs

# Add inputs and outputs from these tool invocations to the build variables 

# All Target
all: ContasDoutorado

# Tool invocations
ContasDoutorado: $(OBJS) $(USER_OBJS)
	@echo 'Building target: $@'
	@echo 'Invoking: GCC C Linker'
	gcc -Lm -o "ContasDoutorado" $(OBJS) $(USER_OBJS) $(LIBS)
	@echo 'Finished building target: $@'
	@echo ' '

# Other Targets
clean:
	-$(RM) $(EXECUTABLES)$(OBJS)$(C_DEPS) ContasDoutorado
	-@echo ' '

.PHONY: all clean dependents

-include ../makefile.targets
