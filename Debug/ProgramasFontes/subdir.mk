################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
C_SRCS += \
../ProgramasFontes/LevyDado.c \
../ProgramasFontes/prl-comment-concentric.c 

OBJS += \
./ProgramasFontes/LevyDado.o \
./ProgramasFontes/prl-comment-concentric.o 

C_DEPS += \
./ProgramasFontes/LevyDado.d \
./ProgramasFontes/prl-comment-concentric.d 


# Each subdirectory must supply rules for building sources it contributes
ProgramasFontes/%.o: ../ProgramasFontes/%.c
	@echo 'Building file: $<'
	@echo 'Invoking: GCC C Compiler'
	gcc -O0 -g3 -Wall -c -fmessage-length=0 -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@)" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


