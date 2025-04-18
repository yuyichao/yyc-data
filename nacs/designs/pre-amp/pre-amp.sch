EESchema Schematic File Version 2
LIBS:power
LIBS:device
LIBS:switches
LIBS:relays
LIBS:motors
LIBS:transistors
LIBS:conn
LIBS:linear
LIBS:regul
LIBS:74xx
LIBS:cmos4000
LIBS:adc-dac
LIBS:memory
LIBS:xilinx
LIBS:microcontrollers
LIBS:dsp
LIBS:microchip
LIBS:analog_switches
LIBS:motorola
LIBS:texas
LIBS:intel
LIBS:audio
LIBS:interface
LIBS:digital-audio
LIBS:philips
LIBS:display
LIBS:cypress
LIBS:siliconi
LIBS:opto
LIBS:atmel
LIBS:contrib
LIBS:valves
LIBS:pre-amp-cache
EELAYER 25 0
EELAYER END
$Descr USLetter 11000 8500
encoding utf-8
Sheet 1 1
Title ""
Date ""
Rev ""
Comp ""
Comment1 ""
Comment2 ""
Comment3 ""
Comment4 ""
$EndDescr
NoConn ~ 4500 2650
NoConn ~ 4500 2850
$Comp
L Conn_Coaxial J1
U 1 1 5AEB2F69
P 3250 2550
F 0 "J1" H 3260 2670 50  0000 C CNN
F 1 "Conn_Coaxial" V 3365 2550 50  0000 C CNN
F 2 "Connectors2:BNC_Socket_Right-Angle_LargePads" H 3250 2550 50  0001 C CNN
F 3 "" H 3250 2550 50  0001 C CNN
	1    3250 2550
	-1   0    0    -1  
$EndComp
$Comp
L GND #PWR01
U 1 1 5AEB2FF5
P 5000 3150
F 0 "#PWR01" H 5000 2900 50  0001 C CNN
F 1 "GND" H 5000 3000 50  0000 C CNN
F 2 "" H 5000 3150 50  0001 C CNN
F 3 "" H 5000 3150 50  0001 C CNN
	1    5000 3150
	1    0    0    -1  
$EndComp
$Comp
L C C1
U 1 1 5AEB3124
P 3800 5400
F 0 "C1" H 3825 5500 50  0000 L CNN
F 1 "0.1u" H 3825 5300 50  0000 L CNN
F 2 "Capacitors_SMD:C_0603" H 3838 5250 50  0001 C CNN
F 3 "" H 3800 5400 50  0001 C CNN
	1    3800 5400
	1    0    0    -1  
$EndComp
$Comp
L C C3
U 1 1 5AEB3174
P 4050 5400
F 0 "C3" H 4075 5500 50  0000 L CNN
F 1 "10u" H 4075 5300 50  0000 L CNN
F 2 "Capacitors_SMD:C_0603" H 4088 5250 50  0001 C CNN
F 3 "" H 4050 5400 50  0001 C CNN
	1    4050 5400
	1    0    0    -1  
$EndComp
$Comp
L -15V #PWR10
U 1 1 5AEB31DF
P 4800 3150
F 0 "#PWR10" H 4800 3250 50  0001 C CNN
F 1 "-15V" H 4800 3300 50  0000 C CNN
F 2 "" H 4800 3150 50  0001 C CNN
F 3 "" H 4800 3150 50  0001 C CNN
	1    4800 3150
	-1   0    0    1   
$EndComp
$Comp
L GND #PWR02
U 1 1 5AEB3300
P 3800 5650
F 0 "#PWR02" H 3800 5400 50  0001 C CNN
F 1 "GND" H 3800 5500 50  0000 C CNN
F 2 "" H 3800 5650 50  0001 C CNN
F 3 "" H 3800 5650 50  0001 C CNN
	1    3800 5650
	1    0    0    -1  
$EndComp
$Comp
L -15V #PWR2
U 1 1 5AEB33C1
P 3800 5150
F 0 "#PWR2" H 3800 5250 50  0001 C CNN
F 1 "-15V" H 3800 5300 50  0000 C CNN
F 2 "" H 3800 5150 50  0001 C CNN
F 3 "" H 3800 5150 50  0001 C CNN
	1    3800 5150
	1    0    0    -1  
$EndComp
$Comp
L +15V #PWR03
U 1 1 5AEB347E
P 4800 2350
F 0 "#PWR03" H 4800 2200 50  0001 C CNN
F 1 "+15V" H 4800 2490 50  0000 C CNN
F 2 "" H 4800 2350 50  0001 C CNN
F 3 "" H 4800 2350 50  0001 C CNN
	1    4800 2350
	1    0    0    -1  
$EndComp
$Comp
L +15V #PWR04
U 1 1 5AEB34BD
P 4300 5150
F 0 "#PWR04" H 4300 5000 50  0001 C CNN
F 1 "+15V" H 4300 5290 50  0000 C CNN
F 2 "" H 4300 5150 50  0001 C CNN
F 3 "" H 4300 5150 50  0001 C CNN
	1    4300 5150
	1    0    0    -1  
$EndComp
$Comp
L C C4
U 1 1 5AEB3510
P 4300 5400
F 0 "C4" H 4325 5500 50  0000 L CNN
F 1 "0.1u" H 4325 5300 50  0000 L CNN
F 2 "Capacitors_SMD:C_0603" H 4338 5250 50  0001 C CNN
F 3 "" H 4300 5400 50  0001 C CNN
	1    4300 5400
	1    0    0    -1  
$EndComp
$Comp
L C C5
U 1 1 5AEB3546
P 4550 5400
F 0 "C5" H 4575 5500 50  0000 L CNN
F 1 "10u" H 4575 5300 50  0000 L CNN
F 2 "Capacitors_SMD:C_0603" H 4588 5250 50  0001 C CNN
F 3 "" H 4550 5400 50  0001 C CNN
	1    4550 5400
	1    0    0    -1  
$EndComp
Wire Wire Line
	4800 3050 4800 3150
Wire Wire Line
	5000 3050 5000 3150
Wire Wire Line
	3800 5550 3800 5650
Wire Wire Line
	4050 5550 4050 5600
Wire Wire Line
	3800 5600 4550 5600
Connection ~ 3800 5600
Wire Wire Line
	4050 5250 4050 5200
Wire Wire Line
	4050 5200 3800 5200
Wire Wire Line
	3800 5150 3800 5250
Connection ~ 3800 5200
Wire Wire Line
	4800 2350 4800 2450
Wire Wire Line
	4300 5150 4300 5250
Wire Wire Line
	4300 5600 4300 5550
Connection ~ 4050 5600
Wire Wire Line
	4300 5200 4550 5200
Wire Wire Line
	4550 5200 4550 5250
Connection ~ 4300 5200
Wire Wire Line
	4550 5600 4550 5550
Connection ~ 4300 5600
$Comp
L R R1
U 1 1 5AEB39CD
P 3700 2550
F 0 "R1" V 3780 2550 50  0000 C CNN
F 1 "1.01k" V 3700 2550 50  0000 C CNN
F 2 "Resistors_SMD:R_0603" V 3630 2550 50  0001 C CNN
F 3 "" H 3700 2550 50  0001 C CNN
	1    3700 2550
	0    1    1    0   
$EndComp
$Comp
L C C2
U 1 1 5AEB3A3A
P 3950 2750
F 0 "C2" H 3975 2850 50  0000 L CNN
F 1 "2200p" H 3975 2650 50  0000 L CNN
F 2 "Capacitors_SMD:C_0603" H 3988 2600 50  0001 C CNN
F 3 "" H 3950 2750 50  0001 C CNN
	1    3950 2750
	1    0    0    -1  
$EndComp
$Comp
L R R2
U 1 1 5AEB3AAC
P 3700 2950
F 0 "R2" V 3780 2950 50  0000 C CNN
F 1 "1.01k" V 3700 2950 50  0000 C CNN
F 2 "Resistors_SMD:R_0603" V 3630 2950 50  0001 C CNN
F 3 "" H 3700 2950 50  0001 C CNN
	1    3700 2950
	0    1    1    0   
$EndComp
Wire Wire Line
	3250 2950 3550 2950
Wire Wire Line
	3850 2550 4500 2550
Wire Wire Line
	4500 2950 3850 2950
Wire Wire Line
	3950 2900 3950 2950
Connection ~ 3950 2950
Wire Wire Line
	3950 2600 3950 2550
Connection ~ 3950 2550
$Comp
L R R5
U 1 1 5AEB3FD7
P 5600 2750
F 0 "R5" V 5680 2750 50  0000 C CNN
F 1 "R" V 5600 2750 50  0000 C CNN
F 2 "Resistors_SMD:R_0603" V 5530 2750 50  0001 C CNN
F 3 "" H 5600 2750 50  0001 C CNN
	1    5600 2750
	0    1    1    0   
$EndComp
$Comp
L C C12
U 1 1 5AEB405F
P 6050 3000
F 0 "C12" H 6075 3100 50  0000 L CNN
F 1 "C" H 6075 2900 50  0000 L CNN
F 2 "Capacitors_SMD:C_0603" H 6088 2850 50  0001 C CNN
F 3 "" H 6050 3000 50  0001 C CNN
	1    6050 3000
	1    0    0    -1  
$EndComp
$Comp
L R R7
U 1 1 5AEB40A1
P 6300 2750
F 0 "R7" V 6380 2750 50  0000 C CNN
F 1 "R" V 6300 2750 50  0000 C CNN
F 2 "Resistors_SMD:R_0603" V 6230 2750 50  0001 C CNN
F 3 "" H 6300 2750 50  0001 C CNN
	1    6300 2750
	0    1    1    0   
$EndComp
$Comp
L C C13
U 1 1 5AEB41BC
P 7000 2150
F 0 "C13" H 7025 2250 50  0000 L CNN
F 1 "C" H 7025 2050 50  0000 L CNN
F 2 "Capacitors_SMD:C_0603" H 7038 2000 50  0001 C CNN
F 3 "" H 7000 2150 50  0001 C CNN
	1    7000 2150
	0    1    1    0   
$EndComp
$Comp
L R R6
U 1 1 5AEB4241
P 6300 2000
F 0 "R6" V 6380 2000 50  0000 C CNN
F 1 "R" V 6300 2000 50  0000 C CNN
F 2 "Resistors_SMD:R_0603" V 6230 2000 50  0001 C CNN
F 3 "" H 6300 2000 50  0001 C CNN
	1    6300 2000
	0    1    1    0   
$EndComp
$Comp
L GND #PWR05
U 1 1 5AEB4363
P 6050 3250
F 0 "#PWR05" H 6050 3000 50  0001 C CNN
F 1 "GND" H 6050 3100 50  0000 C CNN
F 2 "" H 6050 3250 50  0001 C CNN
F 3 "" H 6050 3250 50  0001 C CNN
	1    6050 3250
	1    0    0    -1  
$EndComp
$Comp
L GND #PWR06
U 1 1 5AEB439E
P 6500 3250
F 0 "#PWR06" H 6500 3000 50  0001 C CNN
F 1 "GND" H 6500 3100 50  0000 C CNN
F 2 "" H 6500 3250 50  0001 C CNN
F 3 "" H 6500 3250 50  0001 C CNN
	1    6500 3250
	1    0    0    -1  
$EndComp
Wire Wire Line
	5300 2750 5450 2750
Wire Wire Line
	5750 2750 6150 2750
Wire Wire Line
	6450 2750 6550 2750
Wire Wire Line
	6050 2000 6050 2850
Connection ~ 6050 2750
Wire Wire Line
	6050 2000 6150 2000
Wire Wire Line
	7150 2850 7350 2850
Wire Wire Line
	7250 2000 7250 2850
Wire Wire Line
	7250 2150 7150 2150
Wire Wire Line
	7250 2000 6450 2000
Connection ~ 7250 2150
Wire Wire Line
	6500 2750 6500 2150
Wire Wire Line
	6500 2150 6850 2150
Connection ~ 6500 2750
Wire Wire Line
	6500 3250 6500 2950
Wire Wire Line
	6500 2950 6550 2950
Wire Wire Line
	6050 3250 6050 3150
$Comp
L C C6
U 1 1 5AEB47EF
P 4900 5400
F 0 "C6" H 4925 5500 50  0000 L CNN
F 1 "0.1u" H 4925 5300 50  0000 L CNN
F 2 "Capacitors_SMD:C_0603" H 4938 5250 50  0001 C CNN
F 3 "" H 4900 5400 50  0001 C CNN
	1    4900 5400
	1    0    0    -1  
$EndComp
$Comp
L C C8
U 1 1 5AEB47F5
P 5150 5400
F 0 "C8" H 5175 5500 50  0000 L CNN
F 1 "10u" H 5175 5300 50  0000 L CNN
F 2 "Capacitors_SMD:C_0603" H 5188 5250 50  0001 C CNN
F 3 "" H 5150 5400 50  0001 C CNN
	1    5150 5400
	1    0    0    -1  
$EndComp
$Comp
L GND #PWR07
U 1 1 5AEB47FB
P 4900 5650
F 0 "#PWR07" H 4900 5400 50  0001 C CNN
F 1 "GND" H 4900 5500 50  0000 C CNN
F 2 "" H 4900 5650 50  0001 C CNN
F 3 "" H 4900 5650 50  0001 C CNN
	1    4900 5650
	1    0    0    -1  
$EndComp
$Comp
L -15V #PWR11
U 1 1 5AEB4801
P 4900 5150
F 0 "#PWR11" H 4900 5250 50  0001 C CNN
F 1 "-15V" H 4900 5300 50  0000 C CNN
F 2 "" H 4900 5150 50  0001 C CNN
F 3 "" H 4900 5150 50  0001 C CNN
	1    4900 5150
	1    0    0    -1  
$EndComp
$Comp
L +15V #PWR08
U 1 1 5AEB4807
P 5400 5150
F 0 "#PWR08" H 5400 5000 50  0001 C CNN
F 1 "+15V" H 5400 5290 50  0000 C CNN
F 2 "" H 5400 5150 50  0001 C CNN
F 3 "" H 5400 5150 50  0001 C CNN
	1    5400 5150
	1    0    0    -1  
$EndComp
$Comp
L C C9
U 1 1 5AEB480D
P 5400 5400
F 0 "C9" H 5425 5500 50  0000 L CNN
F 1 "0.1u" H 5425 5300 50  0000 L CNN
F 2 "Capacitors_SMD:C_0603" H 5438 5250 50  0001 C CNN
F 3 "" H 5400 5400 50  0001 C CNN
	1    5400 5400
	1    0    0    -1  
$EndComp
$Comp
L C C11
U 1 1 5AEB4813
P 5650 5400
F 0 "C11" H 5675 5500 50  0000 L CNN
F 1 "10u" H 5675 5300 50  0000 L CNN
F 2 "Capacitors_SMD:C_0603" H 5688 5250 50  0001 C CNN
F 3 "" H 5650 5400 50  0001 C CNN
	1    5650 5400
	1    0    0    -1  
$EndComp
Wire Wire Line
	4900 5550 4900 5650
Wire Wire Line
	5150 5550 5150 5600
Wire Wire Line
	4900 5600 5650 5600
Connection ~ 4900 5600
Wire Wire Line
	5150 5250 5150 5200
Wire Wire Line
	5150 5200 4900 5200
Wire Wire Line
	4900 5150 4900 5250
Connection ~ 4900 5200
Wire Wire Line
	5400 5150 5400 5250
Wire Wire Line
	5400 5600 5400 5550
Connection ~ 5150 5600
Wire Wire Line
	5400 5200 5650 5200
Wire Wire Line
	5650 5200 5650 5250
Connection ~ 5400 5200
Wire Wire Line
	5650 5600 5650 5550
Connection ~ 5400 5600
$Comp
L +15V #PWR09
U 1 1 5AEB4851
P 6750 3250
F 0 "#PWR09" H 6750 3100 50  0001 C CNN
F 1 "+15V" H 6750 3390 50  0000 C CNN
F 2 "" H 6750 3250 50  0001 C CNN
F 3 "" H 6750 3250 50  0001 C CNN
	1    6750 3250
	-1   0    0    1   
$EndComp
$Comp
L -15V #PWR20
U 1 1 5AEB4898
P 6750 2450
F 0 "#PWR20" H 6750 2550 50  0001 C CNN
F 1 "-15V" H 6750 2600 50  0000 C CNN
F 2 "" H 6750 2450 50  0001 C CNN
F 3 "" H 6750 2450 50  0001 C CNN
	1    6750 2450
	1    0    0    -1  
$EndComp
Wire Wire Line
	6750 2450 6750 2550
Wire Wire Line
	6750 3150 6750 3250
$Comp
L R R8
U 1 1 5AEB4F3A
P 7500 2850
F 0 "R8" V 7580 2850 50  0000 C CNN
F 1 "0" V 7500 2850 50  0000 C CNN
F 2 "Resistors_SMD:R_0603" V 7430 2850 50  0001 C CNN
F 3 "" H 7500 2850 50  0001 C CNN
	1    7500 2850
	0    1    1    0   
$EndComp
$Comp
L C C14
U 1 1 5AEB4FCB
P 7750 3100
F 0 "C14" H 7775 3200 50  0000 L CNN
F 1 "0" H 7775 3000 50  0000 L CNN
F 2 "Capacitors_SMD:C_0603" H 7788 2950 50  0001 C CNN
F 3 "" H 7750 3100 50  0001 C CNN
	1    7750 3100
	1    0    0    -1  
$EndComp
$Comp
L GND #PWR010
U 1 1 5AEB505B
P 7750 3350
F 0 "#PWR010" H 7750 3100 50  0001 C CNN
F 1 "GND" H 7750 3200 50  0000 C CNN
F 2 "" H 7750 3350 50  0001 C CNN
F 3 "" H 7750 3350 50  0001 C CNN
	1    7750 3350
	1    0    0    -1  
$EndComp
Connection ~ 7250 2850
Wire Wire Line
	7650 2850 7950 2850
Wire Wire Line
	7750 2850 7750 2950
Wire Wire Line
	7750 3250 7750 3350
$Comp
L Conn_Coaxial J3
U 1 1 5AEB5216
P 8100 2850
F 0 "J3" H 8110 2970 50  0000 C CNN
F 1 "Conn_Coaxial" V 8215 2850 50  0000 C CNN
F 2 "Connectors2:BNC_Socket_Right-Angle_LargePads" H 8100 2850 50  0001 C CNN
F 3 "" H 8100 2850 50  0001 C CNN
	1    8100 2850
	1    0    0    -1  
$EndComp
Connection ~ 7750 2850
Wire Wire Line
	8100 3050 8100 3300
Wire Wire Line
	8100 3300 7750 3300
Connection ~ 7750 3300
$Comp
L AD8429 U1
U 1 1 5AEB57AE
P 4900 2750
F 0 "U1" H 5050 3050 50  0000 C CNN
F 1 "AD8429" H 5050 2950 50  0000 C CNN
F 2 "Housings_SOIC:SOIC-8_3.9x4.9mm_Pitch1.27mm" H 4600 2750 50  0001 C CNN
F 3 "" H 5250 2350 50  0001 C CNN
	1    4900 2750
	1    0    0    -1  
$EndComp
$Comp
L Conn_01x03 J2
U 1 1 5AEB6541
P 4200 6250
F 0 "J2" H 4200 6450 50  0000 C CNN
F 1 "Conn_01x03" H 4200 6050 50  0000 C CNN
F 2 "cmr_link:3-pins_150mil" H 4200 6250 50  0001 C CNN
F 3 "" H 4200 6250 50  0001 C CNN
	1    4200 6250
	1    0    0    -1  
$EndComp
$Comp
L R R3
U 1 1 5AEB681A
P 3750 6250
F 0 "R3" V 3830 6250 50  0000 C CNN
F 1 "0" V 3750 6250 50  0000 C CNN
F 2 "Resistors_SMD:R_0603" V 3680 6250 50  0001 C CNN
F 3 "" H 3750 6250 50  0001 C CNN
	1    3750 6250
	0    1    1    0   
$EndComp
$Comp
L GND #PWR011
U 1 1 5AEB68A7
P 3550 6300
F 0 "#PWR011" H 3550 6050 50  0001 C CNN
F 1 "GND" H 3550 6150 50  0000 C CNN
F 2 "" H 3550 6300 50  0001 C CNN
F 3 "" H 3550 6300 50  0001 C CNN
	1    3550 6300
	1    0    0    -1  
$EndComp
$Comp
L -15V #PWR5
U 1 1 5AEB68FD
P 3900 6500
F 0 "#PWR5" H 3900 6600 50  0001 C CNN
F 1 "-15V" H 3900 6650 50  0000 C CNN
F 2 "" H 3900 6500 50  0001 C CNN
F 3 "" H 3900 6500 50  0001 C CNN
	1    3900 6500
	-1   0    0    1   
$EndComp
$Comp
L +15V #PWR012
U 1 1 5AEB69C1
P 3900 6050
F 0 "#PWR012" H 3900 5900 50  0001 C CNN
F 1 "+15V" H 3900 6190 50  0000 C CNN
F 2 "" H 3900 6050 50  0001 C CNN
F 3 "" H 3900 6050 50  0001 C CNN
	1    3900 6050
	1    0    0    -1  
$EndComp
Wire Wire Line
	3900 6050 3900 6150
Wire Wire Line
	3900 6150 4000 6150
Wire Wire Line
	3900 6250 4000 6250
Wire Wire Line
	4000 6350 3900 6350
Wire Wire Line
	3900 6350 3900 6500
Wire Wire Line
	3600 6250 3550 6250
Wire Wire Line
	3550 6250 3550 6300
$Comp
L GND #PWR013
U 1 1 5AEB6FED
P 4900 6650
F 0 "#PWR013" H 4900 6400 50  0001 C CNN
F 1 "GND" H 4900 6500 50  0000 C CNN
F 2 "" H 4900 6650 50  0001 C CNN
F 3 "" H 4900 6650 50  0001 C CNN
	1    4900 6650
	1    0    0    -1  
$EndComp
$Comp
L -15V #PWR13
U 1 1 5AEB6FF3
P 4900 6150
F 0 "#PWR13" H 4900 6250 50  0001 C CNN
F 1 "-15V" H 4900 6300 50  0000 C CNN
F 2 "" H 4900 6150 50  0001 C CNN
F 3 "" H 4900 6150 50  0001 C CNN
	1    4900 6150
	1    0    0    -1  
$EndComp
$Comp
L +15V #PWR014
U 1 1 5AEB6FF9
P 5400 6150
F 0 "#PWR014" H 5400 6000 50  0001 C CNN
F 1 "+15V" H 5400 6290 50  0000 C CNN
F 2 "" H 5400 6150 50  0001 C CNN
F 3 "" H 5400 6150 50  0001 C CNN
	1    5400 6150
	1    0    0    -1  
$EndComp
$Comp
L CP C10
U 1 1 5AEB7242
P 5400 6400
F 0 "C10" H 5425 6500 50  0000 L CNN
F 1 "CP" H 5425 6300 50  0000 L CNN
F 2 "Capacitors_THT:CP_Radial_D5.0mm_P2.00mm" H 5438 6250 50  0001 C CNN
F 3 "" H 5400 6400 50  0001 C CNN
	1    5400 6400
	1    0    0    -1  
$EndComp
$Comp
L CP C7
U 1 1 5AEB72B5
P 4900 6400
F 0 "C7" H 4925 6500 50  0000 L CNN
F 1 "CP" H 4925 6300 50  0000 L CNN
F 2 "Capacitors_THT:CP_Radial_D5.0mm_P2.00mm" H 4938 6250 50  0001 C CNN
F 3 "" H 4900 6400 50  0001 C CNN
	1    4900 6400
	-1   0    0    1   
$EndComp
Wire Wire Line
	5400 6150 5400 6250
Wire Wire Line
	5400 6550 5400 6600
Wire Wire Line
	5400 6600 4900 6600
Wire Wire Line
	4900 6550 4900 6650
Connection ~ 4900 6600
Wire Wire Line
	4900 6250 4900 6150
$Comp
L OP179GS U2
U 1 1 5AEB86F5
P 6850 2850
F 0 "U2" H 6850 3100 50  0000 L CNN
F 1 "OP211" H 6850 3000 50  0000 L CNN
F 2 "Housings_SOIC:SOIC-8_3.9x4.9mm_Pitch1.27mm" H 6850 2850 50  0001 C CNN
F 3 "" H 7000 3000 50  0001 C CNN
	1    6850 2850
	1    0    0    1   
$EndComp
$Comp
L R R4
U 1 1 5B01BCCF
P 4250 2300
F 0 "R4" V 4330 2300 50  0000 C CNN
F 1 "100k" V 4250 2300 50  0000 C CNN
F 2 "Resistors_SMD:R_0603" V 4180 2300 50  0001 C CNN
F 3 "" H 4250 2300 50  0001 C CNN
	1    4250 2300
	-1   0    0    1   
$EndComp
$Comp
L R R9
U 1 1 5B01BDDF
P 4250 3200
F 0 "R9" V 4330 3200 50  0000 C CNN
F 1 "100k" V 4250 3200 50  0000 C CNN
F 2 "Resistors_SMD:R_0603" V 4180 3200 50  0001 C CNN
F 3 "" H 4250 3200 50  0001 C CNN
	1    4250 3200
	-1   0    0    1   
$EndComp
$Comp
L GND #PWR015
U 1 1 5B01C301
P 4250 3450
F 0 "#PWR015" H 4250 3200 50  0001 C CNN
F 1 "GND" H 4250 3300 50  0000 C CNN
F 2 "" H 4250 3450 50  0001 C CNN
F 3 "" H 4250 3450 50  0001 C CNN
	1    4250 3450
	1    0    0    -1  
$EndComp
$Comp
L GND #PWR016
U 1 1 5B01C363
P 4450 2150
F 0 "#PWR016" H 4450 1900 50  0001 C CNN
F 1 "GND" H 4450 2000 50  0000 C CNN
F 2 "" H 4450 2150 50  0001 C CNN
F 3 "" H 4450 2150 50  0001 C CNN
	1    4450 2150
	1    0    0    -1  
$EndComp
Wire Wire Line
	4250 2450 4250 2550
Connection ~ 4250 2550
Wire Wire Line
	4250 3050 4250 2950
Connection ~ 4250 2950
Wire Wire Line
	4250 3350 4250 3450
Wire Wire Line
	4250 2150 4250 2050
Wire Wire Line
	4250 2050 4450 2050
Wire Wire Line
	4450 2050 4450 2150
Wire Wire Line
	3400 2550 3550 2550
Wire Wire Line
	3250 2750 3250 2950
$EndSCHEMATC
