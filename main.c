#include <msp430.h>
#include <stdlib.h>
#include <stdio.h>
#include <time.h>
#include <errno.h>
#include <HAR.h>
#include <string.h>

#define CS_8MHZ             0
#define CS_16MHZ            1

#define FILE_TRAIN_IMAGE        "C:\\Users\\kisho\\workspace_v9\\lenet2\\train-images-idx3-ubyte"
#define FILE_TRAIN_LABEL        "C:\\Users\\kisho\\workspace_v9\\lenet2\\train-labels-idx1-ubyte"
#define FILE_TEST_IMAGE     "C:\\Users\\kisho\\workspace_v9\\lenet2\\test_single_image"
#define FILE_TEST_LABEL     "C:\\Users\\kisho\\workspace_v9\\lenet2\\test_single_lebel"
#define LENET_FILE      "model.dat"
#define COUNT_TRAIN     60000
#define COUNT_TEST      1

short int i,percent,j,predict;




void configClocks(void)
{
#if CS_8MHZ
    //Ensure that there are 0 wait states enabled
    FRCTL0 = FRCTLPW | NWAITS_0;


    // Clock System Setup
    CSCTL0_H = CSKEY_H;                     // Unlock CS registers
    CSCTL1 = DCOFSEL_0;                     // Set DCO to 1MHz

    // Set SMCLK = MCLK = DCO, ACLK = VLOCLK
    CSCTL2 = SELA__VLOCLK | SELS__DCOCLK | SELM__DCOCLK;

    // Per Errata CS12 set divider to 4 before changing frequency to
    // prevent out of spec operation from overshoot transient
    CSCTL3 = DIVA__4 | DIVS__4 | DIVM__4;   // Set all corresponding clk sources to divide by 4 for errata
    CSCTL1 = DCOFSEL_6;                     // Set DCO to 8MHz

    // Delay by ~10us to let DCO settle. 60 cycles = 20 cycles buffer + (10us / (1/4MHz))
    __delay_cycles(60);

    CSCTL3 = DIVA__1 | DIVS__1 | DIVM__1;   // Set MCLK  to 8MHz, SMCLK to 1MHz

    CSCTL4 = HFXTOFF + LFXTOFF;             // Cancel external clocks
    CSCTL5 &= ~(HFXTOFFG + LFXTOFFG);

    CSCTL0_H = 0;                           // Lock CS Registers

#endif

#if CS_16MHZ
    // Configure one FRAM waitstate as required by the device datasheet for MCLK
    // operation beyond 8MHz _before_ configuring the clock system.
    FRCTL0 = FRCTLPW | NWAITS_1;

    // Clock System Setup
    CSCTL0_H = CSKEY_H;                     // Unlock CS registers
    CSCTL1 = DCOFSEL_0;                     // Set DCO to 1MHz

    // Set SMCLK = MCLK = DCO, ACLK = VLOCLK
    CSCTL2 = SELA__VLOCLK | SELS__DCOCLK | SELM__DCOCLK;

    // Per Errata CS12 set divider to 4 before changing frequency to
    // prevent out of spec operation from overshoot transient
    CSCTL3 = DIVA__4 | DIVS__4 | DIVM__4;   // Set all corresponding clk sources to divide by 4 for errata
    CSCTL1 = DCOFSEL_4 | DCORSEL;           // Set DCO to 16MHz

    // Delay by ~10us to let DCO settle. 60 cycles = 20 cycles buffer + (10us / (1/4MHz))
    __delay_cycles(60);

    CSCTL3 = DIVA__1 | DIVS__1 | DIVM__1;   // Set MCLK to 16MHz, SMCLK to 1MHz

    CSCTL4 = HFXTOFF + LFXTOFF;             // Cancel external clocks
    CSCTL5 &= ~(HFXTOFFG + LFXTOFFG);

    CSCTL0_H = 0;                           // Lock CS Registers
#endif
}


int main(void)
{
	WDTCTL = WDTPW | WDTHOLD;	// stop watchdog timer

	PM5CTL0 &= ~LOCKLPM5;

	configClocks();

	Predict();

    return 0;
}

