/******************************************************************************
*
* (c) Copyright 2010-2013 Xilinx, Inc. All rights reserved.
*
* This file contains confidential and proprietary information of Xilinx, Inc.
* and is protected under U.S. and international copyright and other
* intellectual property laws.
*
* DISCLAIMER
* This disclaimer is not a license and does not grant any rights to the
* materials distributed herewith. Except as otherwise provided in a valid
* license issued to you by Xilinx, and to the maximum extent permitted by
* applicable law: (1) THESE MATERIALS ARE MADE AVAILABLE "AS IS" AND WITH ALL
* FAULTS, AND XILINX HEREBY DISCLAIMS ALL WARRANTIES AND CONDITIONS, EXPRESS,
* IMPLIED, OR STATUTORY, INCLUDING BUT NOT LIMITED TO WARRANTIES OF
* MERCHANTABILITY, NON-INFRINGEMENT, OR FITNESS FOR ANY PARTICULAR PURPOSE;
* and (2) Xilinx shall not be liable (whether in contract or tort, including
* negligence, or under any other theory of liability) for any loss or damage
* of any kind or nature related to, arising under or in connection with these
* materials, including for any direct, or any indirect, special, incidental,
* or consequential loss or damage (including loss of data, profits, goodwill,
* or any type of loss or damage suffered as a result of any action brought by
* a third party) even if such damage or loss was reasonably foreseeable or
* Xilinx had been advised of the possibility of the same.
*
* CRITICAL APPLICATIONS
* Xilinx products are not designed or intended to be fail-safe, or for use in
* any application requiring fail-safe performance, such as life-support or
* safety devices or systems, Class III medical devices, nuclear facilities,
* applications related to the deployment of airbags, or any other applications
* that could lead to death, personal injury, or severe property or
* environmental damage (individually and collectively, "Critical
* Applications"). Customer assumes the sole risk and liability of any use of
* Xilinx products in Critical Applications, subject only to applicable laws
* and regulations governing limitations on product liability.
*
* THIS COPYRIGHT NOTICE AND DISCLAIMER MUST BE RETAINED AS PART OF THIS FILE
* AT ALL TIMES.
*
******************************************************************************/
/*****************************************************************************/
/**
 *
 * @file xaxidma_example_sg_poll.c
 *
 * This file demonstrates how to use the xaxidma driver on the Xilinx AXI
 * DMA core (AXIDMA) to transfer packets in polling mode when the AXIDMA
 * core is configured in Scatter Gather Mode.
 *
 * This code assumes a loopback hardware widget is connected to the AXI DMA
 * core for data packet loopback.
 *
 * To see the debug print, you need a Uart16550 or uartlite in your system,
 * and please set "-DDEBUG" in your compiler options. You need to rebuild your
 * software executable.
 *
 * Make sure that MEMORY_BASE is defined properly as per the HW system. The
 * h/w system built in Area mode has a maximum DDR memory limit of 64MB. In
 * throughput mode, it is 512MB.  These limits are need to ensured for
 * proper operation of this code.
 *
 *
 * <pre>
 * MODIFICATION HISTORY:
 *
 * Ver   Who  Date     Changes
 * ----- ---- -------- -------------------------------------------------------
 * 1.00a jz   05/17/10 First release
 * 2.00a jz   08/10/10 Second release, added in xaxidma_g.c, xaxidma_sinit.c,
 *                     updated tcl file, added xaxidma_porting_guide.h, removed
 *                     workaround for endianness
 * 4.00a rkv  02/22/11 Name of the file has been changed for naming consistency
 *       	       	   Added interrupt support for ARM.
 * 5.00a srt  03/05/12 Added Flushing and Invalidation of Caches to fix CRs
 *		       		   648103, 648701.
 *		       		   Added V7 DDR Base Address to fix CR 649405.
 * 6.00a srt  03/27/12 Changed API calls to support MCDMA driver.
 * 7.00a srt  06/18/12 API calls are reverted back for backward compatibility.
 * 7.01a srt  11/02/12 Buffer sizes (Tx and Rx) are modified to meet maximum
 *		       DDR memory limit of the h/w system built with Area mode
 * 7.02a srt  03/01/13 Updated DDR base address for IPI designs (CR 703656).
 *
 * </pre>
 *
 * ***************************************************************************
 */
/***************************** Include Files *********************************/
#include <stdio.h>
#include <stdlib.h>
#include "xaxidma.h"
#include "xparameters.h"
#include "xdebug.h"
#include "platform.h"
#include "cplx_data.h"
#include "stim.h"
//#include "xgpio.h"

#if (!defined(DEBUG))
extern void xil_printf(const char *format, ...);
#endif

/******************** Constant Definitions **********************************/

/*
 * Device hardware build related constants.
 */

#define DMA_DEV_ID		XPAR_AXIDMA_0_DEVICE_ID

#define MEM_BASE_ADDR		0x10000000

#define TX_BD_SPACE_BASE	(MEM_BASE_ADDR)
#define TX_BD_SPACE_HIGH	(MEM_BASE_ADDR + 0x00000FFF)
#define RX_BD_SPACE_BASE	(MEM_BASE_ADDR + 0x00001000)
#define RX_BD_SPACE_HIGH	(MEM_BASE_ADDR + 0x00001FFF)
#define TX_BUFFER_BASE		(MEM_BASE_ADDR + 0x00100000)
#define RX_BUFFER_BASE		(MEM_BASE_ADDR + 0x00300000)
#define RX_BUFFER_HIGH		(MEM_BASE_ADDR + 0x004FFFFF)

#define MEM_BASE_ADDR2		0x11000000

#define TX_BD_SPACE_BASE2	(MEM_BASE_ADDR2)
#define TX_BD_SPACE_HIGH2	(MEM_BASE_ADDR2 + 0x00000FFF)
#define RX_BD_SPACE_BASE2	(MEM_BASE_ADDR2 + 0x00001000)
#define RX_BD_SPACE_HIGH2	(MEM_BASE_ADDR2 + 0x00001FFF)
#define TX_BUFFER_BASE2		(MEM_BASE_ADDR2 + 0x00100000)
#define RX_BUFFER_BASE2		(MEM_BASE_ADDR2 + 0x00300000)
#define RX_BUFFER_HIGH2		(MEM_BASE_ADDR2 + 0x004FFFFF)

#define MEM_BASE_ADDR3		0x12000000

#define TX_BD_SPACE_BASE3	(MEM_BASE_ADDR3)
#define TX_BD_SPACE_HIGH3	(MEM_BASE_ADDR3 + 0x00000FFF)
#define RX_BD_SPACE_BASE3	(MEM_BASE_ADDR3 + 0x00001000)
#define RX_BD_SPACE_HIGH3	(MEM_BASE_ADDR3 + 0x00001FFF)
#define TX_BUFFER_BASE3		(MEM_BASE_ADDR3 + 0x00100000)
#define RX_BUFFER_BASE3		(MEM_BASE_ADDR3 + 0x00300000)
#define RX_BUFFER_HIGH3		(MEM_BASE_ADDR3 + 0x004FFFFF)

#define FFT_MAX_NUM_PTS      8192
#define MAX_PKT_LEN		0x20

#define TEST_START_VALUE	0xC

// External data Recieve From other file
extern int sig_two_sine_waves[FFT_MAX_NUM_PTS]; // FFT input data


/**************************** Type Definitions *******************************/


/***************** Macros (Inline Functions) Definitions *********************/


/************************** Function Prototypes ******************************/
static int RxSetup(XAxiDma * AxiDmaInstPtr);
static int TxSetup(XAxiDma * AxiDmaInstPtr);
static int SendPacket(XAxiDma * AxiDmaInstPtr,int c);
static int CheckData(void);
static int CheckDmaResult(XAxiDma * AxiDmaInstPtr);
static void Setdata(void);

/************************** Variable Definitions *****************************/
/*
 * Device instance definitions
 */
XAxiDma AxiDma;
XAxiDma AxiDma2;
XAxiDma AxiDma3;
//XGpio	*gpio_inst;
/*
 * Buffer for transmit packet. Must be 32-bit aligned to be used by DMA.
 */
u32 *Sig = (u32 *) TX_BUFFER_BASE;
u32 *SigRef = (u32 *) TX_BUFFER_BASE2;
u32 *Disp = (u32 *) TX_BUFFER_BASE3;

/*****************************************************************************/
/**
*
* Main function
*
* This function is the main entry of the tests on DMA core. It sets up
* DMA engine to be ready to receive and send packets, then a packet is
* transmitted and will be verified after it is received via the DMA loopback
* widget.
*
* @param	None
*
* @return
*		- XST_SUCCESS if test passes
*		- XST_FAILURE if test fails.
*
* @note		None.
*
******************************************************************************/
int main(void)
{

	int Status;
	XAxiDma_Config *Config;
	XAxiDma_Config *Config2;
	XAxiDma_Config *Config3;

	init_platform();
	    // Fill stimulus buffer with some signal
int i;
	for (i = 0; i < 2048; ++i) {
if(i <= 1023){
	Sig[i] = (u32)(sigL_img[i] & 1023);
	Sig[i] = Sig[i]<<10 | (u32)(sigL_real[i] & 1023);
	Disp[i] = (u32)(Dispa[i] & 262143);

}
SigRef[i] = (u32)(sigR_img[i] & 1023);
SigRef[i] = SigRef[i]<<10 | (u32)(sigR_Real[i] & 1023);


	}

	   // memcpy(Sig, sig_two_sine_waves, sizeof(u32)*FFT_MAX_NUM_PTS);
	    //memcpy(SigRef, sig_two_sine_waves2, sizeof(u32)*FFT_MAX_NUM_PTS);
	    //memcpy(Disp, sig_two_sine_waves3, sizeof(u32)*FFT_MAX_NUM_PTS);
	xil_printf("\r\n--- Entering main() --- \r\n");

	Config = XAxiDma_LookupConfig(DMA_DEV_ID);
	if (!Config) {
		xil_printf("No config found for %d\r\n", DMA_DEV_ID);

		return XST_FAILURE;
	}
	Config2 = XAxiDma_LookupConfig(XPAR_AXIDMA_1_DEVICE_ID);
		if (!Config) {
			xil_printf("No config found for %d\r\n", DMA_DEV_ID);

			return XST_FAILURE;
		}
		Config3 = XAxiDma_LookupConfig(XPAR_AXIDMA_2_DEVICE_ID);
			if (!Config) {
				xil_printf("No config found for %d\r\n", DMA_DEV_ID);

				return XST_FAILURE;
			}


	/* Initialize DMA engine */
	Status = XAxiDma_CfgInitialize(&AxiDma, Config);
	if (Status != XST_SUCCESS) {
		xil_printf("Initialization failed %d\r\n", Status);
		return XST_FAILURE;
	}
	Status = XAxiDma_CfgInitialize(&AxiDma2, Config2);
		if (Status != XST_SUCCESS) {
			xil_printf("Initialization failed %d\r\n", Status);
			return XST_FAILURE;
		}
		Status = XAxiDma_CfgInitialize(&AxiDma3, Config3);
			if (Status != XST_SUCCESS) {
				xil_printf("Initialization failed %d\r\n", Status);
				return XST_FAILURE;
			}

	if(!XAxiDma_HasSg(&AxiDma)) {
		xil_printf("Device configured as Simple mode \r\n");

		return XST_FAILURE;
	}
	if(!XAxiDma_HasSg(&AxiDma2)) {
			xil_printf("Device configured as Simple mode \r\n");

			return XST_FAILURE;
		}
	if(!XAxiDma_HasSg(&AxiDma3)) {
			xil_printf("Device configured as Simple mode \r\n");

			return XST_FAILURE;
		}

	Status = TxSetup(&AxiDma);
	if (Status != XST_SUCCESS) {
		return XST_FAILURE;
	}


	Status = RxSetup(&AxiDma);
	if (Status != XST_SUCCESS) {
		return XST_FAILURE;
	}

	xil_printf("checkpoint 1\n\r");

	/* Send a packet */
	int sendC =0;
	for (sendC = 0; sendC < 1; sendC++) {
		SendPacket(&AxiDma,sendC);
		if (Status != XST_SUCCESS) {
				return XST_FAILURE;
			}

	}


//	xil_printf("checkpoint 2 %d",XPAR_AXI_GPIO_0_DEVICE_ID);
//	XGpio_Initialize(gpio_inst, 2);
//	if (Status != XST_SUCCESS)
//		{
//			xil_printf("ERROR! Initialization of AXI GPIO instance failed.\n\r");
//
//		}
//	int reg= 1;
//	xil_printf("checkpoint 3");
//	XGpio_DiscreteWrite(gpio_inst, 1, reg);
//	//XGpio_DiscreteWrite()

xil_printf("\nfinish!!!!!!!");
	//Status = SendPacket(&AxiDma);

	/* Check DMA transfer result */
	Status = CheckDmaResult(&AxiDma);

	xil_printf("AXI DMA SG Polling Test %s\r\n",
		(Status == XST_SUCCESS)? "passed":"failed");

	xil_printf("--- Exiting main() --- \r\n");

	if (Status != XST_SUCCESS) {
		return XST_FAILURE;
	}

	return XST_SUCCESS;
}

/*****************************************************************************/
/**
*
* This function sets up RX channel of the DMA engine to be ready for packet
* reception
*
* @param	AxiDmaInstPtr is the pointer to the instance of the DMA engine.
*
* @return	XST_SUCCESS if the setup is successful, XST_FAILURE otherwise.
*
* @note		None.
*
******************************************************************************/
static int RxSetup(XAxiDma * AxiDmaInstPtr)
{
	XAxiDma_BdRing *RxRingPtr;
	XAxiDma_BdRing *RxRingPtr2;
	XAxiDma_BdRing *RxRingPtr3;
	int Delay = 0;
	int Coalesce = 1;
	int Status;
	XAxiDma_Bd BdTemplate;
	XAxiDma_Bd BdTemplate2;
	XAxiDma_Bd BdTemplate3;
	XAxiDma_Bd *BdPtr;
	XAxiDma_Bd *BdPtr2;
	XAxiDma_Bd *BdPtr3;
	XAxiDma_Bd *BdCurPtr;
	XAxiDma_Bd *BdCurPtr2;
	XAxiDma_Bd *BdCurPtr3;
	u32 BdCount;
	u32 BdCount2;
	u32 BdCount3;
	u32 FreeBdCount;
	u32 FreeBdCount2;
	u32 FreeBdCount3;
	u32 RxBufferPtr;
	u32 RxBufferPtr2;
	u32 RxBufferPtr3;
	int Index;

	RxRingPtr = XAxiDma_GetRxRing(&AxiDma);
	RxRingPtr2 = XAxiDma_GetRxRing(&AxiDma2);
	RxRingPtr3 = XAxiDma_GetRxRing(&AxiDma3);
	/* Disable all RX interrupts before RxBD space setup */

	XAxiDma_BdRingIntDisable(RxRingPtr, XAXIDMA_IRQ_ALL_MASK);
	XAxiDma_BdRingIntDisable(RxRingPtr2, XAXIDMA_IRQ_ALL_MASK);
	XAxiDma_BdRingIntDisable(RxRingPtr3, XAXIDMA_IRQ_ALL_MASK);
	/* Set delay and coalescing */
	XAxiDma_BdRingSetCoalesce(RxRingPtr, Coalesce, Delay);
	XAxiDma_BdRingSetCoalesce(RxRingPtr2, Coalesce, Delay);
	XAxiDma_BdRingSetCoalesce(RxRingPtr3, Coalesce, Delay);
	/* Setup Rx BD space */
	BdCount = XAxiDma_BdRingCntCalc(XAXIDMA_BD_MINIMUM_ALIGNMENT,
				RX_BD_SPACE_HIGH - RX_BD_SPACE_BASE + 1);
	BdCount2 = XAxiDma_BdRingCntCalc(XAXIDMA_BD_MINIMUM_ALIGNMENT,
				RX_BD_SPACE_HIGH2 - RX_BD_SPACE_BASE2 + 1);
	BdCount3 = XAxiDma_BdRingCntCalc(XAXIDMA_BD_MINIMUM_ALIGNMENT,
				RX_BD_SPACE_HIGH3 - RX_BD_SPACE_BASE3 + 1);

	Status = XAxiDma_BdRingCreate(RxRingPtr, RX_BD_SPACE_BASE,
				RX_BD_SPACE_BASE,
				XAXIDMA_BD_MINIMUM_ALIGNMENT, BdCount);
	Status = XAxiDma_BdRingCreate(RxRingPtr2, RX_BD_SPACE_BASE2,
				RX_BD_SPACE_BASE2,
				XAXIDMA_BD_MINIMUM_ALIGNMENT, BdCount2);
	Status = XAxiDma_BdRingCreate(RxRingPtr3, RX_BD_SPACE_BASE3,
				RX_BD_SPACE_BASE3,
				XAXIDMA_BD_MINIMUM_ALIGNMENT, BdCount3);

	if (Status != XST_SUCCESS) {
		xil_printf("RX create BD ring failed %d\r\n", Status);

		return XST_FAILURE;
	}

	/*
	 * Setup an all-zero BD as the template for the Rx channel.
	 */
	XAxiDma_BdClear(&BdTemplate);
	XAxiDma_BdClear(&BdTemplate2);
	XAxiDma_BdClear(&BdTemplate3);

	Status = XAxiDma_BdRingClone(RxRingPtr, &BdTemplate);
	Status = XAxiDma_BdRingClone(RxRingPtr2, &BdTemplate2);
	Status = XAxiDma_BdRingClone(RxRingPtr3, &BdTemplate3);
	if (Status != XST_SUCCESS) {
		xil_printf("RX clone BD failed %d\r\n", Status);

		return XST_FAILURE;
	}

	/* Attach buffers to RxBD ring so we are ready to receive packets */

	FreeBdCount = XAxiDma_BdRingGetFreeCnt(RxRingPtr);
	FreeBdCount2 = XAxiDma_BdRingGetFreeCnt(RxRingPtr2);
	FreeBdCount3 = XAxiDma_BdRingGetFreeCnt(RxRingPtr3);
	xil_printf("RxFreeBdCount %d  %d  %d", FreeBdCount, FreeBdCount2, FreeBdCount3);
	Status = XAxiDma_BdRingAlloc(RxRingPtr, FreeBdCount, &BdPtr);
	Status = XAxiDma_BdRingAlloc(RxRingPtr2, FreeBdCount2, &BdPtr2);
	Status = XAxiDma_BdRingAlloc(RxRingPtr3, FreeBdCount3, &BdPtr3);
	if (Status != XST_SUCCESS) {
		xil_printf("RX alloc BD failed %d\r\n", Status);

		return XST_FAILURE;
	}
	//int i = 0;
	BdCurPtr = BdPtr;
	BdCurPtr2 = BdPtr2;
	BdCurPtr3 = BdPtr3;
	RxBufferPtr = RX_BUFFER_BASE;
	RxBufferPtr2 = RX_BUFFER_BASE2;
	RxBufferPtr3 = RX_BUFFER_BASE3;
	for (Index = 0; Index < FreeBdCount; Index++) {
		Status = XAxiDma_BdSetBufAddr(BdCurPtr, RxBufferPtr);
		Status = XAxiDma_BdSetBufAddr(BdCurPtr2, RxBufferPtr2);
		Status = XAxiDma_BdSetBufAddr(BdCurPtr3, RxBufferPtr3);
		if (Status != XST_SUCCESS) {
			xil_printf("Set buffer addr %x on BD %x failed %d\r\n",
			    (unsigned int)RxBufferPtr,
			    (unsigned int)BdCurPtr, Status);

			return XST_FAILURE;
		}
// test_fix
		Status = XAxiDma_BdSetLength(BdCurPtr, MAX_PKT_LEN*2,
				RxRingPtr->MaxTransferLen);
		Status = XAxiDma_BdSetLength(BdCurPtr2, MAX_PKT_LEN*2,
						RxRingPtr2->MaxTransferLen);
		Status = XAxiDma_BdSetLength(BdCurPtr3, MAX_PKT_LEN*2,
						RxRingPtr3->MaxTransferLen);
		if (Status != XST_SUCCESS) {
			xil_printf("Rx set length %d on BD %x failed %d\r\n",
			    MAX_PKT_LEN*2, (unsigned int)BdCurPtr, Status);

			return XST_FAILURE;
		}

		/* Receive BDs do not need to set anything for the control
		 * The hardware will set the SOF/EOF bits per stream status
		 */
		XAxiDma_BdSetCtrl(BdCurPtr, 0);
		XAxiDma_BdSetCtrl(BdCurPtr2, 0);
		XAxiDma_BdSetCtrl(BdCurPtr3, 0);
		XAxiDma_BdSetId(BdCurPtr, RxBufferPtr);
		XAxiDma_BdSetId(BdCurPtr2, RxBufferPtr2);
		XAxiDma_BdSetId(BdCurPtr3, RxBufferPtr3);

		RxBufferPtr += MAX_PKT_LEN*2;
		RxBufferPtr2 += MAX_PKT_LEN*2;
		RxBufferPtr3 += MAX_PKT_LEN*2;
		BdCurPtr = (XAxiDma_Bd *)((void *)XAxiDma_BdRingNext(RxRingPtr, BdCurPtr));
		BdCurPtr2 = (XAxiDma_Bd *)((void *)XAxiDma_BdRingNext(RxRingPtr2, BdCurPtr2));
		BdCurPtr3 = (XAxiDma_Bd *)((void *)XAxiDma_BdRingNext(RxRingPtr3, BdCurPtr3));
	}

	/* Clear the receive buffer, so we can verify data
	 */
	memset((void *)RX_BUFFER_BASE, 0, MAX_PKT_LEN*128);
	memset((void *)RX_BUFFER_BASE2, 0, MAX_PKT_LEN*128);
	memset((void *)RX_BUFFER_BASE3, 0, MAX_PKT_LEN*128);
//	u32 *RxPacket = RX_BUFFER_BASE;
//					int t = 0;
//					for (t = 0; t <16*32 ; t++) {
//						xil_printf("num[%d] = %d \n\r",i,RxPacket[t]);
//						i++;
//					}
	Status = XAxiDma_BdRingToHw(RxRingPtr, FreeBdCount,
						BdPtr);
	if (Status != XST_SUCCESS) {
			xil_printf("RX submit hw failed %d\r\n", Status);

			return XST_FAILURE;
		}
	Status = XAxiDma_BdRingToHw(RxRingPtr2, FreeBdCount2,
							BdPtr2);
	if (Status != XST_SUCCESS) {
			xil_printf("RX// submit hw failed %d\r\n", Status);

			return XST_FAILURE;
		}
	Status = XAxiDma_BdRingToHw(RxRingPtr3, FreeBdCount3,
							BdPtr3);
	if (Status != XST_SUCCESS) {
		xil_printf("RX/// submit hw failed %d\r\n", Status);

		return XST_FAILURE;
	}

	/* Start RX DMA channel */
	Status = XAxiDma_BdRingStart(RxRingPtr);
	Status = XAxiDma_BdRingStart(RxRingPtr2);
	Status = XAxiDma_BdRingStart(RxRingPtr3);
	if (Status != XST_SUCCESS) {
		xil_printf("RX start hw failed %d\r\n", Status);

		return XST_FAILURE;
	}

	return XST_SUCCESS;
}

/*****************************************************************************/
/**
*
* This function sets up the TX channel of a DMA engine to be ready for packet
* transmission
*
* @param	AxiDmaInstPtr is the instance pointer to the DMA engine.
*
* @return	XST_SUCCESS if the setup is successful, XST_FAILURE otherwise.
*
* @note		None.
*
******************************************************************************/
static int TxSetup(XAxiDma * AxiDmaInstPtr)
{

	XAxiDma_BdRing *TxRingPtr;
	XAxiDma_Bd BdTemplate;
	XAxiDma_BdRing *TxRingPtr2;
	XAxiDma_Bd BdTemplate2;
	XAxiDma_BdRing *TxRingPtr3;
	XAxiDma_Bd BdTemplate3;

	int Delay = 0;
	int Coalesce = 1;
	int Status;
	u32 BdCount;
	u32 BdCount2;
	u32 BdCount3;

	TxRingPtr = XAxiDma_GetTxRing(&AxiDma);
	TxRingPtr2 = XAxiDma_GetTxRing(&AxiDma2);
	TxRingPtr3 = XAxiDma_GetTxRing(&AxiDma3);
	/* Disable all TX interrupts before TxBD space setup */

	XAxiDma_BdRingIntDisable(TxRingPtr, XAXIDMA_IRQ_ALL_MASK);
	XAxiDma_BdRingIntDisable(TxRingPtr2, XAXIDMA_IRQ_ALL_MASK);
	XAxiDma_BdRingIntDisable(TxRingPtr3, XAXIDMA_IRQ_ALL_MASK);
	/* Set TX delay and coalesce */
	XAxiDma_BdRingSetCoalesce(TxRingPtr, Coalesce, Delay);
	XAxiDma_BdRingSetCoalesce(TxRingPtr2, Coalesce, Delay);
	XAxiDma_BdRingSetCoalesce(TxRingPtr3, Coalesce, Delay);

	/* Setup TxBD space  */
	BdCount = XAxiDma_BdRingCntCalc(XAXIDMA_BD_MINIMUM_ALIGNMENT,
				TX_BD_SPACE_HIGH - TX_BD_SPACE_BASE + 1);
	BdCount2 = XAxiDma_BdRingCntCalc(XAXIDMA_BD_MINIMUM_ALIGNMENT,
				TX_BD_SPACE_HIGH2 - TX_BD_SPACE_BASE2 + 1);
	BdCount3 = XAxiDma_BdRingCntCalc(XAXIDMA_BD_MINIMUM_ALIGNMENT,
				TX_BD_SPACE_HIGH3 - TX_BD_SPACE_BASE3 + 1);


	Status = XAxiDma_BdRingCreate(TxRingPtr, TX_BD_SPACE_BASE,
				TX_BD_SPACE_BASE,
				XAXIDMA_BD_MINIMUM_ALIGNMENT, BdCount);
	Status = XAxiDma_BdRingCreate(TxRingPtr2, TX_BD_SPACE_BASE2,
				TX_BD_SPACE_BASE2,
				XAXIDMA_BD_MINIMUM_ALIGNMENT, BdCount2);
	Status = XAxiDma_BdRingCreate(TxRingPtr3, TX_BD_SPACE_BASE3,
				TX_BD_SPACE_BASE3,
				XAXIDMA_BD_MINIMUM_ALIGNMENT, BdCount3);
	if (Status != XST_SUCCESS) {
		xil_printf("failed create BD ring in txsetup\r\n");

		return XST_FAILURE;
	}
//xil_printf("BdCount = %d\r\n",BdCount);
	/*
	 * We create an all-zero BD as the template.
	 */
	XAxiDma_BdClear(&BdTemplate);
	XAxiDma_BdClear(&BdTemplate2);
	XAxiDma_BdClear(&BdTemplate3);

	Status = XAxiDma_BdRingClone(TxRingPtr, &BdTemplate);
	Status = XAxiDma_BdRingClone(TxRingPtr2, &BdTemplate2);
	Status = XAxiDma_BdRingClone(TxRingPtr3, &BdTemplate3);
	if (Status != XST_SUCCESS) {
		xil_printf("failed bdring clone in txsetup %d\r\n", Status);

		return XST_FAILURE;
	}

	/* Start the TX channel */
	Status = XAxiDma_BdRingStart(TxRingPtr);
	Status = XAxiDma_BdRingStart(TxRingPtr2);
	Status = XAxiDma_BdRingStart(TxRingPtr3);
	if (Status != XST_SUCCESS) {
		xil_printf("failed start bdring txsetup %d\r\n", Status);

		return XST_FAILURE;
	}

	return XST_SUCCESS;
}

/*****************************************************************************/
/**
*
* This function transmits one packet non-blockingly through the DMA engine.
*
* @param	AxiDmaInstPtr points to the DMA engine instance
*
* @return	- XST_SUCCESS if the DMA accepts the packet successfully,
*		- XST_FAILURE otherwise.
*
* @note     None.
*
******************************************************************************/
static int SendPacket(XAxiDma * AxiDmaInstPtr,int c)
{
	XAxiDma_BdRing *TxRingPtr; //TxRing
	XAxiDma_BdRing *TxRingPtr2; //TxRing
	XAxiDma_BdRing *TxRingPtr3; //TxRing
	u32 *TxPacket;	   //PacketPointer
	u32 *TxPacket2;
	u32 *TxPacket3;

	XAxiDma_Bd *BdPtr;			//Start BDPointer
	XAxiDma_Bd *BdPtr2;
	XAxiDma_Bd *BdPtr3;
	int Status;
	int Index;
	TxRingPtr = XAxiDma_GetTxRing(&AxiDma);
	TxRingPtr2 = XAxiDma_GetTxRing(&AxiDma2);
	TxRingPtr3 = XAxiDma_GetTxRing(&AxiDma3);


	TxPacket = 	(u32 *)TX_BUFFER_BASE;	//setPointerForPacket
	TxPacket2 = 	(u32 *)TX_BUFFER_BASE2;
	TxPacket3 = 	(u32 *)TX_BUFFER_BASE3;
	/* Flush the SrcBuffer before the DMA transfer, in case the Data Cache
	 * is enabled
	 */
	Xil_DCacheFlushRange((u32)TxPacket, MAX_PKT_LEN*128); //8 point each packet
	Xil_DCacheFlushRange((u32)TxPacket2, MAX_PKT_LEN*256);
	Xil_DCacheFlushRange((u32)TxPacket3, MAX_PKT_LEN*128);


	int FreeBdCount = XAxiDma_BdRingGetFreeCnt(TxRingPtr);
	int FreeBdCount2 = XAxiDma_BdRingGetFreeCnt(TxRingPtr2);
	int FreeBdCount3 = XAxiDma_BdRingGetFreeCnt(TxRingPtr3);
		xil_printf("TxFreeBdCount %d  %d  %d", FreeBdCount ,FreeBdCount2 ,FreeBdCount3);
	/* Allocate a BD */
	Status = XAxiDma_BdRingAlloc(TxRingPtr, FreeBdCount, &BdPtr);
	Status = XAxiDma_BdRingAlloc(TxRingPtr2, FreeBdCount2, &BdPtr2);
	Status = XAxiDma_BdRingAlloc(TxRingPtr3, FreeBdCount3, &BdPtr3);

	if (Status != XST_SUCCESS) {
		return XST_FAILURE;
	}
	xil_printf("check1\n\r");
	XAxiDma_Bd *BdCurPtr;
	XAxiDma_Bd *BdCurPtr2;
	XAxiDma_Bd *BdCurPtr3;
	u32 TxBufferPtr;
	u32 TxBufferPtr2;
	u32 TxBufferPtr3;

	BdCurPtr = BdPtr;		// Set Current BDptr
	BdCurPtr2 = BdPtr2;
	BdCurPtr3 = BdPtr3;
	TxBufferPtr = (u32 )TX_BUFFER_BASE;
	TxBufferPtr2 = (u32 )TX_BUFFER_BASE2;
	TxBufferPtr3 = (u32 )TX_BUFFER_BASE3;
	int i = 0;
	for (Index = 0; Index < FreeBdCount; Index++) {
		TxPacket = (u32*)TxBufferPtr;
		int t = 0;
		for (t = 0; t <16 ; t++) {
			xil_printf("num[%d] = %d  %d  %d\n\r",i,TxPacket[t],TxPacket2[t],TxPacket3[t]);
			i++;
		}


		Status = XAxiDma_BdSetBufAddr(BdCurPtr, TxBufferPtr);
		Status = XAxiDma_BdSetBufAddr(BdCurPtr2, TxBufferPtr2);
		Status = XAxiDma_BdSetBufAddr(BdCurPtr3, TxBufferPtr3);
		if (Status != XST_SUCCESS) {
		xil_printf("Tx set buffer addr %x on BD %x failed %d\r\n",
		    (unsigned int)TxBufferPtr, (unsigned int)BdCurPtr, Status);

		return XST_FAILURE;
	}


	Status = XAxiDma_BdSetLength(BdCurPtr, MAX_PKT_LEN*4,
				TxRingPtr->MaxTransferLen);
	Status = XAxiDma_BdSetLength(BdCurPtr2, MAX_PKT_LEN*8,
				TxRingPtr2->MaxTransferLen);
	Status = XAxiDma_BdSetLength(BdCurPtr3, MAX_PKT_LEN*4,
				TxRingPtr3->MaxTransferLen);
	if (Status != XST_SUCCESS) {
		xil_printf("Tx set length %d on BD %x failed %d\r\n",
		    MAX_PKT_LEN*2, (unsigned int)BdCurPtr, Status);

		return XST_FAILURE;
	}
	#if (XPAR_AXIDMA_0_SG_INCLUDE_STSCNTRL_STRM == 1)
	Status = XAxiDma_BdSetAppWord(BdCurPtr,
	    XAXIDMA_LAST_APPWORD, MAX_PKT_LEN*2);
	Status = XAxiDma_BdSetAppWord(BdCurPtr2,
		    XAXIDMA_LAST_APPWORD, MAX_PKT_LEN*4);
	Status = XAxiDma_BdSetAppWord(BdCurPtr3,
		    XAXIDMA_LAST_APPWORD, MAX_PKT_LEN*2);

	/* If Set app length failed, it is not fatal
	 */
	if (Status != XST_SUCCESS) {
		xil_printf("Set app word failed with %d\r\n", Status);
	}
	#endif

	/* For single packet, both SOF and EOF are to be set
	 */
	XAxiDma_BdSetCtrl(BdCurPtr, XAXIDMA_BD_CTRL_TXEOF_MASK |
						XAXIDMA_BD_CTRL_TXSOF_MASK);
	XAxiDma_BdSetCtrl(BdCurPtr2, XAXIDMA_BD_CTRL_TXEOF_MASK |
						XAXIDMA_BD_CTRL_TXSOF_MASK);
	XAxiDma_BdSetCtrl(BdCurPtr3, XAXIDMA_BD_CTRL_TXEOF_MASK |
						XAXIDMA_BD_CTRL_TXSOF_MASK);
	XAxiDma_BdSetId(BdCurPtr, (u32) TxBufferPtr);
	XAxiDma_BdSetId(BdCurPtr2, (u32) TxBufferPtr2);
	XAxiDma_BdSetId(BdCurPtr3, (u32) TxBufferPtr3);
	TxBufferPtr += MAX_PKT_LEN*4;
	TxBufferPtr2 += MAX_PKT_LEN*8;
	TxBufferPtr3 += MAX_PKT_LEN*4;
	BdCurPtr = (XAxiDma_Bd *)((void *)XAxiDma_BdRingNext(TxRingPtr, BdCurPtr));
	BdCurPtr2 = (XAxiDma_Bd *)((void *)XAxiDma_BdRingNext(TxRingPtr2, BdCurPtr2));
	BdCurPtr3 = (XAxiDma_Bd *)((void *)XAxiDma_BdRingNext(TxRingPtr3, BdCurPtr3));
	}

	/* Give the BD to DMA to kick off the transmission. */
	Status = XAxiDma_BdRingToHw(TxRingPtr, FreeBdCount, BdPtr);
	Status = XAxiDma_BdRingToHw(TxRingPtr2, FreeBdCount2, BdPtr2);
	Status = XAxiDma_BdRingToHw(TxRingPtr3, FreeBdCount3, BdPtr3);
	if (Status != XST_SUCCESS) {

		xil_printf("to hw failed %d\r\n", Status);
		return XST_FAILURE;
	}
	return XST_SUCCESS;
}

/*****************************************************************************/
/*
*
* This function checks data buffer after the DMA transfer is finished.
*
* @param	None
*
* @return	- XST_SUCCESS if validation is successful
*		- XST_FAILURE if validation is failure.
*
* @note		None.
*
******************************************************************************/
static int CheckData(void)
{

	//cplx_data_t *RxPacket;
	u32 *RxPacket;
	int Index = 0;
	//char str[25];


	RxPacket =  (u32*)RX_BUFFER_BASE;
	//Value = TEST_START_VALUE;

	/* Invalidate the DestBuffer before receiving the data, in case the
	 * Data Cache is enabled
	 */
	Xil_DCacheInvalidateRange((u32)RxPacket, MAX_PKT_LEN*64);
	short tmp;
	int i = 0;
	for(Index = 0; Index < MAX_PKT_LEN*16; Index++) {
		RxPacket =  (u32*)(RX_BUFFER_BASE+i*4);
		if(*RxPacket>>19==1){
					 *RxPacket = (4095<<20)| (*RxPacket);
				}
		else *RxPacket = *RxPacket&1048575;
//		for (tmp = 31; tmp >=0; tmp--) {
//			if((tmp+1)%4==0) xil_printf(" ");
//				xil_printf("%d", (*RxPacket>>tmp)&1);
//			}
//			xil_printf("   No. %d \n\r",i);
			i++;
		xil_printf("%d \n\r",*RxPacket );

			//xil_printf("%d\n\r", ( unsigned int)(RxPacket[Index])&65535*2);
		//cplx_data_get_string(str, RxPacket[Index]);
		//xil_printf("xn(%d) = %s\n\r", Index, str);
//		if (RxPacket[Index] != Value) {
//			xil_printf("Data error %d: %x/%x\r\n",
//			    Index, (unsigned int)RxPacket[Index],
//			    (unsigned int)Value);
//
//			return XST_FAILURE;
		//}
		//Value = (Value + 1) & 0xFF;
	}

	return XST_SUCCESS;
}

/*****************************************************************************/
/**
*
* This function waits until the DMA transaction is finished, checks data,
* and cleans up.
*
* @param	None
*
* @return	- XST_SUCCESS if DMA transfer is successful and data is correct,
*		- XST_FAILURE if fails.
*
* @note		None.
*
******************************************************************************/
static int CheckDmaResult(XAxiDma * AxiDmaInstPtr)
{
	XAxiDma_BdRing *TxRingPtr;
	XAxiDma_BdRing *RxRingPtr;
	XAxiDma_BdRing *TxRingPtr2;
	XAxiDma_BdRing *RxRingPtr2;
	XAxiDma_BdRing *TxRingPtr3;
	XAxiDma_BdRing *RxRingPtr3;
	XAxiDma_Bd *BdPtr;
	XAxiDma_Bd *BdPtr2;
	XAxiDma_Bd *BdPtr3;
	int ProcessedBdCount;
	int FreeBdCount;
	int ProcessedBdCount2;
	int FreeBdCount2;
	int ProcessedBdCount3;
	int FreeBdCount3;
	int Status;

	TxRingPtr = XAxiDma_GetTxRing(&AxiDma);
	RxRingPtr = XAxiDma_GetRxRing(&AxiDma);
	TxRingPtr2 = XAxiDma_GetTxRing(&AxiDma2);
	RxRingPtr2 = XAxiDma_GetRxRing(&AxiDma2);
	TxRingPtr3 = XAxiDma_GetTxRing(&AxiDma3);
	RxRingPtr3 = XAxiDma_GetRxRing(&AxiDma3);
	/* Wait until the one BD TX transaction is done */
	while ((ProcessedBdCount = XAxiDma_BdRingFromHw(TxRingPtr,
						       XAXIDMA_ALL_BDS,
						       &BdPtr)) == 0) {
		xil_printf("wait for Transmit");
	}
	while ((ProcessedBdCount2 = XAxiDma_BdRingFromHw(TxRingPtr2,
							       XAXIDMA_ALL_BDS,
							       &BdPtr2)) == 0) {
			xil_printf("wait for Transmit2");
		}
	while ((ProcessedBdCount3 = XAxiDma_BdRingFromHw(TxRingPtr3,
							       XAXIDMA_ALL_BDS,
							       &BdPtr3)) == 0) {
			xil_printf("wait for Transmit3");
		}
		//xil_printf("Transmit ProcessedBdCount = %d\n",ProcessedBdCount);
		//xil_printf("Transmit ProcessedBdCount = %d\n",ProcessedBdCount);
	/* Free all processed TX BDs for future transmission */
	Status = XAxiDma_BdRingFree(TxRingPtr, ProcessedBdCount, BdPtr);
	Status = XAxiDma_BdRingFree(TxRingPtr2, ProcessedBdCount2, BdPtr2);
	Status = XAxiDma_BdRingFree(TxRingPtr3, ProcessedBdCount3, BdPtr3);
	if (Status != XST_SUCCESS) {
		xil_printf("Failed to free %d tx BDs %d\r\n",
		    ProcessedBdCount, Status);
		return XST_FAILURE;
	}
	xil_printf("CheckedFreeTxBd");

	/* Wait until the data has been received by the Rx channel */
	while ((ProcessedBdCount = XAxiDma_BdRingFromHw(RxRingPtr,
						       XAXIDMA_ALL_BDS,
						       &BdPtr)) == 0) {
		xil_printf("wait for Recieve");
	}
	xil_printf("Recieve ProcessedBdCount = %d\n\r",ProcessedBdCount);
	/* Check received data */
	if (CheckData() != XST_SUCCESS) {

		return XST_FAILURE;
	}

	/* Free all processed RX BDs for future transmission */
	Status = XAxiDma_BdRingFree(RxRingPtr, ProcessedBdCount, BdPtr);
	if (Status != XST_SUCCESS) {
		xil_printf("Failed to free %d rx BDs %d\r\n",
		    ProcessedBdCount, Status);
		return XST_FAILURE;
	}

	/* Return processed BDs to RX channel so we are ready to receive new
	 * packets:
	 *    - Allocate all free RX BDs
	 *    - Pass the BDs to RX channel
	 */
	FreeBdCount = XAxiDma_BdRingGetFreeCnt(RxRingPtr);
	Status = XAxiDma_BdRingAlloc(RxRingPtr, FreeBdCount, &BdPtr);
	if (Status != XST_SUCCESS) {
		xil_printf("bd alloc failed\r\n");
		return XST_FAILURE;
	}

	Status = XAxiDma_BdRingToHw(RxRingPtr, FreeBdCount, BdPtr);
	if (Status != XST_SUCCESS) {
		xil_printf("Submit %d rx BDs failed %d\r\n", FreeBdCount, Status);
		return XST_FAILURE;
	}
	return XST_SUCCESS;
}

static void Setdata(void){
//		int x = 1024;
//		int y = 1080;
//		int n = 3;
//		FILE * fi = fopen("raw_img.dat","rb");
//		printf("F : %x\n", fi);
//		char *data = (unsigned char*)malloc( x*y*n * sizeof(unsigned char));
//		fread( data, sizeof(char) , n*x*y, fi);
//		fclose(fi);
}
