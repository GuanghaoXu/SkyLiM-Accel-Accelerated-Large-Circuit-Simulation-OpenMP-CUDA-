// skylim.cpp : This file contains the 'main' function. Program execution begins and ends there.
//
/*  Program to simulate large circuits
Started Tuesday February 9,1999, 5:40 PM
Copyright Jose Schutt-Aine
University of Illinois    */
// Modified to work on Visual Studio. 
// July 10, 2018

// Modified to use CUDA parallel computing by student Guanghao Xu
// Dec 2025

/* This program uses only one current source */
#define _CRT_SECURE_NO_WARNINGS

#include <stdio.h>
#include <math.h>
#include <time.h>
#include <string.h>
#include <stdlib.h>
#include <string.h>
#include <iostream>
#include <sys/types.h>
#include <sys/stat.h>
#include <chrono>  // [Added] include <chrono> for high-precision timing (no algorithm changes; timing only)

#include "lg.h"

// -----------------------------------------------------------------------------------------------------------
//                                 CUDA Header Starter
// -----------------------------------------------------------------------------------------------------------
#include <cuda_runtime.h>
// 引入 atomic functions
// #include <device_functions.h>
#define checkCudaErrors(val) check_cuda((val), #val, __FILE__, __LINE__)
// ------------------------------------ CUDA Header End ------------------------------------------------------

// -----------------------------------------------------------------------------------------------------------
//                                 CUDA Helper Function
// -----------------------------------------------------------------------------------------------------------
void check_cuda(cudaError_t result, char const* const func, const char* const file, int const line) {
	if (result) {
		fprintf(stderr, "CUDA error at %s:%d code=%d (%s) \"%s\" \n", file, line, static_cast<int>(result), cudaGetErrorString(result), func);
		exit(EXIT_FAILURE);
	}
}
// ---------------------------------- CUDA Helper Function END ---------------------------------------------
int* Nadi, * Nado, ** ingo, ** outgo;
branch* ckbrch;

// ---------------------------------------------------------------------------------------------------------
//                                 CUDA Kernel Definitions
// ---------------------------------------------------------------------------------------------------------

// Device version of comp function（Define before Kernel）
__device__ void comp_device(double t, double a, double ris, double fal, double wid, double* att)
{
	double b, c, d, p;
	b = a + ris;
	c = b + wid;
	d = c + fal;
	p = t;
	/* function is zero before a*/
	if (p < a) *att = 0;
	else if (p < b)  *att = 1 / (b - a) * (p - a);
	/* rising edge of pulse*/
	/* top of pulse (flat portion)*/
	else if (p < c) *att = 1;
	/* falling edge of pulse*/
	else if (p < d)  *att = 1 / (c - d) * (p - d);
	else  *att = 0;
	return;
}

// Branch Loop Kernel
__global__ void branch_kernel(branch* d_ckbrch, node* d_cknod, int Nb, double dt, int imethod) {
	int i = blockIdx.x * blockDim.x + threadIdx.x + 1;

	if (i <= Nb) {
		const double i_old = d_ckbrch[i].cur;
		const double l12 = d_ckbrch[i].lind;
		const double r12 = d_ckbrch[i].res;
		const int    n1 = d_ckbrch[i].n1;
		const int    n2 = d_ckbrch[i].n2;

		const double v1 = d_cknod[n1].vol;
		const double v2 = d_cknod[n2].vol;
		const double Gb = d_ckbrch[i].gb;
		const double Tb = d_ckbrch[i].tb;
		const double Gni = d_cknod[n1].gn, Tni = d_cknod[n1].tn;
		const double Gnj = d_cknod[n2].gn, Tnj = d_cknod[n2].tn;

		// Equivalent "snapshot" ixum
		const double ixum1 = d_cknod[n1].xum - i_old;
		const double ixum2 = d_cknod[n2].xum + i_old;

		double i_new;
		if (imethod == 1) {
			i_new = i_old + dt / l12 * (v1 - v2 - r12 * i_old);
		}
		else {
			i_new = (Tb * Gni * v1 - Tb * Tni * ixum1 - Tb * Gnj * v2 + Tb * Tnj * ixum2 + Gb * i_old)
				/ (1.0 + Tb * Tni + Tb * Tnj);
		}

		d_ckbrch[i].cur = i_new;

		// // Atomically update xum
		double delta1 = (-i_old + i_new);
		double delta2 = (+i_old - i_new);

		// Perform a thread-safe update using atomicAdd(ptr, val)
		atomicAdd(&(d_cknod[n1].xum), delta1);
		atomicAdd(&(d_cknod[n2].xum), delta2);
	}
}

// Node Loop Kernel
__global__ void node_kernel(node* d_cknod, int Nn, double tim, double dt, int imethod, int idvs, double magi, double del, double ris, double fal, double wid) {
	int i = blockIdx.x * blockDim.x + threadIdx.x + 1;

	if (i <= Nn) {
		const double g1 = d_cknod[i].gond;
		const double c1 = d_cknod[i].cap;
		const double v_old = d_cknod[i].vol;
		const double Gni = d_cknod[i].gn;
		const double Tni = d_cknod[i].tn;
		const double xum = d_cknod[i].xum;
		const int    nid = d_cknod[i].id;

		double zeti = 0.0;

		if (nid == idvs) {
			double zet;
			// Call the device version of the comp function
			comp_device(tim, del, ris, fal, wid, &zet);
			zeti = magi * zet;
		}

		const double isum = xum - zeti;

		double v_new;
		if (imethod == 1) {
			v_new = ((v_old * c1) / dt - isum) / (c1 / dt + g1);
		}
		else {
			v_new = Gni * v_old - Tni * isum;
		}

		d_cknod[i].vol = v_new;
	}
}

// -------------------------------------------- CUDA Kernel Definitions END --------------------------------------------

FILE* ofpg;
void stop();
void thexyplot(int mpl, double* timpl, double* vpl1, double* vpl2, double* vpl3, double* vpl4, char Mtitle[], char Xtitle[], char Ytitle[], char lg1[], char lg2[], char lg3[], char lg4[]);

int main(int argc, char* argv[])
{
	/* Declare variables */
	int i, j, k, idvs, t, tstart, tstop;
	int Nb, Nn, Ni;
	int mt, ps1, ps2, ps3, ps4;
	int kount, ivl, ovl;
	int nx1, nx2, fcnt, flag;
	int adjmax, imethod, kmaxx, ndinter, zer;
	int nodid, nodtyp, brid, brn1, brn2, brtyp, idumi;
	branch brch;
	node* cknod;
	node nod;
	time_t start, ending, td_init, td_begin, td_ended;
	double td_du;
	double ris, fal, wid, del, zet, sw, simstep, simdur;
	double dur, tim, dt, slc, ls, cs, zo;
	double i12, l12, r12, v1, v2, rtot, gtot;
	double inew, vnew, isum, g1, c1, magi, zeti;
	double iscp, ccp, gcp, rcp, vscp, lcp, esr, si, hsr;
	double nodgond, nodcap, brres, brlind;
	double ixum1, ixum2, dmaxsq, dmax, ii, ij;
	double* lps, * cps, hs1, per, vlow, vmag;
	double* vy1, * vy2, * vy3, * vy4, * timex;
	double Gb, Tb, Gni, Tni, Gnj, Tnj;
	double nodrg, gfict, cfict, brling, nodcag, nano, rshfict;
	FILE* ofp1, * ofp2, * inf1, * outf1, * outf2;

	// CUDA: Device pointers
	branch* d_ckbrch;
	node* d_cknod;
	size_t Nb_size;
	size_t Nn_size;

	// simulator monitors probe nodes (ps1, ps2, ps3, ps4) pinned Local cache
	// 1-based indexing,ignore index 0
	node* h_monitored_nodes = nullptr;
	checkCudaErrors(cudaMallocHost(&h_monitored_nodes, 5 * sizeof(node)));

	/* Open files for input and output */
	inf1 = fopen("lout_jsa", "r");
	outf1 = fopen("pdnjsa.cir", "w");
	//outf2 = fopen("lgan", "w");
	sw = 0.0;
	nano = 1.0e-09;
	fcnt = 10; /* scaling factor for numerical stability */
	kmaxx = 50;
	dmaxsq = 0.000100 * 0.250000;
	dmax = sqrt(dmaxsq);
	vlow = 0.0;
	vmag = 1.0;
	ris = 1; /* rise time */
	fal = 1; /* fall time */
	wid = 20.0; /* pulse width */
	del = 1.2; /* pulse delay */
	simstep = 0.02;
	simdur = 60.0;
	imethod = 2;
	adjmax = 200;
	idvs = 1; /* node where current source is located */
	dur = 3.0 * (del + ris + fal + wid); /* simulation duration */
	per = 1.5 * dur;

	Ni = 1;
	Nb = 169997; /* Number of branches */
	Nn = 100000; /* Number of nodes */
	mt = 20000 * fcnt; /* number of simulation time points */
	printf("MT is %10d \n", mt);
	dt = 50 * dur / mt; /* time step */
	magi = 0.02; /* magnitude of current source */
	printf("TSTEP IS: %.5e DMAX IS: %.5e \n", dt, dmax);
	//	scanf("%ld \n", &idumi);

		/* Allocate branch and node structure arrays */
	ckbrch = BranchVectorSpace(1, Nb);
	cknod = NodeVectorSpace(1, Nn);
	cps = VectorSpace(1, Nn);
	lps = VectorSpace(1, Nb);
	Nadi = DVectorSpace(1, Nn);
	Nado = DVectorSpace(1, Nn);
	ingo = DMatrixSpace(1, Nn, 1, adjmax);
	outgo = DMatrixSpace(1, Nn, 1, adjmax);

	/* Assign node numbers */
	ps1 = 2;
	ps2 = 3;
	ps3 = 4;
	ps4 = 5;
	gfict = 0.000000117939;
	cfict = 0.0001 * nano;
	rshfict = 1.0 / gfict;

	printf("READING BRANCHES... \n");
	/* Read in branch information */
	/* Data is arranged as follows:
	Branch ID - origin node - destination node - resistance (ohms) - inductance (nH) - Type
*/

	zer = 0;
	for (i = 1; i <= Nn; i++) {
		Nadi[i] = 0;
		Nado[i] = 0;
	}
	ndinter = Nn;

	inf1 = fopen("lout_jsa", "r");
	for (i = 1; i <= Nb; i++) {

		fscanf(inf1, "%ld %ld %ld %lf %lf %ld \n", &brid, &brn1, &brn2, &brres, &brlind, &brtyp);
		brtyp = 0;


		ckbrch[i].id = brid;
		ckbrch[i].n1 = brn1;
		ckbrch[i].n2 = brn2;
		ckbrch[i].res = brres;
		ckbrch[i].lind = brlind;
		ckbrch[i].typ = brtyp;
		ckbrch[i].cur = 0.0;
		ckbrch[i].is = 1.0e-15;

		Nadi[brn2]++;
		Nado[brn1]++;
		outgo[brn1][Nado[brn1]] = brid;
		ingo[brn2][Nadi[brn2]] = brid;

	}

	/* Read in node information */
	/* Data is read as follows:
	Node ID - conductance (mhos) - capacitance (nF)  - Type */

	printf("READING NODES... \n");

	for (i = 1; i <= Nn; i++) {

		nodtyp = 0;
		fscanf(inf1, "%ld %lf %lf %ld \n", &nodid, &nodgond, &nodcap, &nodtyp);
		//	printf("%5d %.5e %.5e %5d \n", &nodid, &nodgond, &nodcap, &nodtyp);

		cknod[i].id = nodid;
		cknod[i].cap = nodcap;
		cknod[i].gond = nodgond;
		cknod[i].typ = nodtyp;
		cknod[i].is = 1.0e-15;
		cknod[i].xum = 0.0;

	}


	for (i = 1; i <= Nb; i++) {
		lps[i] = ckbrch[i].lind;
		ckbrch[i].gb = (ckbrch[i].lind / dt - ckbrch[i].res / 2.0) / (ckbrch[i].lind / dt + ckbrch[i].res / 2.0);
		ckbrch[i].tb = 1.0 / (ckbrch[i].lind / dt + ckbrch[i].res / 2.0);
	}
	cknod[idvs].gond = 0.02;

	for (i = 1; i <= Nn; i++) {
		cps[i] = cknod[i].cap;
		cknod[i].gn = (cknod[i].cap / dt - cknod[i].gond / 2.0) / (cknod[i].cap / dt + cknod[i].gond / 2.0);
		cknod[i].tn = 1.0 / (cknod[i].cap / dt + cknod[i].gond / 2.0);
	}

	// --------------------------------------------------------------------------------------------------------------------
	// CUDA: Host-to-Device (Initialization)
	// --------------------------------------------------------------------------------------------------------------------
	Nb_size = (Nb + 1) * sizeof(branch); // +1 because of 1-based indexing, need to allocate extra space
	Nn_size = (Nn + 1) * sizeof(node);

	// Malloc device memory
	checkCudaErrors(cudaMalloc((void**)&d_ckbrch, Nb_size));
	checkCudaErrors(cudaMalloc((void**)&d_cknod, Nn_size));

	// Copy from Host to Device
	// Skip d_ckbrch[0] to maintain 1-based indexing
	checkCudaErrors(cudaMemcpy(d_ckbrch + 1, ckbrch + 1,
		Nb * sizeof(branch), cudaMemcpyHostToDevice));

	checkCudaErrors(cudaMemcpy(d_cknod + 1, cknod + 1,
		Nn * sizeof(node), cudaMemcpyHostToDevice));

	// CUDA Stream and Event Initialization
	cudaStream_t stream = NULL;
	checkCudaErrors(cudaStreamCreate(&stream)); // Create a new CUDA stream

	cudaEvent_t node_complete_event = NULL;
	checkCudaErrors(cudaEventCreate(&node_complete_event)); // Create an event to mark the completion of the NODE kernel

	// -------------------------------------------------- CUDA: Host-to-Device End ------------------------------------------

	timex = VectorSpace(1, mt);
	vy1 = VectorSpace(1, mt);
	vy2 = VectorSpace(1, mt);
	vy3 = VectorSpace(1, mt);
	vy4 = VectorSpace(1, mt);
	printf("STARTING SIMULATION \n");

	tstart = 0;
	tstop = mt;
	kount = 0;
	tim = 0.0;
	t = -1;
	printf("TIME STEP IS: %.5e \n", dt);
	
	// ====== [Added] higher-precision total and phase timing variables (no change to original logic) ======
	using Clock = std::chrono::steady_clock;                        // [Added] use a steady clock for timing
	std::chrono::duration<double> acc_branch(0.0), acc_node(0.0);       // [Added] accumulators: BRANCH and NODE phase times


	// ============================ Version Choice Loop ================================
	while (true) {
		int version = -1;
		std::cout << "Select which version of execution( 1:Origin; 2:CUDA; 3:Quit;): ";
		std::cin >> version;
		if (version != 1 && version != 2 && version != 3) {
			std::cout << "Value invalid, please enter again. ";
			continue;
		}
		if (version == 3) break;

		// ===== [Added] Reset timer accumulators =====
		acc_branch = std::chrono::duration<double>::zero();
		acc_node = std::chrono::duration<double>::zero();

		// ===== [Added] Reset simulation =====
		t = -1;
		tim = 0.0;
		kount = 0;

		// [Added] Clear node voltages and branch currents
		for (int i = 1; i <= Nn; ++i) {
			cknod[i].xum = 0.0;
			cknod[i].vol = 0.0;
		}
		for (int i = 1; i <= Nb; ++i) {
			ckbrch[i].cur = 0.0;
		}

		// CUDA: Copy the reset host state to the device
		if (version == 2) {
			checkCudaErrors(cudaMemcpy(d_ckbrch + 1, ckbrch + 1,
				Nb * sizeof(branch),
				cudaMemcpyHostToDevice));

			checkCudaErrors(cudaMemcpy(d_cknod + 1, cknod + 1,
				Nn * sizeof(node),
				cudaMemcpyHostToDevice));
		}

		time(&td_begin); // original time function call
		auto t_all_begin = Clock::now();                                // [Added] mark simulation start time

		//---------------Staring time loop -----------------------------------
		while (tim <= dur) {
			t++;
			tim = dt * t;
			timex[t] = tim;

			/*-------------------------------------------------------*/
			/*                   BRANCH LOOP*/
			/*-------------------------------------------------------*/
			auto t_beg = Clock::now();  // [Added] mark start of BRANCH loop for this step
			if (version == 1){
				for (i = 1; i <= Nb; i++) {
					brid = ckbrch[i].id;
					i12 = ckbrch[i].cur;
					l12 = ckbrch[i].lind;
					r12 = ckbrch[i].res;
					nx1 = ckbrch[i].n1;
					nx2 = ckbrch[i].n2;
					v1 = cknod[nx1].vol;
					v2 = cknod[nx2].vol;
					brch = ckbrch[i];
					esr = 0.0;
					Gb = ckbrch[i].gb;
					Tb = ckbrch[i].tb;
					Gni = cknod[nx1].gn;
					Tni = cknod[nx1].tn;
					Gnj = cknod[nx2].gn;
					Tnj = cknod[nx2].tn;
					cknod[ckbrch[i].n1].xum -= ckbrch[i].cur;
					cknod[ckbrch[i].n2].xum += ckbrch[i].cur;
					ixum1 = cknod[ckbrch[i].n1].xum;
					ixum2 = cknod[ckbrch[i].n2].xum;

					if (imethod == 1)
						ckbrch[i].cur = i12 + dt / l12 * (v1 - v2 - r12 * i12);  //Old Method
					else
						ckbrch[i].cur = (Tb * Gni * v1 - Tb * Tni * ixum1 - Tb * Gnj * v2 + Tb * Tnj * ixum2 + Gb * ckbrch[i].cur) / (1.0 + Tb * Tni + Tb * Tnj); // New Method
					cknod[ckbrch[i].n1].xum += ckbrch[i].cur;
					cknod[ckbrch[i].n2].xum -= ckbrch[i].cur;
				}
			} else { 
				// CUDA: BRANCH Kernel Launch
				int threadsPerBlock = 64;
				int blocksPerGrid = (Nb + threadsPerBlock - 1) / threadsPerBlock;

				branch_kernel << < blocksPerGrid, threadsPerBlock, 0, stream >> > (d_ckbrch, d_cknod, Nb, dt, imethod);
				checkCudaErrors(cudaPeekAtLastError());
			}
			acc_branch += Clock::now() - t_beg;  // [Added] accumulate BRANCH loop time for this step

			/*-------------------------------------------------------*/
			/*		                NODE LOOP*/
			/*-------------------------------------------------------*/
			auto t_n_beg = Clock::now();  // [Added] mark start of NODE loop for this step
			if (version == 1) {
				for (i = 1; i <= Nn; i++) {
					g1 = cknod[i].gond;
					c1 = cknod[i].cap;
					v1 = cknod[i].vol;
					Gni = cknod[i].gn;
					Tni = cknod[i].tn;
					nod = cknod[i];
					hsr = 0.0;
					isum = 0.0;

					isum = cknod[i].xum;

					if (cknod[i].id == idvs) {
						comp(tim, del, ris, fal, wid, &zet);
						zeti = magi * zet;
						hs1 = -zeti;
						isum -= zeti;
					}
					else hs1 = 0.0;
					si = isum;
					if (imethod == 1)
						vnew = (v1 * c1 / dt - isum) / (c1 / dt + g1); //Old Method
					else
						vnew = Gni * v1 - Tni * isum; // New Method

					cknod[i].vol = vnew;
				}
			} else {
				// CUDA: NODE Kernel Launch
				int threadsPerBlock = 64;
				int blocksPerGrid = (Nn + threadsPerBlock - 1) / threadsPerBlock;

				// All params (Including: tim, del, ris, fal, wid, magi)
				node_kernel << < blocksPerGrid, threadsPerBlock, 0, stream >> > (d_cknod, Nn, tim, dt, imethod, idvs, magi, del, ris, fal, wid);
				checkCudaErrors(cudaPeekAtLastError());

				// Memcpy monitored nodes back to Host
				// Enqueue four independent memory copy operations in the same stream, which will start after all preceding tasks complete
				// Since ps1, ps2, ps3, and ps4 are contiguous constants, these copies can be issued in parallel in a single operation
				// Asynchronous operations allow the CPU to continue execution, improving overall efficiency
				cudaMemcpyAsync(&h_monitored_nodes[1], &d_cknod[ps1], 4 * sizeof(node), cudaMemcpyDeviceToHost, stream);

				// cudaEventRecord acts like stamping a marker in the stream, indicating that all preceding work (node_kernel) must complete before this point is considered reached
				checkCudaErrors(cudaEventRecord(node_complete_event, stream));
			} 

			/*---------------------------------------------- END OF NODE LOOP -----------------------------------------------------*/
			
			// CUDA: Before accessing the monitored node data on the host, ensure that the asynchronous copy has completed
			if (version == 2) {
				// Block the current thread until all tasks in the GPU stream launched before (and including) the node_complete_event have completed
				checkCudaErrors(cudaEventSynchronize(node_complete_event));
			}

			acc_node += Clock::now() - t_n_beg;  // [Added] accumulate NODE loop time for this step

			// CUDA: Choosing monitored nodes voltage values depending on version (h_monitored_nodes = local cache)
			double v1_val = (version == 2 ? h_monitored_nodes[1].vol : cknod[ps1].vol);
			double v2_val = (version == 2 ? h_monitored_nodes[2].vol : cknod[ps2].vol);
			double v3_val = (version == 2 ? h_monitored_nodes[3].vol : cknod[ps3].vol);
			double v4_val = (version == 2 ? h_monitored_nodes[4].vol : cknod[ps4].vol);

			if (kount == 0)
				printf("%10.5f \t %10.5f \t %10.5f \t %10.5f \t %10.5f \n", tim, v1_val, v2_val, v3_val, v4_val); 

			vy1[t] = v1_val;
			vy2[t] = v2_val;
			vy3[t] = v3_val;
			vy4[t] = v4_val;
			kount++;

			if (kount == kmaxx)
				kount = 0;
		}

		// CUDA: Device-to-Host (Finalization)
		if (version == 2) {
			checkCudaErrors(cudaMemcpy(ckbrch + 1, d_ckbrch + 1,
									   Nb * sizeof(branch),
									   cudaMemcpyDeviceToHost));

			checkCudaErrors(cudaMemcpy(cknod + 1, d_cknod + 1,
									   Nn * sizeof(node),
									   cudaMemcpyDeviceToHost));
		}

		time(&td_ended);
		auto t_all_end = Clock::now();                                              // [Added] mark simulation end time
		td_du = difftime(td_ended, td_begin);
		printf("DONE with Simulations...%.5e seconds with %5d points \n", td_du, t);

		// ====== [Added] print higher-precision total time and phase breakdown (no change to original outputs) ======
		double Ttot = std::chrono::duration<double>(t_all_end - t_all_begin).count(); // [Added] total time (s)
		double Tbr = acc_branch.count();                                           // [Added] accumulated BRANCH time (s)
		double Tnd = acc_node.count();                                             // [Added] accumulated NODE time (s)
		double Toth = Ttot - Tbr - Tnd;                                             // [Added] other time (init, finalization, plotting prep, etc.)
		printf("\n==== Timing (Release x64) ====\n");                               // [Added] header
		printf("Total:   %.3f s\n", Ttot);                                          // [Added] total
		printf("Branch:  %.3f s (%.1f%%)\n", Tbr, 100.0 * Tbr / Ttot);              // [Added] BRANCH breakdown
		printf("Node:    %.3f s (%.1f%%)\n", Tnd, 100.0 * Tnd / Ttot);              // [Added] NODE breakdown
		printf("Other:   %.3f s (%.1f%%)\n\n", Toth, 100.0 * Toth / Ttot);          // [Added] others
		// =====================================================================

		char lgd1[62] = "V1"; char lgd2[62] = "V2"; char lgd3[62] = "V3"; char lgd4[62] = "V4";
		char Mttle[124] = "PDN"; char Yttle[24] = "Voltage (V)"; char Xttle[24] = "Time (ns)";

		thexyplot(mt, timex, vy1, vy2, vy3, vy4, Mttle, Xttle, Yttle, lgd1, lgd2, lgd3, lgd4);
	}
	// Release CUDA resources
	checkCudaErrors(cudaStreamDestroy(stream));
	checkCudaErrors(cudaEventDestroy(node_complete_event));

	// CUDA: Free Device Memory
	checkCudaErrors(cudaFree(d_ckbrch));
	checkCudaErrors(cudaFree(d_cknod));

	checkCudaErrors(cudaFreeHost(h_monitored_nodes));

	return 0;
}

/*--------------------------------------------------------------*/
/*        *DVectorSpace                                       */
/*--------------------------------------------------------------*/
int* DVectorSpace(int nl, int nh)
{
	int* v;

	v = (int*)calloc(((unsigned)(nh - nl + 1)), sizeof(int));
	return v - nl;
}
/*--------------------------------------------------------------*/
/*       FreeDVector                                          */
/*--------------------------------------------------------------*/
void  FreeDVector(int* v, int nl, int nh)
{
	free((char*)(v + nl));
	return;
}
/*--------------------------------------------------------------*/
/*        *BranchVectorSpace                                */
/*--------------------------------------------------------------*/
branch* BranchVectorSpace(int nl, int nh)
{
	branch* v;

	v = (branch*)calloc(((unsigned)(nh - nl + 1)), sizeof(branch));
	return v - nl;
}
/*--------------------------------------------------------------*/
/*        *NodeVectorSpace                                */
/*--------------------------------------------------------------*/
node* NodeVectorSpace(int nl, int nh)
{
	node* v;

	v = (node*)calloc(((unsigned)(nh - nl + 1)), sizeof(node));
	return v - nl;
}
/*--------------------------------------------------------------*/
/*       FreeBranchVector                                           */
/*--------------------------------------------------------------*/
void  FreeBranchVector(branch* v, int nl, int nh)
{
	free((char*)(v + nl));
	return;
}
/*--------------------------------------------------------------*/
/*       FreeNodeVector                                           */
/*--------------------------------------------------------------*/
void  FreeNodeVector(node* v, int nl, int nh)
{
	free((char*)(v + nl));
	return;
}
/*-------------------------------------------------------------------*/
/*                              comp                                 */
/*-------------------------------------------------------------------*/
void comp(double t, double a, double ris, double fal, double wid, double* att)
{
	double b, c, d, p;
	b = a + ris;
	c = b + wid;
	d = c + fal;
	p = t;
	/*    function is zero before a*/
	if (p < a) *att = 0;
	else if (p < b)  *att = 1 / (b - a) * (p - a);
	/*    rising edge of pulse*/
	/*    top of pulse (flat portion)*/
	else if (p < c) *att = 1;
	/*    falling edge of pulse*/
	else if (p < d)  *att = 1 / (c - d) * (p - d);
	else  *att = 0;
	return;
}


/*--------------------------------------------------------------*/
/*        *VectorSpace                                        */
/*--------------------------------------------------------------*/
double* VectorSpace(int nl, int nh)
{
	double* v;

	v = (double*)calloc(((unsigned)(nh - nl + 1)), sizeof(double));
	return v - nl;
}
/*--------------------------------------------------------------*/
/*       FreeVector                                           */
/*--------------------------------------------------------------*/
void  FreeVector(double* v, int nl, int nh)
{
	free((char*)(v + nl));
	return;
}

/*--------------------------------------------------------------*/
/*        **DMatrixSpace                                       */
/*--------------------------------------------------------------*/
int** DMatrixSpace(int nrl, int nrh, int ncl, int nch)
{
	int i;
	int** a;

	a = (int**)calloc(((unsigned)(nrh - nrl + 1)), sizeof(int*));
	a -= nrl;
	for (i = nrl; i <= nrh; i++)
	{
		a[i] = (int*)calloc(((unsigned)(nch - ncl + 1)), sizeof(int));
		a[i] -= ncl;
	}
	return a;
}
/*--------------------------------------------------------------*/
/*       FreeDMatrix                                           */
/*--------------------------------------------------------------*/
void  FreeDMatrix(int** a, int nrl, int nrh, int ncl, int nch)
{
	int i;
	for (i = nrh; i >= nrl; i--)
		free((char*)(a[i] + ncl));
	free((char*)(a + nrl));
	return;
}


//-------------------------------------------------------------
//                       STOP
//---------------------------------------------------------------
void stop()
{
	int idumi;
	getchar();
	//	scanf("%ld", &idumi);
	return;
}
