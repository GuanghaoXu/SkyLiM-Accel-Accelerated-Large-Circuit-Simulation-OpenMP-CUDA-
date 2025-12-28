// skylim.cpp : This file contains the 'main' function. Program execution begins and ends there.
//
/*  Program to simulate large circuits
Started Tuesday February 9,1999, 5:40 PM
Copyright Jose Schutt-Aine
University of Illinois    */
// Modified to work on Visual Studio. 
// July 10, 2018

// Modified to use OpenMP parallel computing by student Guanghao Xu
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
#include <omp.h>    

#include "lg.h"
#include <windows.h>
int* Nadi, * Nado, ** ingo, ** outgo;
branch* ckbrch;

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

	// ====== OpenMP tuning settings ======
	int num_threads = omp_get_num_procs();                  // get number of available processor cores
	omp_set_num_threads(num_threads);                       // set number of threads
	printf("OMP threads = %d\n", omp_get_max_threads());    // print thread count for confirmation
	
	// ====== [Added] higher-precision total and phase timing variables (no change to original logic) ======
	using Clock = std::chrono::steady_clock;                        // [Added] use a steady clock for timing
	std::chrono::duration<double> acc_branch(0.0), acc_node(0.0);       // [Added] accumulators: BRANCH and NODE phase times


	// ============================ Version Choice Loop ================================
	while (true) {
		int version = -1;
		std::cout << "Select which version of execution( 1:Origin; 2:Optimized; 3:Quit): ";
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
			}
			else {
#pragma omp parallel for schedule(dynamic,128)
				for (int i = 1; i <= Nb; ++i) {
					// --- Read globals once, copy to locals ---
					const double i_old = ckbrch[i].cur;
					const double l12 = ckbrch[i].lind;
					const double r12 = ckbrch[i].res;
					const int    n1 = ckbrch[i].n1;
					const int    n2 = ckbrch[i].n2;
					const double v1 = cknod[n1].vol;
					const double v2 = cknod[n2].vol;
					const double Gb = ckbrch[i].gb;
					const double Tb = ckbrch[i].tb;
					const double Gni = cknod[n1].gn, Tni = cknod[n1].tn;
					const double Gnj = cknod[n2].gn, Tnj = cknod[n2].tn;

					// Equivalent "snapshot" ixum
					const double ixum1 = cknod[n1].xum - i_old;
					const double ixum2 = cknod[n2].xum + i_old;

					// compute i_new
					double i_new;
					if (imethod == 1)
						i_new = i_old + dt / l12 * (v1 - v2 - r12 * i_old);
					else
						i_new = (Tb * Gni * v1 - Tb * Tni * ixum1 - Tb * Gnj * v2 + Tb * Tnj * ixum2 + Gb * i_old)
						/ (1.0 + Tb * Tni + Tb * Tnj);

					ckbrch[i].cur = i_new;

					// Key: update node with net increment to preserve equivalence to sequential order
					cknod[n1].xum += (-i_old + i_new);
					cknod[n2].xum += (+i_old - i_new);
				}
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
			}
			else {
#pragma omp parallel for schedule(dynamic,256)
				for (int i = 1; i <= Nn; ++i) {
					// --- Read only once ---
					const double g1 = cknod[i].gond;
					const double c1 = cknod[i].cap;
					const double v_old = cknod[i].vol;
					const double Gni = cknod[i].gn;
					const double Tni = cknod[i].tn;
					const double xum = cknod[i].xum;
					const int    nid = cknod[i].id;

					double zeti = 0.0;
					double hs1 = 0.0;

					if (nid == idvs) {
						double zet;
						comp(tim, del, ris, fal, wid, &zet);
						zeti = magi * zet;
						hs1 = -zeti;
					}

					const double isum = xum - zeti;

					// calculate v_new
					double v_new;
					if (imethod == 1) {
						v_new = ((v_old * c1) / dt - isum) / (c1 / dt + g1);
					}
					else {
						v_new = Gni * v_old - Tni * isum;
					}
					cknod[i].vol = v_new;
				}
			}
			acc_node += Clock::now() - t_n_beg;  // [Added] accumulate NODE loop time for this step

			/*---------------------------------------------- END OF NODE LOOP -----------------------------------------------------*/

					/* printf("%10.5f %10.5f %10.5f \n",tim,cknod[ps1].vol,cknod[ps2].vol); */

			if (kount == 0)

				printf("%10.5f \t %10.5f \t %10.5f \t %10.5f \t %10.5f \n", tim,
					cknod[ps1].vol, cknod[ps2].vol, cknod[ps3].vol, cknod[ps4].vol);

			vy1[t] = cknod[ps1].vol;
			vy2[t] = cknod[ps2].vol;
			vy3[t] = cknod[ps3].vol;
			vy4[t] = cknod[ps4].vol;
			kount++;
			if (kount == kmaxx)
				kount = 0;
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
