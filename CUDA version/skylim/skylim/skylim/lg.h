/*  Program to generate mesh  
Started Tuesday January 19,1999, 11:34 AM
   Copyright Jose Schutt-Aine
   University of Illinois    */
   
#include <stdio.h>
#include <math.h>
#include <stdlib.h>

typedef struct branch {
                        int id; /* ID number of branch */
                        int n1; /* Node from which current is assumed to originate */
                        int n2; /* Node toward which current is directed */
                        int typ; /* type of branch 0: linear, 1: forward-biased
                                  nonlinear element, -1: reversed biased 
                                  nonlinear element. All is when current is flowing
                                  from n1 to n2 */
                        double res; /* resistance of branch */
                        double lind;/* inductance of branch */
                        double cur; /* value of current in branch (from n1 to n2) */
                        double curp;/* current value from previous time step */
                        double curpp;/* current value from 2 time steps ago */
                        double cs; /* not used */
                        double vtrap; /* not used */
                        double is; /* saturation current for diode
                                      nonlinear element */
						double gb;
						double tb;
                      } branch;
                      
typedef struct node {
                        int id; /* ID number of node */
                        int typ; /* type ofnode 0: linear, 1: forward-biased
                                  nonlinear element, -1: reversed biased 
                                  nonlinear element. All is when current is flowing
                                  to ground */
                        double gond; /* conductance from node to ground */ 
                        double cap; /* capacitance from node to ground */ 
                        double vol; /* node voltage */ 
                        double volp;/* node voltage from previous time step */
                        double volpp;/* node voltage from 2 time steps ago */
                        double vs; /* not used */
                        double itrap; /* not used */
                        double is; /* saturation current for diode
                                      nonlinear element */
						double xum;
						double gn;
						double tn;
                      } node;
        
extern branch *BranchVectorSpace(int nl, int nh);
extern node *NodeVectorSpace(int nl, int nh);
extern void  FreeBranchVector(branch *v, int nl, int nh);
extern void  FreeNodeVector(node *v, int nl, int nh);
extern int *DVectorSpace(int nl, int nh);
extern void  FreeDVector(int *v, int nl, int nh);
extern void comp(double t,double a,double ris,double fal,
double wid,double *att);
extern double getslc(double *lps,double *cps, int Nb, int Nn);
extern void  FreeVector(double *v, int nl, int nh);
extern double *VectorSpace(int nl, int nh);
extern void setbranch(int ident,int i1,int i2,branch *br);
extern void setnode(int i,node *ckn);
extern double genrand(int seed,double vmin,double vmax);
extern int isstable(node ckn,int i,double tim,double dt);
extern int **DMatrixSpace(int nrl, int nrh, int ncl, int nch);
extern void  FreeDMatrix(int **a, int nrl, int nrh, int ncl, int nch);
extern double citerate(branch brnl,double dt,double vi,double vj,double esr);
extern double viterate(node nd,double dt,double si,double hsr);

