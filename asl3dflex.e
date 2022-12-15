@inline epic.h

@global
/*********************************************************************
 *                        GLOBAL SECTION                             *
 *                                                                   *
 * Common code shared between the Host and IPG PSD processes.  This  *
 * section contains all the #define's, global variables and function *
 * declarations (prototypes).                                        *
 *********************************************************************/

/* include all basic libraries */
#include <stdio.h>
#include <string.h>

#include "em_psd_ermes.in"
#include "stddef_ep.h"
#include "epicconf.h"
#include "pulsegen.h"
#include "epic_error.h"
#include "epicfuns.h"
#include "epic_loadcvs.h"
#include "InitAdvisories.h"
#include "psdiopt.h"
#include "psdutil.h"
#include "psd_proto.h"
#include "epic_iopt_util.h"
#include "filter.h"

/* include grad rf globals header file */
#include "grad_rf_asl3dflex.globals.h"

/* inline intwave.h */
@inline intwave.h

/* set preprocessor globals (#define functions) */
#define MAXWAVELEN      50000 /* Maximum waveform array size */
#define MAXNUMECHOES    1000 /* Maximum number of echoes */
#define MAXITR          50 /* Maximum number of iterations in iterative processes */
#define TSP_GRAD        4 /* Scanner gradient sampling rate (us) */
#define TSP_RF          2 /* Scanner rf sampling rate (us) */
#define GOLDENANGLE     2.3999632297286531 /* Golden angle (rad) */
#define GAMMA           26754 /* Gyromagnetic ratio (rad/s/G) */

/* inline the prescan global section */
@inline Prescan.e PSglobal
int debugstate = 1; /* included from grass.e */


@ipgexport
/*********************************************************************
 *                      IPGEXPORT SECTION                            *
 *                                                                   *
 * Standard C variables of _any_ type common for both the Host and   *
 * IPG PSD processes. Declare here all the complex type, e.g.,       *
 * structures, arrays, files, etc.                                   *
 *                                                                   *
 * NOTE FOR Lx:                                                      *
 * Since the architectures between the Host and the IPG sides are    *
 * different, the memory alignment for certain types varies. Hence,  *
 * the following types are "forbidden": short, char, and double.     *
 *********************************************************************/

/* inline the prescan ipgexport section */
@inline Prescan.e PSipgexport

/* declare the rf pulse info object */
RF_PULSE_INFO rfpulseinfo[RF_FREE];

/* initialize non-cv global variables accessible by host and ipg processes*/


@cv
/*********************************************************************
 *                           CV SECTION                              *
 *                                                                   *
 * Standard C variables of _limited_ types common for both the Host  *
 * and IPG PSD processes. Declare here all the simple types, e.g,    *
 * int, float, and C structures containing the min and max values,   *
 * and ID description, etc.                                          *
 *                                                                   *
 * NOTE FOR Lx:                                                      *
 * Since the architectures between the Host and the IPG sides are    *
 * different, the memory alignment for certain types varies. Hence,  *
 * the following types are "forbidden": short, char, and double.     *
 *********************************************************************/

/* inline the system cvs */
@inline loadrheader.e rheadercv
@inline vmx.e SysCVs

/* initialize back-end cvs */
int obl_debug = 0 with {0, 1, 0, INVIS, "On(=1) to print messages for obloptimize",};
int obl_method = 0 with {0, 1, 0, INVIS, "On(=1) to optimize the targets based on actual rotation matrices",};

/* initialize readout cvs */
int nleaves = 1 with {1, 128, 1, VIS, "Number of interleaves in spiral readout",};
int nframes = 2 with {1, , 2, VIS, "Number of temporal frames",};
int nM0frames = 2 with {0, , 2, VIS, "Number of M0 frames (no ASL prep)",};
int nextra = 2 with {0, , 2, VIS, "Number of extra shots before acquisition",};
int do_shortrf = 1 with {0, 1, 1, VIS, "Use short RF pulses (1) or not (0)",};
float r_accel = 0.7 with {0.1, , 0.7, VIS, "Spiral radial acceleration factor",};
float theta_accel = 0.7 with {0.1, , 0.7, VIS, "Spiral angular acceleration factor",};

/* initialize asl cvs */
int do_bkgsup = 0 with {0, 1, 0, VIS, "Do background suppression (1) or not (0)",};
int do_artsup = 0 with {0, 1, 0, VIS, "Do arterial suppression (1) or not (0)",};

/* initialize other cvs */
int reconscript = 0 with {0, , 0, VIS, "Recon script number (0 = none)",};

/* inline the prescan cvs */
@inline Prescan.e PScvs


@host
/*********************************************************************
 *                          HOST SECTION                             *
 *                                                                   *
 * Write here the code unique to the Host PSD process. The following *
 * functions must be declared here: cvinit(), cveval(), cvcheck(),   *
 * and predownload().                                                *
 *                                                                   *
 *********************************************************************/

/* include basic libraries needed for host functions */
#include "support_func.host.h"  /* new for 28x */
#include <float.h>
#include <math.h>
#include <stdlib.h>  
#include "epic_iopt_util.h"
#include "psd.h"
#include "psdIF.h"
#include "psdopt.h"
#include "psdutil.h"
#include "psd_receive_chain.h"
#include "rfsspsummary.h"
#include "sar_burst_api.h"
#include "sar_display_api.h"
#include "sar_limit_api.h"
#include "sar_pm.h"
#include "support_func.host.h"

/* include field strength dependency library */
#include "sysDep.h"
#include "sysDepSupport.h"

/* include grad rf header file */
#include "grad_rf_asl3dflex.h"

/* inline the prescan and r header host sections */
@inline Prescan.e PShostVars
@inline loadrheader.e rheaderhost

/* load the psd header */
abstract("ASL sequence with flexible 3D FSE spiral readout");
psdname("asl3dflex");

/* declare echo and auto prescan filters */
FILTER_INFO echo1_filt;
FILTER_INFO aps2_filt;

/* include fudgetargets and make routine failure message */
extern "C" {
#include "fudgetargets.c"
}
static char supfailfmt[] = "Support routine %s exploded!";


STATUS cvinit( void )
{
/************************************************************************/
/*                              CVINIT                                  */
/* Invoked once (& only once) when the PSD host process is started up.  */
/* Code which is independent of any OPIO button operation is put here.  */
/************************************************************************/

	/* print out update */
	fprintf(stderr, "\nStarting cvinit...");

	/* configure system and optimize gradients */
	configSystem();
	EpicConf();
	inittargets(&loggrd, &phygrd);
	fudgetargets(&loggrd, &phygrd, rtimescale);	

	/* inline the prescan cvinit section */
@inline Prescan.e PScvinit
#include "cvinit.in"

	/* initialize flip andgle and turn off button */
	opflip = 90.0
	cvdef(opflip, 90.0);
	pifanub = 0;

	/* fix reciever bw and turn off button */
	oprbw = 500.0 / (float)TSP_GRAD;
	pircbnub = 0;

	/* initialize fov and its buttons */
	opfov = 240; /* mm */
	pifovnub = 5;
	pifovval2 = 200;
	pifovval3 = 220;
	pifovval4 = 240;
	pifovval5 = 260;
	pifovval6 = 280;

	/* initialize tr and its buttons */
	optr = 4500ms;
	pitrnub = 3;
	pitrval2 = PSD_MINIMUMTR;
	pitrval3 = 4500ms;
	pitrval4 = 5000ms;

	/* initialize tr and its buttons */
	opte = PSD_MINFULLTE;
	pitrnub = 2;
	pite1val2 = PSD_MINFULLTE;
	piteval3 = 100ms;

	/* initialize frequency (xres) and turn off button */
	cvmin(opxres, 32);
	cvmax(opxres, 512);
	cvdef(opxres, 64);
	opxres = 64;
	pixresnub = 0;

	/* hide phase (yres) option */
	piyresnub = 0;
	
	/* hide second bandwidth option */
	pircb2nub = 0;

	/* hide nex stuff */
	piechnub = 0;
	pinexnub = 0;

	return SUCCESS;
}

/* inline the InitAdvPnlCVs code */
@inline InitAdvisories.e InitAdvPnlCVs


STATUS cveval( void )
{
/************************************************************************/
/*                              CVEVAL                                  */
/* Called w/ every OPIO button push which has a corresponding CV.       */
/* CVEVAL should only contain code which impacts the advisory panel--   */
/* put other code in cvinit or predownload                              */
/************************************************************************/

	/* print out update */
	fprintf(stderr,"\nStarting cveval...");

	/* configure system and initiate advanced cvs panel */
	configSystem();
	InitAdvPnlCVs();
	
	/* set psd_rf_wait and psd_grd_wait */
	if (_psd_rf_wait.fixedflag == 0)  { 
                if (setsysparms() == FAILURE)  { 
                        epic_error(use_ermes,"Support routine setsysparams failed", 
                                        EM_PSD_SUPPORT_FAILURE,1, STRING_ARG,"setsysparms"); 
                        return FAILURE; 
                } 
        }	

	/* optimize targets */
	if (obloptimize(&loggrd, &phygrd, scan_info, exist(opslquant), exist(opplane),
		exist(opcoax), obl_method, obl_debug, &opnewgeo, cfsrmode)==FAILURE) {
		return FAILURE;
	}

	/* set up opuser (advanced cvs) tab */
	pititle = 1;
	cvdesc(pititle, "asl3dflex User CV Page");

	/* opuser0 = nleaves */
	cvdesc(opuser0, "Number of interleaves in spiral readout");
	cvdef(opuser0, 1);
	opuser0 = 1;
	cvmin(opuser0, 1);
	cvmax(opuser0, 128);
	nleaves = (int)opuser0;
	piuset = use0;

	/* opuser1 = nframes */
	cvdesc(opuser1, "Number of temporal frames");
	cvdef(opuser1, 2);
	opuser1 = 2;
	cvmin(opuser1, 1);
	nframes = (int)opuser1;
	piuset += use1;

	/* opuser2 = nM0frames */
	cvdesc(opuser2, "Number of M0 frames (no ASL prep)");
	cvdef(opuser2, 2);
	opuser1 = 2;
	cvmin(opuser1, 0);
	nM0frames = (int)opuser2;
	piuset += use2;

	/* opuser3 = recon script number */
	cvdesc(opuser3, "Recon script number (0 = none)");
	cvdef(opuser3, 0);
	opuser3 = 0;
	cvmin(opuser3, 0);
	reconscript = (int)opuser3;
	piuset += use3;

	/* opuser4 = number of extra shots */
	cvdesc(opuser4, "Number of extra shots before acquisition");
	cvdef(opuser4, 2);
	opuser4 = 2;
	cvmin(opuser4, 0);
	nextra = (int)opuser4;
	piuset += use4;
	
	/* opuser5 = use short rf pulses */
	cvdesc(opuser5, "Use short RF pulses (1) or not (0)");
	cvdef(opuser5, 1);
	opuser5 = 1;
	cvmin(opuser5, 0);
	cvmax(opuser5, 1);
	do_shortrf = (int)opuser5;
	piuset += use5;

	/* opuser6 = radial acceleration factor */
	cvdesc(opuser6, "Spiral radial acceleration factor");
	cvdef(opuser6, 0.7);
	opuser6 = 0.7;
	cvmin(opuser6, 0.1);
	r_accel = opuser6;
	piuser += use6;

	/* opuser7 = angular acceleration factor */
	cvdesc(opuser7, "Spiral angular acceleration factor");
	cvdef(opuser7, 0.7);
	opuser7 = 0.7;
	cvmin(opuser7, 0.1);
	theta_accel = opuser7;
	piuset += use7;

	/* opuser8 = do background suppression */
	cvdesc(opuser8, "Do background supression (1) or not (0)");
	cvdef(opuser8, 0);
	opuser8 = 0;
	cvmin(opuser8, 0);
	cvmax(opuser8, 1);
	do_bkgsup = (int)opuser8;
	piuset += use8;

	/* opuser9 = do arterial suppression */
	cvdesc(opuser9, "Do arterial suppression (1) or not (0)");
	cvdef(opuser9, 0);
	opuser9 = 0;
	cvmin(opuser9, 0);
	cvmax(opuser9, 1);
	do_artsup = (int)opuser9;
	piuset += use9;

	/* inline the prescan cveval section */
@inline Prescan.e PScveval

	return SUCCESS;
}


/* Declare APx functions for DV26 */
void getAPxParam(optval *min, optval *max, optdelta *delta,
	optfix *fix, float coverage, int algorithm)
{
}
int getAPxAlgorithm(optparam *optflag, int *algorithm)
{
	return APX_CORE_NONE;
}


STATUS cvcheck( void )
{
/************************************************************************/
/*                              CVCHECK                                 */
/* Executed on each 'next page' to ensure prescription can proceed      */
/* to the next page.                                                    */
/************************************************************************/

	return SUCCESS;
}


STATUS predownload( void )
{
/************************************************************************/
/*                          PRE-DOWNLOAD                                */
/* Executed prior to a download--all operations not needed for the      */
/* advisory panel results.  Execute the pulsegen macro expansions for   */
/* the predownload section here.  All internal amps, slice ordering,    */
/* prescan slice calc., and SAT placement calculations are performed    */
/* in this section.  Time anchor settings for pulsegen are done in this */
/* section too.                                                         */
/************************************************************************/

	return SUCCESS;
}
