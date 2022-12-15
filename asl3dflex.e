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

/* set preprocessor globals (#define functions) */

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

/* initialize pulse sequence cvs */

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
	fprintf(stderr, "\nStarting cvinit");

	/* configure system and optimize gradients */
	configSystem();
	EpicConf();
	inittargets(&loggrd, &phygrd);
	fudgetargets(&loggrd, &phygrd, rtimescale);
	
	/* optimize targets */
	if (obloptimize(&loggrd, &phygrd, scan_info, exist(opslquant), exist(opplane),
		exist(opcoax), obl_method, obl_debug, &opnewgeo, cfsrmode)==FAILURE) {
		return FAILURE;
	}

	/* set psd_rf_wait and psd_grd_wait */
	if (_psd_rf_wait.fixedflag == 0)  { 
                if (setsysparms() == FAILURE)  { 
                        epic_error(use_ermes,"Support routine setsysparams failed", 
                                        EM_PSD_SUPPORT_FAILURE,1, STRING_ARG,"setsysparms"); 
                        return FAILURE; 
                } 
        }	

	/* inline the prescan cvinit section */
@inline Prescan.e PScvinit
#include "cvinit.in"

	/* initialize high-level op cvs */

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

	/* configure system and initiate advanced cvs panel */
	configSystem();
	InitAdvPnlCVs();

	/* inline the prescan cveval section */
@inline Prescan.e PScveval

	return SUCCESS;
}
