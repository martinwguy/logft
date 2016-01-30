/* LOGFT.C
Takes as infile a soundfile of shorts (16 bit ints);
throws away the first 1024 bytes assuming they're a soundfile header.
See dumphdr() to get rid of this feature.


adaptation of dft.c using log channels separated by
/*       a semitone with resolution of a semitone
/*       later modified; two chnls/semitone with variable resolution
/*       calculates discrete Fourier transform of binary short sound samples
/*  (from stdin) and writes binary float coefficients (mag or db) to stdout.
/*       9/90 include option so in pitch tracking doesn't report the channels "between         freqs corresp to musical notes

/* usage: logft [options] <infile >outfile
/*   -hHANNING(1) Use Hanning or hamming window(1) or (0) rectangular.
/*   -HHAMMING(0) For Hamming instead of Hanning
/*   -dDB (0)     Calculate logft in decibels.
/*   -rRES(17)    Resolution f/deltaf.
/*   -FFRMSZ(500) Number of samples analyzed per frame.
/*   -pNCHNLS(156)Number of frequencies for which logft calculated per frame
/*   -fMINHZ(174.6)Lowest frequency value = 7.05 = F3 = midi 53
/*   -aHALFSEMI(1) Calculate with 2 frequency values per semitone.
/*   -TTUNCORREC   Follow with tuning correction for sounds digitized at
                  ems(1.04) implemented 7/20/88
/*   -WWINDSIZ-1(0)Use windsiz - 1 in hanning(mm) window
/*   -SSRATE       Sample rate (32K)
/*   -AWRITEHDR(0)  Write a header on the file if = 1.
/*   -dropmidchnls  Drop channels between those tuned to notes
/***************************************************************************/

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <fcntl.h> /* defines O_RDONLY */
#include "./frmhdr.h"

#define FRMSZ 500
#define WINDMAX 20000
/*#define MINHZ 174.6  /* F3 = 7.05  midi=53 */
#define MINHZ 130.80 /* c3 midi = 48 */
/*#define MAXHZ 1600  */
#define RES 34 /* f/deltaf = 1/.06 */
#define MAXCHNLS 200
#define NCHNLS 156
#define SRATE  32000.
#define TWOPI  6.2831854
#define PI  3.1415927
#define PMODE 0644
#define HANMAX 700000 /*ems can't handle 700000 /* works with 500000 */
#define MINLAM 4                       /* at 10k 6  */
#define MAXLAM 250  /* for c3 ; was 195                      /* at 10k 75 */

static struct  frmhdr hdr;
static short   dumphdr[1024];                  /* dump header */
static short   sampbuf[WINDMAX];               /* buffer for input samples        */
static float   outbuf[MAXCHNLS];     /* frame of out coefs in mag or db */
static float   srate = SRATE, midtuncor=1.;

main(argc,argv)
 int argc;
 char **argv;
{  /* begin main */
        double  *hanfilrd;
        double  theta, a, b, c, db, one06, one03;
        double  onedws, twopidws, RES2pidws, alpha,phase;
        int     frmsz = FRMSZ, windbytes, framebytes, nchnls=NCHNLS;
        int     hanning = 1, dbout = 0, frmcnt = 0;
        int     windsiz[MAXCHNLS];
        int     flag, windmaxi, sumwind=0, res = RES;
        int     nw;
        int     halfsemi=1;
        int     maxreq=0;
        int     hamming = 1, windsiz_1 = 0, calcphase=0;
        int     writehdr=0, dropmidchnls=0;
        char    *cp;
        float   windmax, minhz= MINHZ, tuncorrec=1.00;
        float   windsizf[MAXCHNLS];
        register int    n, k;
        register short  *samp, *samp2;
        register double *hanp, *hanfilrdp;

        argc--; argv++;                         /* align onto command args */
        while ((cp = *argv) && *cp++ =='-' && (flag = *cp++))
        {
                switch (flag){
                        case 'h': sscanf(cp, "%d", &hanning);
                                break;
                        case 'H': sscanf(cp, "%d", &hamming);
                                break;
                        case 'd': sscanf(cp, "%d", &dbout);
                                break;
                        case 'f': sscanf(cp, "%f", &minhz);
                                break;
                        case 'F': sscanf(cp, "%d", &frmsz);
                                break;
                        case 'r': sscanf(cp, "%d", &res);/* resolution f/delf*/
                                break;
                        case 'p':sscanf(cp,"%d",&nchnls);/*num points in ft*/
                                break;
                        case 'a': sscanf(cp, "%d", &halfsemi);
                                break;
                        case 'W':  sscanf(cp, "%d", &windsiz_1);
                                break;
                        case 'M': sscanf(cp, "%d", &maxreq);
                                break;
                        case 'T': sscanf(cp, "%f", &tuncorrec);
                                break;
                        case 'S': sscanf(cp, "%f", &srate);
                                break;
                        case 's': sscanf(cp, "%d", &calcphase);
                                break;
                        case 'Q': sscanf(cp, "%f", &midtuncor);
                                break;
                        case 'A': sscanf(cp, "%d", &writehdr);
                                break;
                        case 'D': sscanf(cp, "%d", &dropmidchnls);
                                break;

                        default: fprintf(stderr, "unknown option\n");
                                break;
                }
                argc--; argv++;
        }

        hanfilrd = (double *)malloc(HANMAX*8);
        if (hanfilrd == NULL)
          die ("memory allocation failure");

        if(hamming) {
                 alpha=(double)25./(double)46.;
                 fprintf(stderr,"Using Hamming window\n");
        }
        else alpha = (double).5; /* hanning */
        if(nchnls > MAXCHNLS)die("Too many channels");

        minhz /= tuncorrec;

        if (writehdr){
           fprintf(stderr, "Writing a header on outfile\n");
           hdr.frm_magic = FRM_MAGIC;
           hdr.frm_width = nchnls;
           hdr.frm_datatype = FRM_DATA_FLOAT;
           hdr.frm_scale = FRM_SCALE_LOG;
           hdr.frm_frm_size = frmsz;
           hdr.frm_sr = srate;
           hdr.frm_frm_p_s = srate/frmsz;
           hdr.frm_minfreq = minhz;
           if((nw = write(1, &hdr, sizeof hdr)) != sizeof hdr)
	     die("Error writing file header");
        }
        one06 = pow(2.,1./12.);
        one03 = pow(2.,1./24.); /* for stepping by 1/2 semitone */
                                     /* freqrat = (float)MAXHZ /minhz; */
                                /*nchnls = (int)(log(freqrat)/log(one06)) + 1*/
        windmax = (float)(res*srate)/minhz; /* res determines harmonic
                    number which is const=res; freq is varied with windsiz[k]
                    and equals  srate*res/windsiz[k] */
        if (n= read(0, dumphdr, 1024) != 1024)
          die("something wrong with header dump\n");

        k=nchnls;       n=0;
        while (k--){
           if(halfsemi==1){  /* 2 chnls per semitone */
             windsizf[n] = (float)windmax/pow(one03,(float)n);
           }
           else{
             windsizf[n] = (float)windmax/pow(one06,(float)n);
           }
           windsiz[n] = (int)windsizf[n];
           sumwind += windsiz[n];        /* total window space needed
                                         for sin tables = sum of all windows */
           ++n;
        }
	if(2*sumwind > HANMAX)die("2*Sumwind is > HANMAX.");
	windmaxi = (int)windmax;     /* no of samples to read in is the size
					of the largest window (lowest freq*/
	if(windmaxi > WINDMAX)die("Windmaxi larger than WINDMAX");
	windbytes = windmaxi *2;     /*windmaxi is old variable windsiz in dft*/
	framebytes = frmsz * 2;

	/* Calculate window (or rect) values */
	hanp = hanfilrd;
	for (k=0; k < nchnls; ++k){
	    if(windsiz_1)twopidws = TWOPI/(windsiz[k] - 1);
	    else twopidws = TWOPI/(windsiz[k]);
	    RES2pidws = res* TWOPI/windsiz[k];
	    onedws = 1./windsiz[k];

	    for (n=0; n < windsiz[k]; ++n){
	       a=onedws*((!hanning)?1.:alpha-((1-alpha)*cos(n*twopidws)));
	       theta = n * RES2pidws;
	       *hanp++ = a * sin(theta);
	       *hanp++ = a * cos(theta);
	    }
	}

        n = read(0, sampbuf, windbytes);      /* init samp window w. input */
        if (n != windbytes)
              die("premature end of infile");
frame:

	for (k=0,hanfilrdp=hanfilrd; k<nchnls; k++) { /* for one frame: */
	   a = 0.0;
	   b = 0.0;
	   samp =sampbuf;
	   for (n=0; n<windsiz[k]; n++) { /*  calculate coefs  */
	      a += *samp * *hanfilrdp++;
	      b += *samp++ * *hanfilrdp++; /*wrote sin then cos*/
	   }
	   if(calcphase){
	      if ( a < .01) phase = 0. ;
	      else {
	  	phase = atan(a/b);
	  	outbuf[k] = phase;
	      }
	   } else {
	      c = sqrt( a*a + b*b );
	      if (!dbout)
	 	outbuf[k] = c;       /* stor as magnitude */
	      else outbuf[k] = db = 20. * log10(c); /*       or db        */
	   }
	}                            /* end calc from table read in */

       if((nw = write(1, outbuf, 4*nchnls)) < 0)/*print outbuf to file*/
          fprintf(stderr, "nw = %d  ERROR\n", nw);

       if (nw < 4*nchnls){
            fprintf(stderr, "error writing from outbuf\n wrote %d bytes\
                and should have written %d\n", nw, 4*nchnls);
       }

       ++frmcnt;

       if ((n = windmaxi - frmsz) > 0) {       /* nxt: for step < windowsiz */
                samp = sampbuf;
                samp2 = sampbuf + frmsz;
                while (n--)
                        *samp++ = *samp2++;     /*      slide the sample buf */
                n = read(0, samp, framebytes);  /*     & refill to windmaxi  */
                if (n == framebytes){
                        goto frame;
                }
        } else {
                for (n = -n; n>=windmaxi; n-=windmaxi)  /* else waste any   */
                        read(0, sampbuf, windbytes);    /*   extra samps    */
                if (n)  read(0, sampbuf, n*2);
                n = read(0, sampbuf, windbytes); /* then rd whole new window */
                if (n == windbytes){

                        goto frame;             /*      & go calc new frame */
                }
        }
} /* end main */

die(s)
 char *s;
{
        fprintf(stderr,"DIE, HELAS!  %s\n",s);
        exit(1);
}
