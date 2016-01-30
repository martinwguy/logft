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

/* usage: logft -hHANNING -PPRINT -dDB -bDEBUG -rRES -FFRMSZ -pNCHNLS -fMINHZ
              -iHANFILE -jHANWRITTEN -lLINE -oNOOUT -aHALFSEMI -vVARRES
              -RFLTRESP -Ddropmidchnls
      -xXCORRI {-nNHARM -cCOMBPRINT-qBACKZERO -gGRAPH -zNRMXCORR -tXCORRPRNT
                -wLOGFTVALIN -kPIKPEAK {-mGRPHPEAK -ePRNTPIK -yMIDIFILE }}
       <infile >outfile
/*   -hHANNING(1) Use Hanning or hamming window(1) or (0) rectangular.
/*   -HHAMMING(0) For Hamming instead of Hanning
/*   -PPRINT(0)   Print values of logft to stderr.
/*   -dDB (0)     Calculate logft in decibels.
/*   -rRES(17)    Resolution f/deltaf.
/*   -FFRMSZ(500) Number of samples analyzed per frame.
/*   -pNCHNLS(156)Number of frequencies for which logft calculated per frame
/*   -fMINHZ(174.6)Lowest frequency value = 7.05 = F3 = midi 53
/*   -iHANFILE    Print array of sin & cos times hanning values to "hanfile".
/*   -jHANWRITTEN(0) Hanfile already written.
/*   -lLINE(0)       Just print the first frame.
/*   -oNOOUT(0)      Don't print an output file; use for writing number to a
                       file with >& .
/*   -aHALFSEMI(1) Calculate with 2 frequency values per semitone.
/*   -vVARRES(1)   Makes resolution 4 time original for midinotes > 90 (G6).
                  7/19/88 change so makes resolutin VARRES times orig w deflt 2
/*   -TTUNCORREC   Follow with tuning correction for sounds digitized at
                  ems(1.04) implemented 7/20/88
/*   -WWINDSIZ-1(0)Use windsiz - 1 in hanning(mm) window
/*   -SSRATE       Sample rate (32K)
/*   -AWRITEHDR(0)  Write a header on the file if = 1.
/*   -dropmidchnls  Drop channels between those tuned to notes
/*   -xXCORRI(0)   Perform cross correlation.
         -nNHARM(6)     Number of harmonics in xcorr fn.
         -cCOMBPRINT(0) Print cross correlation function.
         -qBACKZERO(1)  Cross corr fn = 0 between harmonics if backzero=1.
                          Otherwise the average of the x corr fn is 0.
         -gGRAPH(0)     Will print a frame of the x corr fn after each frame
                          of logft values.
         -zNRMXCORR(0)  Will normalize the frame of xcorr values so the max
                          of the xcorr fn = the max of the logft values.
         -tXCORRPRNT(0) Will print values of the cross correlation function.
         -wLOGFTVALIN   If logft values prev written to "logftvalin" will
                          read this so can try peak picker w/o recalculating.
         -MMAXREQ(0)    Require that logft be at a max to add into cross
                          correlation.
         -RNEGVAL       Sets alternate values of harmonic comb = this value
                           and sets doubcomb = 1 to use double comb code
         -kPIKPEAK(0)   Will find peak of cross corr fn if pikpeak = 1.
               -mGRPHPEAK(0)   Will print frame with the channel of the"picked
                                   peak" = 1.
               -ePRNTPIK(0)    Will print value of the channel of "picked peak"
               -yMIDIFILE      Will output a file of shorts of the midivalues
                                   corresponding to the "picked peak" for
                                   graphing with ~brown/bin/sndout+
are just using the 1st value from miditable[lam] so pretty stupid to calculate
the whole thing
/***************************************************************************/

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
#define NHARM 6
#define PMODE 0644
#define HANMAX 700000 /*ems can't handle 700000 /* works with 500000 */
#define MINLAM 4                       /* at 10k 6  */
#define MAXLAM 250  /* for c3 ; was 195                      /* at 10k 75 */

struct  frmhdr hdr;
short   dumphdr[1024];                  /* dump header */
short   sampbuf[WINDMAX];               /* buffer for input samples        */
int miditable[MAXLAM + 1];
float   outbuf[MAXCHNLS], *outbufp;     /* frame of out coefs in mag or db */
char    midifile[20];
char    hanfile[20] ;
char    logftvalin[20];
/* double  hanfilrd[HANMAX]; */
float   srate = SRATE, midtuncor=1.;
main(argc,argv)
 int argc;
 char **argv;
{  /* begin main */
        double  *hanfilrd;
        double  theta, a, b, c, db, *wsin, *wcos, *hanwindow, one06, one03, y;
        double  hansum[100], hspacd[21];
        double  onedws, pidws, twopidws, RES2pidws, freqrat, alpha,phase;
        int     frmsz = FRMSZ, windbytes, framebytes, ncoefs,nchnls=NCHNLS;
        int     p, hanning = 1, dbout = 0, print = 0,frmcnt = 0,combprint=0,l;
        int     windsiz[MAXCHNLS], hspac[41], *hspacptr, xcorri=0,graph=0;
        int     flag, windmaxi, sumwind=0, res = RES, nharm=NHARM;
        int     points, m, xcorrprnt=0, nw,rempoints, pikpeak=0,hisamp,sampno;
        int     nextsamp, grphpeak=0, prntpik=0, ofp, hanwritten=0, flag2,ifp;
        int     minmidi, minlam, line=0, noout=0, nrmxcorr=0, halfsemi=1;
        int     midifreq, midi[MAXCHNLS], varres=2, flaghi = 1, flagm = 0;
        int     ofpm, backzero=1, nologftcalc=0, ifp2, fltresp=0, maxreq=0;
        int     hamming = 1, windsiz_1 = 0,siz, calcphase=0, doubcomb=0;
        int     rlcmbg, writehdr=0, dropmidchnls=0;
        short   pitch;
        char    *cp;
        float   windmax, twopi = TWOPI, pi = PI, minhz= MINHZ, tuncorrec=1.00;
        float   neg, *combptr, comb[MAXCHNLS], xcorr[MAXCHNLS],windsizf[MAXCHNLS];
        float   dltn, dltnm1, xn, xnm1, *xcorrp, hival,nextval;
        float   *peakpicp, peakpic[MAXCHNLS], hivaloutbuf, hivalxcorr,normxc;
        float   negval = -.02;
        register int    n, k;
        register short  *samp, *samp2;
        register double *sinp, *cosp, *hanp, *hanfilrdp;

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
                        case 'P': sscanf(cp, "%d", &print);
                                break;
                        case 'f': sscanf(cp, "%f", &minhz);
                                break;
                        case 'F': sscanf(cp, "%d", &frmsz);
                                break;
                        case 'r': sscanf(cp, "%d", &res);/* resolution f/delf*/
                                break;
                        case 'p':sscanf(cp,"%d",&nchnls);/*num points in ft*/
                                break;
                        case 'c': sscanf(cp, "%d", &combprint);
                                break;
                        case 'n': sscanf(cp, "%d", &nharm);
                                break;
                        case 'g': sscanf(cp, "%d", &graph);
                                break;
                        case 'x': sscanf(cp, "%d", &xcorri);
                                break;
                        case 't': sscanf(cp, "%d", &xcorrprnt);
                                break;
                        case 'm': sscanf(cp, "%d", &grphpeak);
                                break;
                        case 'k': sscanf(cp, "%d", &pikpeak);
                                break;
                        case 'e': sscanf(cp, "%d", &prntpik);
                                break;
                        case 'i': sscanf(cp, "%s", hanfile);
                                flag2=1; hanwritten=0;
                                break;
                        case 'j': sscanf(cp, "%d", &hanwritten);
                                break;
                        case 'l': sscanf(cp, "%d", &line);
                                break;
                        case 'o': sscanf(cp, "%d", &noout);
                                break;
                        case 'z': sscanf(cp, "%d", &nrmxcorr);
                                break;
                        case 'a': sscanf(cp, "%d", &halfsemi);
                                break;
                        case 'v': sscanf(cp, "%d", &varres);
                                break;
                        case 'y': sscanf(cp, "%s", midifile);
                                flagm=1;
                                break;
                        case 'q': sscanf(cp, "%d", &backzero);
                                break;
                        case 'w':  sscanf(cp, "%s", logftvalin);
                                nologftcalc=1;
                                break;
                        case 'W':  sscanf(cp, "%d", windsiz_1);
                                break;
                        case 'M': sscanf(cp, "%d", &maxreq);
                                break;
                        case 'T': sscanf(cp, "%f", &tuncorrec);
                                break;
                        case 'S': sscanf(cp, "f", &srate);
                                break;
                        case 's': sscanf(cp, "%d", &calcphase);
                                break;
                        case 'R': sscanf(cp, "%f", &negval);
                                doubcomb = 1;
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
        if(!nologftcalc){
            hanfilrd = (double *)malloc(HANMAX*8);
            if (hanfilrd == NULL)
              die ("memory allocation failure");
        }
        if(hamming) {
                 alpha=(double)25./(double)46.;
                 fprintf(stderr,"Using Hamming window\n");
        }
        else alpha = (double).5; /* hanning */
        makemiditable();
        minlam = (int)((srate/minhz)+ .5);
        minmidi = miditable[minlam];
        if(nchnls > MAXCHNLS)die("Too many channels");

        if(doubcomb) nharm *= 2;
        minhz /= tuncorrec;

        if (flag2){
          ofp = creat(hanfile,PMODE);
          fprintf(stderr,"Prnting hanning values to file %s\n", hanfile);
        }
        if (flagm){
          ofpm = creat(midifile,PMODE);
          fprintf(stderr, "Prnting midi values to file %s\n", midifile);
        }
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
        if(!nologftcalc){
           if (hanwritten){
              if((ifp = open("hanning",O_RDONLY)) < 0)
                die("can't open hanning");
              if((n=read(ifp, hanfilrd, 8*HANMAX))<0)
                die("can't read hanning");
              hanfilrdp = hanfilrd;
           }
        }
        if (nologftcalc){
           if((ifp2 = open(logftvalin,O_RDONLY)) < 0)
             die("Can't open logftvalin");
        }
        one06 = pow(2.,1./12.);
        one03 = pow(2.,1./24.); /* for stepping by 1/2 semitone */
                                     /* freqrat = (float)MAXHZ /minhz; */
                                /*nchnls = (int)(log(freqrat)/log(one06)) + 1*/
        windmax = (float)(res*srate)/minhz; /* res determines harmonic
                    number which is const=res; freq is varied with windsiz[k]
                    and equals  srate*res/windsiz[k] */
        if(!nologftcalc){
           if (n= read(0, dumphdr, 1024) != 1024)
             die("something wrong with header dump\n");
        }
        p = nharm;  l=1;
        while(p--){ /*calc spacing of harmonics hspac[l] in channels
                      for nharm harmonic components */
           if(halfsemi==1)hspacd[l] = log((double)l)/log(one03);
           else hspacd[l] = log((double)l)/log(one06);
           hspac[l] = (int)(hspacd[l] + .5);
           if(combprint)
             fprintf(stderr,
               "nharm %d hspacd[%d]= %.2f\tdiff=%.2f (int)hspac= %d\n"
                ,l,l,hspacd[l],hspacd[l]-hspacd[l-1], hspac[l]);
           ++l;
        }
        hspacptr = hspac; ++hspacptr;  /* harmonic spacing pointer */
       k = nchnls; l = 0;
        if (backzero) neg = 0;
        else
          neg = -(float)nharm/(nchnls - nharm); /* makes ave of comb 0 */
        while (k--){
                   if(doubcomb){
                      m = hspacptr - hspac;
                      comb[l] = (l == *hspacptr ? (m % 2 == 0 ? 1. : negval) : neg);
                   }
                   else comb[l] = (l == *hspacptr ? 1. : neg); /* comb = 1 at
                                    a harmonic position and neg elsewhere*/
                   if ((comb[l]==1. || comb[l]==negval) && l < hspac[nharm]) ++hspacptr;
                   if (combprint)
                       fprintf(stderr, "comb[%d] = %.1f\n", l, comb[l]);
                   ++l;
        }
        k=nchnls;       n=0;
        while (k--){
           midi[n]= minmidi + ((halfsemi==1) ? (n >> 1) : n) ;
           if(halfsemi==1){  /* 2 chnls per semitone */
             windsizf[n] = (float)windmax/pow(one03,(float)n);
           }
           else{
             windsizf[n] = (float)windmax/pow(one06,(float)n);
           }
           if(varres && (midi[n]>90))windsizf[n] = varres * windsizf[n];
                         /* for notes over G6 use varres * resolution */
           windsiz[n] = (int)windsizf[n];
           sumwind += windsiz[n];        /* total window space needed
                                         for sin tables = sum of all windows */
           ++n;
         }
        if(!nologftcalc){
           if(2*sumwind > HANMAX)die("2*Sumwind is > HANMAX.");
           windmaxi = (int)windmax;     /* no of samples to read in is the size
                                           of the largest window (lowest freq*/
           if(windmaxi > WINDMAX)die("Windmaxi larger than WINDMAX");
           windbytes = windmaxi *2;     /*windmaxi is old variable windsiz in dft*/
           framebytes = frmsz * 2;

           if (!hanwritten){            /* Window (or rect) values  not
                       previously written to file so must calculate them */
                 fprintf(stderr,"logft: %s window, %s out, making tables ..\n",
                     (hanning) ? "hanning":"rect", (dbout) ? "db":"magnitude");
                 hanp = hanfilrd;
                 for (k=0; k < nchnls; ++k){
                     if(windsiz_1)twopidws = twopi/(windsiz[k] - 1);
                     else twopidws = twopi/(windsiz[k]);
                     if(varres && (midi[k]>90))
                          RES2pidws = (float)varres* res* twopi/windsiz[k];
                     else RES2pidws = res* twopi/windsiz[k];
                     onedws = 1./windsiz[k];

                     for (n=0; n < windsiz[k]; ++n){
                        a=onedws*((!hanning)?1.:alpha-((1-alpha)*cos(n*twopidws)));
                        theta = n * RES2pidws;
                        *hanp++ = a * sin(theta);
                        *hanp++ = a * cos(theta);
                        /* if want Hanning file written */
/*                      if(flag2){
                          write(ofp, sinp, 8); write(ofp, cosp, 8);
                        } *//* ofp points to hanfile which can be input */

                     }
                 }
          } /* end calc of Hanning window factor */
          if(flag2){
             siz = hanp - hanfilrd;
            write(ofp,hanfilrd,siz*8);
          }
          n = read(0, sampbuf, windbytes);      /* init samp window w. input */
          if (n != windbytes)
                die("premature end of infile");
      }
frame:
        if(!nologftcalc){

           for (k=0,hanfilrdp=hanfilrd; k<nchnls; k++) { /* for one frame: */
              a = 0.0;
              b = 0.0;
              samp =sampbuf;
              if(varres && (midi[k]>90) && flaghi){res=varres*res; flaghi=0;}
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
              }
              else{
                 c = sqrt( a*a + b*b );
                 if (!dbout)
                   outbuf[k] = c;       /* stor as magnitude */
                 else outbuf[k] = db = 20. * log10(c); /*       or db        */
              }
              if (print){
                 midifreq = miditable[(int)(windsiz[k]/res)];
                 if (!dbout)fprintf(stderr, /* print for mag coeff */
                                    "chnl %d midi=%d freq %.1f coef %.1f\n",
                                    k,midifreq, srate*(midi[k]>90?varres*res:res)/windsiz[k], c );

                 else                   /* print for db coeff */
                   fprintf(stderr,"chnl %d midi=%d freq %.1f db=%.1f\n",
                           k, midifreq, srate*(midi[k]>90?varres*res:res)/windsiz[k], db );
              }



           }                            /* end calc from table read in */
        }
       if(nologftcalc){
          if((n=read(ifp2, outbuf, nchnls*4))< nchnls*4)
            exit(-1);
       }
       if(xcorri){ /* cross correlate   */ /* note must change for one06
                                     or not use doubcomb */

          points = nchnls - hspac[nharm];/* do xcorrelation on "points"points*/
                                    /*stop before hspac[] reaches beyond EOF*/
       if(doubcomb){ /* do "for ..." for doubcomb case */
          rlcmbg = (int) (log((double)2.)/log(one03));
          for (k=0; k < points; ++k){   /* begin for of xcorr values */
             xcorr[k] = 0.;             /* initialize cross corr for kth point */
             l = k;
             m = 0;
             if (k < rlcmbg){
                combptr = comb + rlcmbg;
                p = nchnls - k;
                m = k;
                while(p--)
                  xcorr[k] += outbuf[m++] * (*combptr++);
             }
             else {
                combptr = comb;
                /*zzz*/                p = nchnls - k + rlcmbg;
                m = k -rlcmbg;
                while(p--)
                   xcorr[k] += outbuf[m++] * (*combptr++);
             }
          }
       }
      else{
          for (k=0; k < points; ++k){
             xcorr[k] = 0;              /* initialize cross corr for kth point */
             l = k;
             m = 0;
             while(l--)
               xcorr[k] += outbuf[m++] * neg; /* mult points before kth
                                                 by neg*/
             p = nchnls - k;            /* p points from k to nchnls */
             combptr = comb;            /* start comb on the kth point*/
             while(p--){
                if(!maxreq || (*combptr && outbuf[m-1]<outbuf[m] &&
                               outbuf[m]>outbuf[m+1])){ /* if not requiring max
                                     /* or if at a max of logft */
                       xcorr[k] += outbuf[m++] * (*combptr++);

                     }
                     else {  /* max req but no max in fn  or not at the
                                position of a harmonic */
                       xcorr[k] += outbuf[m++] * neg;
                       ++combptr;
                     }
                   }
                 } /* end for calc of points values of xcorr[k] */
          }/* end else for ordinary comb */
          /*both can use from here */
          rempoints = nchnls - points;/*remaining points for number of points
                                         in frame = nchnls*/
          while(rempoints--)
                   xcorr[k++] = 0.;
          /* to print to stdout */ k = nchnls; l=0;
          if(xcorrprnt){
            while(k--)
              fprintf(stderr,"xcorr[%d]=%.2f\n", l++, xcorr[l]);
          }

       }/* end xcorrelation */
if(!nologftcalc){
       if (noout==0){/*if noout==1 don't print to output file; want numbers */
              if((nw = write(1, outbuf, 4*nchnls)) < 0)/*print outbuf to file*/
                fprintf(stderr, "nw = %d  ERROR\n", nw);
       }
       if(line==1)exit(0);  /* Stop if want just one frame */
       if (nw < 4*nchnls && noout ==0){
            fprintf(stderr, "error writing from outbuf\n wrote %d bytes\
                and should have written %d\n", nw, 4*nchnls);
       }
}
        /* if want normalized graph of xcorrelation */
        if(nrmxcorr && xcorri){
                 hivaloutbuf = 0.; xn = 0.;
                 l = nchnls;
                 outbufp = outbuf;
                 while(l--){ /*find max of outbuf for current frame */
                          xnm1 = xn;
                          xn = *outbufp++;
                          dltnm1 = dltn;
                          dltn = xn - xnm1;
                          if(((dltn*dltnm1)<0) && (dltn < 0)){ /* positive pk*/
                                if(xnm1 > hivaloutbuf){
                                         hivaloutbuf = xnm1;
                                }

                          }
                 }/* end of calc max of outbuf */

                 hivalxcorr = 0.; xn = 0.;
                 l = nchnls;
                 xcorrp = xcorr;
                 while(l--){ /*find max of xcorr for current frame */
                          xnm1 = xn;
                          xn = *xcorrp++;
                          dltnm1 = dltn;
                         dltn = xn - xnm1;
                          if(((dltn*dltnm1)<0) && (dltn < 0)){ /* positive pk*/
                                if(xnm1 > hivalxcorr){
                                         hivalxcorr = xnm1;
                                }
                          }
                 }/* end find max for xcorr for this frame */
                 normxc = hivaloutbuf/hivalxcorr;
                 l= nchnls; xcorrp = xcorr;
                 while(l--)/* normalize xcorrp to have same max as outbuf */
                   *xcorrp = normxc* *xcorrp++;
       } /* end normalization of xcorr */

       if(xcorri==1 && graph==1){  /*if want graph of cross correlation */
          if ((nw = write(1, xcorr, 4*nchnls)) < 0){
            fprintf(stderr, "nw = %d  ERROR writing comb\n", nw);exit(1);
          }
          if (nw < 4*nchnls){
            fprintf(stderr, "error writing from comb\n wrote %d bytes\
                and should have written %d\n", nw, 4*nchnls); exit(1);
          }
       }  /*end if want graph of cross correlation */
       ++frmcnt;

        /* peak picker for cross correlation */
        if(pikpeak && xcorri){
                 sampno = 0; hival = 0.; xn = 0.;
                 l = nchnls;
                 xcorrp = xcorr;
                 while(l--){ /*search for peaks in  current frame */
                          xnm1 = xn;
                          xn = *xcorrp++;
                          dltnm1 = dltn;
                          dltn = xn - xnm1;
                          if(((dltn*dltnm1)<0) && (dltn < 0)){ /* positive pk*/
                                if(xnm1 > hival){
                                         nextval = hival;
                                         hival = xnm1;
                                         nextsamp = hisamp;
                                         hisamp = sampno-1;
                                }
                                else if(xnm1 > nextval){
                                         nextval = xnm1;
                                         nextsamp = sampno-1;
                                }
                          }
                          ++sampno;
                 }/* end search for peaks */
                 if(flagm){ /* print midi value picked to file for graphing */
                   pitch = (short)(((halfsemi==1) ? (hisamp >> 1) : hisamp)
                                   + minmidi);
                   if (dropmidchnls)/* throw out channels half way bet tuned
                                       notes */
                     if ((pitch - (short)minmidi) % 2 != 0)
                       pitch = 0;
                   write(ofpm, &pitch,2);
                 }

                 if(prntpik)
                   fprintf(stderr,"\t%d=midip (%d=chnl max %.0f)\t%d=nextmidip (%d=chnl max %.0f)",((halfsemi==1) ? (hisamp >> 1) : hisamp
) +minmidi, hisamp,hival,((halfsemi==1) ? (nextsamp >> 1) : nextsamp) +minmidi, nextsamp,
                           nextval);
                 if(grphpeak){ /* if want graph of result of peakpicker*/
                          l = nchnls;
                          peakpicp = peakpic; n=0;
                          while(l--){
                                   if(n++==hisamp)*peakpicp++ = hival;
                                   else *peakpicp++ = 0. ;
                          }
                          write(1, peakpic, 4*nchnls);
                 }/* end if want graph of result of peakpicker*/
        } /*end peak picker  */

if(!nologftcalc){
       if ((n = windmaxi - frmsz) > 0) {       /* nxt: for step < windowsiz */
                samp = sampbuf;
                samp2 = sampbuf + frmsz;
                while (n--)
                        *samp++ = *samp2++;     /*      slide the sample buf */
                n = read(0, samp, framebytes);  /*     & refill to windmaxi  */
                if (n == framebytes){
                        goto frame;
                }
        }
        else {
                for (n = -n; n>=windmaxi; n-=windmaxi)  /* else waste any   */
                        read(0, sampbuf, windbytes);    /*   extra samps    */
                if (n)  read(0, sampbuf, n*2);
                n = read(0, sampbuf, windbytes); /* then rd whole new window */
                if (n == windbytes){

                        goto frame;             /*      & go calc new frame */
                }
        }
}
        if(nologftcalc)
                goto frame;
} /* end main */

die(s)
 char *s;
{
        fprintf(stderr,"DIE, HELAS!  %s\n",s);
        exit(1);
}

makemiditable() {  /* returns midi value for given period in samples */
         float A, B, freq;
         int lam, freqi;
         A = (float)(12.0/log(2.0));
         B = (float)(69.0 - A* log(440.0));
         for (lam = MINLAM; lam <= MAXLAM; ++lam) {
/*                freq = (srate/lam)*1.04;    /* 4% flat tuning  */
                  freq = (srate/lam)*midtuncor;
                  miditable[lam] = (int)(A * log(freq) + B + 0.5);
                  freqi = (int)freq;
         }
}
