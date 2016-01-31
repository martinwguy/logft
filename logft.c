/* LOGFT.C
 * Takes as infile a soundfile using libsndfile
 *
 * adaptation of dft.c using log channels separated by
 * a semitone with resolution of a semitone
 * later modified; two chnls/semitone with variable resolution
 * calculates discrete Fourier transform of sound samples
 * and writes a PNG file to the named output file.
 */

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <fcntl.h> /* defines O_RDONLY */
#include "./frmhdr.h"

#include "sndfile.h"
#include "png.h"

#define FRMSZ 500
/*#define MINHZ 174.6  /* F3 = 7.05  midi=53 */
#define MINHZ 130.80 /* c3 midi = 48 */
/*#define MAXHZ 1600  */
#define RES 34 /* f/deltaf = 1/.06 */
#define NCHNLS 156
#define DYNRANGE 100

#define TWOPI  6.2831854
#define PI  3.1415927
#define PMODE 0644
#define MINLAM 4                       /* at 10k 6  */
#define MAXLAM 250  /* for c3 ; was 195                      /* at 10k 75 */

/* Forward function declarations */
static void output_frame(png_structp png_ptr, float *outbuf,
			 png_uint_32 png_width, float maxAmp, float dynRange);
static void show_min_max(void);

main(argc,argv)
 int argc;
 char **argv;
{  /* begin main */
	float   *sampbuf;        /* buffer for input samples        */
	float   *outbuf;         /* frame of out coefs in mag or db [nchnls] */
	int     srate = 0;
	float   midtuncor=1.;
        double  *hanfilrd;
        double  theta, a, b, c, db, one06, one03;
        double  onedws, twopidws, RES2pidws, alpha;
        int     frmsz = FRMSZ, nchnls=NCHNLS;
        int     hann = 1, frmcnt = 0;
        int     *windsiz;	/* [nchnls] */
        float   *windsizf;	/* [nchnls] */
        int     flag, windmaxi, sumwind=0, res = RES;
        int     halfsemi=1;
        int     maxreq=0;
        int     hamming = 1, windsiz_1 = 0;
        char    *cp;
        float   windmax, minhz=MINHZ, tuncorrec=1.00;
        register int    n, k;
        register float  *samp, *samp2;
        register double *hanp, *hanfilrdp;
	float  maxAmp=0.0, dynRange=(float)DYNRANGE;

	/* Stuff for reading sound file */
	char *fileNameIn;
	SNDFILE *sndfile;
	SF_INFO sfinfo;

	/* Stuff for writing PNG file */
	char *fileNameOut;
	png_structp png_ptr;
	png_infop png_info;
	png_uint_32 png_width;

        argc--; argv++;                         /* align onto command args */
        while ((cp = *argv) && *cp++ =='-' && (flag = *cp++))
        {
                switch (flag){
		case 'h': sscanf(cp, "%d", &hann);
			break;
		case 'H': sscanf(cp, "%d", &hamming);
			break;
		case 'f': sscanf(cp, "%f", &minhz);
			break;
		case 'F': sscanf(cp, "%d", &frmsz);
			break;
		case 'r': sscanf(cp, "%d", &res);/* resolution f/delf*/
			break;
		case 'p':sscanf(cp,"%d", &nchnls);/*num points in ft*/
			break;
		case 'a': sscanf(cp, "%d", &halfsemi);
			break;
		case 'A': sscanf(cp, "%f", &maxAmp);	/* value to display as white (usually negative) */
			break;
		case 'D': sscanf(cp, "%f", &dynRange);	/* maxAmp-dynRange displays as black */
			break;
		case 'W':  sscanf(cp, "%d", &windsiz_1);
			break;
		case 'M': sscanf(cp, "%d", &maxreq);
			break;
		case 'T': sscanf(cp, "%f", &tuncorrec);
			break;
		case 'S': sscanf(cp, "%i", &srate);
			break;
		case 'Q': sscanf(cp, "%f", &midtuncor);
			break;

		default: fprintf(stderr, "unknown option \"%c\"\n", flag);

usage:
fprintf(stderr, "usage: logft [options] infile.wav outfile.png\n\
-h0           Use rectangular window instead of Hann/Hamming\n\
-H1           Use Hamming window instead of Hann\n\
-rRES(17)     Resolution f/deltaf.\n\
-FFRMSZ(500)  Number of samples analyzed per frame.\n\
-pNCHNLS(156) Number of frequencies for which logft calculated per frame\n\
-fMINHZ(174.6)Lowest frequency value = 7.05 = F3 = midi 53\n\
-aHALFSEMI(1) Calculate with 2 frequency values per semitone.\n\
-TTUNCORREC   Follow with tuning correction for sounds digitized at\n\
              ems(1.04) implemented 7/20/88\n\
-WWINDSIZ-1(0)Use windsiz - 1 in hann(mm) window\n\
-SSRATE       Sample rate (default: auto-detected)\n\
-AmaxAmp(0)   Output value to display as white (negative values brighten output)\n\
-DdynRange(%d) Dynamic range of output.\n\
", DYNRANGE);
			exit(1);
			break;
                }
                argc--; argv++;
        }

	if (argc != 2) goto usage;
	fileNameIn  = argv[0];
	fileNameOut = argv[1];

	if (!hann)
            fprintf(stderr,"Using rectangular window\n");
	else {
            if(hamming) {
                alpha=(double)25./(double)46.;
                fprintf(stderr,"Using Hamming window\n");
            } else {
	        alpha = (double).5; /* hann */
                fprintf(stderr,"Using Hann window\n");
	    }
	}

	windsiz = calloc(nchnls, sizeof(*windsiz));
	windsizf = calloc(nchnls, sizeof(*windsizf));
	outbuf = calloc(nchnls, sizeof(*outbuf));
	if (!windsiz || !windsizf || !outbuf)
	    die("Not enough memory");

        minhz /= tuncorrec;

	/* Open sound file for reading */
	{
	    memset(&sfinfo, 0, sizeof(SF_INFO));
	    sndfile = sf_open(fileNameIn, SFM_READ, &sfinfo);
	    if (!sndfile) {
		    fprintf(stderr, "Failed to open input file \"%s\": %s",
			    fileNameIn, sf_strerror(sndfile));
		    exit(1);
	    }
	    if (sfinfo.channels != 1) {
		    fprintf(stderr, "I can only process mono files as yet\n.");
		    exit(1);
	    }
	    if (srate == 0) srate = sfinfo.samplerate;
	}

	/* Prepare the calculation subsystem */

        one06 = pow(2.,1./12.);
        one03 = pow(2.,1./24.); /* for stepping by 1/2 semitone */
                                     /* freqrat = (float)MAXHZ /minhz; */
                                /*nchnls = (int)(log(freqrat)/log(one06)) + 1*/
        windmax = (float)(res*srate)/minhz; /* res determines harmonic
                    number which is const=res; freq is varied with windsiz[k]
                    and equals  srate*res/windsiz[k] */

        k=nchnls;       n=0;
fputs("Window sizes:", stderr);
        while (k--){
           if(halfsemi==1){  /* 2 chnls per semitone */
             windsizf[n] = (float)windmax/pow(one03,(float)n);
           } else {
             windsizf[n] = (float)windmax/pow(one06,(float)n);
           }
           windsiz[n] = (int)windsizf[n];
           sumwind += windsiz[n];        /* total window space needed
                                         for sin tables = sum of all windows */
fprintf(stderr, " %d", windsiz[n]);
           ++n;
        }
fprintf(stderr, "\n");

        hanfilrd = (double *)malloc(2*sumwind*sizeof(double));
        if (hanfilrd == NULL)
          die ("memory allocation failure");
		
	windmaxi = (int)windmax;     /* no of samples to read in is the size
					of the largest window (lowest freq*/

	sampbuf = calloc(windmaxi, sizeof(*sampbuf));
	if (!sampbuf) die("Not enough memory");

	/* Calculate window (or rect) values */
	hanp = hanfilrd;
	for (k=0; k < nchnls; ++k){
	    if(windsiz_1)twopidws = TWOPI/(windsiz[k] - 1);
	    else twopidws = TWOPI/(windsiz[k]);
	    RES2pidws = res* TWOPI/windsiz[k];
	    onedws = 1./windsiz[k];

	    for (n=0; n < windsiz[k]; ++n){
	       a=onedws*((!hann)?1.:alpha-((1-alpha)*cos(n*twopidws)));
	       theta = n * RES2pidws;
	       *hanp++ = a * sin(theta);
	       *hanp++ = a * cos(theta);
	    }
	}

	/* Read a full window of samples from the audio file */
	n = sf_readf_float(sndfile, sampbuf, windmaxi);
        if (n != windmaxi)
              die("premature end of infile");

	/* Create PNG file for writing */
	{
	    FILE *png_fp = fopen(fileNameOut, "wb");
	    png_uint_32 png_height;

	    if (!png_fp) {
		fputs("Cannot create ", stderr);
		perror(fileNameOut);
		exit(1);
	    }
	    png_ptr = png_create_write_struct(
		PNG_LIBPNG_VER_STRING, NULL, NULL, NULL);
	    if (!png_ptr) {
		fputs("Cannot create PNG writer.\n", stderr);
		exit(1);
	    }
	    png_info = png_create_info_struct(png_ptr);
	    if (!png_info) {
		fputs("Cannot create PNG info.\n", stderr);
		exit(1);
	    }
	    if (setjmp(png_jmpbuf(png_ptr)))
	    {
		fputs("Something went wrong in the PNG-writing subsystem. Quitting...\n", stderr);
		png_destroy_write_struct(&png_ptr, &png_info);
		fclose(png_fp);
		exit(1);
	    }
	    png_width = (png_uint_32) nchnls;
	    png_height = ((sfinfo.frames - windmaxi) / frmsz) + 1;
	    png_init_io(png_ptr, png_fp);
	    png_set_IHDR(png_ptr, png_info,
		png_width, png_height,
		8, // bit depth
		PNG_COLOR_TYPE_GRAY,
		PNG_INTERLACE_NONE,
		PNG_COMPRESSION_TYPE_DEFAULT,
		PNG_FILTER_TYPE_DEFAULT);
	    png_write_info(png_ptr, png_info);
	}

frame:
	for (k=0,hanfilrdp=hanfilrd; k<nchnls; k++) { /* for one frame: */
	   a = 0.0;
	   b = 0.0;
	   samp =sampbuf;
	   for (n=0; n<windsiz[k]; n++) { /*  calculate coefs  */
	      a += *samp * *hanfilrdp++;
	      b += *samp++ * *hanfilrdp++; /*wrote sin then cos*/
	   }
	   c = sqrt( a*a + b*b );
	   outbuf[k] = db = 20. * log10(c);
	}                            /* end calc from table read in */

	output_frame(png_ptr, outbuf, png_width, maxAmp, -dynRange);

        ++frmcnt;

        if ((n = windmaxi - frmsz) > 0) {       /* nxt: for step < windowsiz */
                samp = sampbuf;
                samp2 = sampbuf + frmsz;
                while (n--)
                        *samp++ = *samp2++;     /*      slide the sample buf */
                n = sf_readf_float(sndfile, samp, frmsz);  /*     & refill to windmaxi  */
                if (n == frmsz) goto frame;
        } else {
		                                /* else waste any extra samps */
                for (n = -n; n>=windmaxi; n-=windmaxi)
                        sf_readf_float(sndfile, sampbuf, windmaxi);
                if (n)  sf_readf_float(sndfile, sampbuf, n);
		                                /* then rd whole new window */
                n = sf_readf_float(sndfile, sampbuf, windmaxi);
                if (n == windmaxi) goto frame;  /* & go calc new frame */
        }

	/* Finish writing the PNG file */
	png_write_end(png_ptr, NULL);

	show_min_max();

	exit(0);

} /* end main */

static float minval, maxval;
static int initval=0;	/* Have we initialized minval and maxval? */

static void
show_min_max()
{
    printf("Minimum output value is %f\n", minval);
    printf("Maximum output value is %f\n", maxval);
}

static void
output_frame(png_structp png_ptr, float *outbuf, png_uint_32 png_width,
	     float maxAmp, float dynRange)
{
    static png_bytep png_row = NULL;
    unsigned i;

    if (png_row == NULL) png_row = (png_bytep) malloc(png_width);
    if (!png_row) die("Not enough memory");

    if (!initval) {
	minval = maxval = outbuf[0];
	initval = 1;
    }

    /* Convert floats to 0-255 */
    for (i = 0; i < png_width; i++ ) {
        double value = outbuf[i];
	if (value < minval) minval=value;
        if (value > maxval) maxval=value;
	value -= maxAmp;
	value /= dynRange; /* values are <=0 and dynRange is negative */
        if (value < 0.0) { value = 0.0; }
        if (value > 1.0) { value = 1.0; }
        png_row[i] = ((1.0-value) * 255) + 0.5;
    }
    png_write_row(png_ptr, png_row);
}

die(s)
 char *s;
{
        fprintf(stderr,"DIE, HELAS!  %s\n",s);
        exit(1);
}
