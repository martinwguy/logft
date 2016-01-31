/* LOGFT.C
 * Reads a soundfile writes a log-frequncy-axis spectrogram.
 *
 * Usage: logft [options] infile.wav outfile.png
 *
 * Adapted by Martin Guy <martinwguy@gmail.com> in January 2016 from
 * Judith Brown's original C code from 1988-90.
 */

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <fcntl.h> /* defines O_RDONLY */
#include "./frmhdr.h"

#include "sndfile.h"
#include "png.h"

/* Default values for parameters */
#define PPSEC 88.2     /* Number of analyses (hence output pixels) per second */
#define MINHZ 130.80 /* c3 midi = 48 */
#define MAXHZ 11442.0
#define PPSEMI 2
#define DYNRANGE 100

#define PI  M_PI
#define TWOPI  (2.0 * PI)

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
        double  *hanfilrd;
        double  onedws, twopidws, RES2pidws, alpha;
	double	delta_f_over_f;
	float	ppsec=PPSEC;	/* Output columns per second */
        int     frmsz, nchnls;
        int     hann = 1, frmcnt = 0;
        int     *windsiz;	/* [nchnls] */
        float   *windsizf;	/* [nchnls] */
        int     flag, windmaxi, sumwind=0;
	float	res = 0.;	/* f/deltaf = 1/.06 */
	int	ppsemi=PPSEMI;	/* Number of buckets per semitone */
        int     hamming = 1;
        char    *cp;
        float   windmax, minhz=MINHZ, maxhz=MAXHZ;
        register int    n, k, y;
        register double *hanfilrdp;
	float  maxAmp=0.0, dynRange=(float)DYNRANGE;

	/* Stuff for reading sound file */
	char *fileNameIn;
	SNDFILE *sndfile;
	SF_INFO sfinfo;

	/* Stuff for writing PNG file */
	char *fileNameOut;
	png_structp png_ptr;
	png_infop png_info;
	png_uint_32 png_width, png_height;

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
		case 'x': sscanf(cp, "%f", &maxhz);
			break;
		case 'p': sscanf(cp, "%d", &ppsemi);
			break;
		case 's': sscanf(cp, "%f", &ppsec);
			break;
		case 'A': sscanf(cp, "%f", &maxAmp);/* value to display as white */
			break;
		case 'D': sscanf(cp, "%f", &dynRange);	/* maxAmp-dynRange displays as black */
			break;
		case 'S': sscanf(cp, "%i", &srate);
			break;

		default: fprintf(stderr, "unknown option \"%c\"\n", flag);

usage:
fprintf(stderr, "usage: logft [options] infile.wav outfile.png\n\
-h0           Use rectangular window instead of Hann/Hamming\n\
-H1           Use Hamming window instead of Hann\n\
-fMINHZ       Lowest frequency value (%g Hz)\n\
-xMAXHZ       Highest frequency value (%g Hz)\n\
-pPPSEMI      Number of frequency values per semitone (%d)\n\
-sPPSEC       Number of output columns per second (%g)\n\
-SSRATE       Sample rate (default: auto-detected)\n\
-AmaxAmp(0)   Output value to display as white (negative values brighten output)\n\
-DdynRange(%d) Dynamic range of output.\n\
", MINHZ, MAXHZ, PPSEMI, PPSEC, DYNRANGE);
			exit(1);
			break;
                }
                argc--; argv++;
        }

	if (argc != 2) goto usage;
	fileNameIn  = argv[0];
	fileNameOut = argv[1];

	if (hann) {
            if(hamming) alpha=(double)25./(double)46.;
            else alpha = (double).5; /* hann */
	}

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

	frmsz = (int)(srate / ppsec + 0.5);

	delta_f_over_f = pow(2.,1./(12.*ppsemi));
	if (res == 0.0) res = 1.0 / (delta_f_over_f - 1.0);

        nchnls = (int)(log(maxhz/minhz)/log(delta_f_over_f)) + 1;
        windmax = (float)(res*srate)/minhz; /* res determines harmonic
                    number which is const=res; freq is varied with windsiz[k]
                    and equals  srate*res/windsiz[k] */

	windsiz = calloc(nchnls, sizeof(*windsiz));
	windsizf = calloc(nchnls, sizeof(*windsizf));
	outbuf = calloc(nchnls, sizeof(*outbuf));
	if (!windsiz || !windsizf || !outbuf)
	    die("Not enough memory");

	/* Calculate window size for each bucket */
	for (k=0; k < nchnls; ++k) {
           windsizf[k] = (float)windmax/pow(delta_f_over_f,(double)k);
           windsiz[k] = (int)windsizf[k];
           sumwind += windsiz[k];        /* total window space needed
                                         for sin tables = sum of all windows */
        }

	/* Allocate memory for sin/cos pairs */
        hanfilrd = (double *)malloc(2*sumwind*sizeof(double));
        if (hanfilrd == NULL)
          die ("memory allocation failure");
		
	windmaxi = (int)windmax;     /* no of samples to read in is the size
					of the largest window (lowest freq)*/

	sampbuf = calloc(windmaxi, sizeof(*sampbuf));
	if (!sampbuf) die("Not enough memory");

	/* Calculate window (or rect) values */
	hanfilrdp = hanfilrd;
	for (k=0; k < nchnls; ++k){
	    twopidws = TWOPI/(windsiz[k]);
	    RES2pidws = res* TWOPI/windsiz[k];
	    onedws = 1./windsiz[k];

	    for (n=0; n < windsiz[k]; ++n){
	       double a = onedws*((!hann)?1.:alpha-((1-alpha)*cos(n*twopidws)));
	       double theta = n * RES2pidws;
	       *hanfilrdp++ = a * sin(theta);
	       *hanfilrdp++ = a * cos(theta);
	    }
	}

	/* Create PNG file for writing */
	{
	    FILE *png_fp = fopen(fileNameOut, "wb");

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
	    png_height = ppsec * sfinfo.frames / srate;
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

	/* The analysis is centred on the centre of the window, so to get an
	 * analysis of the whole audio file, we start with half a window of
	 * silence and half a window of audio */
	memset(sampbuf, 0, windmaxi/2);
	n = sf_readf_float(sndfile, sampbuf+windmaxi/2, windmaxi-windmaxi/2);
        if (n < windmaxi-windmaxi/2) {
	    /* Fill the remaining space with silence */
            memset(sampbuf+windmaxi/2+n, 0, windmaxi-windmaxi/2-n);
	}

	for (y=0; y < png_height; y++) {
	    hanfilrdp = hanfilrd;
	    for (k=0; k<nchnls; k++) { /* for one frame: */
	       double a = 0.0, b = 0.0, c;
	       /* Center all buckets on the central sample of the window
		* to avoid skewing the output */
	       float *samp = sampbuf + (windmaxi - windsiz[k])/2;
	       for (n=0; n<windsiz[k]; n++, samp++) { /*  calculate coefs  */
		  a += *samp * *hanfilrdp++; /* sin */
		  b += *samp * *hanfilrdp++; /* cos */
	       }
	       c = sqrt( a*a + b*b );
	       outbuf[k] = 20. * log10(c);
	    }                            /* end calc from table read in */

	    output_frame(png_ptr, outbuf, png_width, maxAmp, -dynRange);

	    ++frmcnt;

	    if ((n = windmaxi - frmsz) > 0) {  /* nxt: for step < windowsiz */
		float *samp = sampbuf;
		float *samp2 = sampbuf + frmsz;
		while (n--)
			*samp++ = *samp2++;    /* slide the sample buf */
		n = sf_readf_float(sndfile, samp, frmsz);  /* & refill to windmaxi  */
		/* At EOF, fill with silence */
		if (n < frmsz) memset(samp+n, 0, frmsz-n);
	    } else {
						/* else waste any extra samps */
		for (n = -n; n>=windmaxi; n-=windmaxi)
			sf_readf_float(sndfile, sampbuf, windmaxi);
		if (n)  sf_readf_float(sndfile, sampbuf, n);
						/* then rd whole new window */
		n = sf_readf_float(sndfile, sampbuf, windmaxi);
		/* At EOF, fill with silence */
		if (n < windmaxi) memset(sampbuf+n, 0, windmaxi-n);
	    }
	}

	/* Finish writing the PNG file */
	png_write_end(png_ptr, NULL);

	show_min_max();

	exit(0);

} /* end main */

/* Accumulate the minimum and maximum output values so that
 * they can normalize the brghtness if they want */
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
