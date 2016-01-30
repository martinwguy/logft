/* frmhdr.h */
/* Data header for frame data */

/* In general, zero for any field indicates the value is unknown */

union maxmin {
   short s;
   int i;
   float f;
};

struct frmhdr {
   int frm_magic;                       /* 440 decimal */

   int frm_width;       /*nchnls*/      /* F-domain values per frame. */
   int frm_count;                       /* Number of frames in file,
                                           or zero if unknown. */
   int frm_datatype;    /* int,float..  /* The datatype of the data. */
   int frm_scale;       /* lin(1)or log(2)/* The scaling of the data,
                                           e.g. linear or log. */
   int frm_frm_size;    /* frmsz */     /* T-domain samples per frame. */

/**/float frm_sr;       /*srate*/       /* Original T-domain sample rate. */
/**/float frm_frm_p_s;  /* sr/frmsz */  /* Frames per second. */
   float frm_minfreq;    /* minhz*/     /* Lowest freq calc   */

   union maxmin max;                    /* Maximum data value contained. */
   union maxmin min;                    /* Minimum ditto. */

   int spare[19];
};

#define FRM_MAGIC (440)

#define FRM_DATA_SHORT  (1)
#define FRM_DATA_INT    (2)
#define FRM_DATA_FLOAT  (3)

#define FRM_SCALE_LINEAR        (1)
#define FRM_SCALE_LOG           (2)

/**/ /* NB these are ints in some progs and floats in others */
