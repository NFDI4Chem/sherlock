#define SIZE_CUR_STRING 5
#define MAX_SIZE_FRAG_NAME 31

/* opcode definition, double ended queue of integer pairs */
typedef struct Opcode {
/* pointer to next opcode */
	struct Opcode *next ;
/* pointer to previous opcode */
	struct Opcode *prev ;
/* pointer to nexp opcode for JZ an JNZ opcodes,
to be used when the test succeeds */
	struct Opcode *jump ;
/* opcode mnemonic, as defined by enum opcodes_mn */
	int mnemo ;
/* address, used for jump, label and fragment opcodes */
	int addr ;
/* value for fragment opcodes, caches fragment search */
	int value ;
/* adress of a value, points to the single value that is
concerned by a given fragment addr(ess) */
	int *pvalue ;
} opcode ;
typedef opcode *popcode ;
#define ONE ((popcode)1)

typedef struct Lvm {
/* the string that constains the logical expression
that will be parsed before sub-structure search.
Il must be initialized before anything is started. */
	char *parsed_string ;
	int parsed_string_index ;

/* value of last parsed token id */
	int cur_token ;

/* stores current string for fragment id and logical operators */
	char cur_string[MAX_SIZE_FRAG_NAME] ;

/* pointers to first and last opcodes */
	popcode root ;
	popcode last ;

/* number of independent fragments */
	int fragc ;

/* the vector of fragment (numeric) names */
	int *fragv ;

#ifdef DEBUG
/* counts calls to get_token() */
	int ngt ;
#endif
} lvm ;
typedef lvm *plvm ;

/* Evaluation function */

typedef int (*lvm_eval)(int) ;

/* function prototypes */

extern plvm lvm_new(char *s) ;
extern int lvm_error(plvm pl) ;
extern void lvm_vars(plvm pl, int *pfragc, int **pfragv) ;
extern int lvm_run(plvm pl, lvm_eval f) ;
extern void lvm_clean(plvm pl) ;
